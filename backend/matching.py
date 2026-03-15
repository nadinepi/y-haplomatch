import sqlite3
import logging
import re


SNP24_POSITION_RE = re.compile(r'^snp_24_(\d+)(?:_|$)', re.IGNORECASE)

def _strip_plink_suffix(snp_name):
    # Normalize SNP names for matching.
    if snp_name.lower().startswith("snp_24_"):
        pos = _extract_marker_position(snp_name)
        if pos is not None:
            return f"snp_24_{pos}"

    # remove _1 .1 etc
    if "_" in snp_name:
        snp_name = snp_name.rsplit("_", 1)[0]
    if "." in snp_name:
        snp_name = snp_name.split(".", 1)[0]

    return snp_name


def _extract_marker_position(snp_name):
    # extract Y-position from markers
    m = SNP24_POSITION_RE.match(snp_name)
    if m:
        return int(m.group(1))
    return None


def resolve_haplogroup_simple(db_path, user_input):
    # figure out haplogroup from user input

    h = (user_input or '').strip().rstrip('?*~')
    if not h:
        return None, None, False

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # try using SNP to select haplogroup from db
    candidate_snp = h.split('-', 1)[1] if '-' in h else h
    candidate_snp = candidate_snp.strip()
    for snp_name in (candidate_snp, f'{candidate_snp}.1', f'{candidate_snp}.2', f'{candidate_snp}.3'):
        cur.execute(
            "SELECT haplogroup FROM snp_reference WHERE UPPER(snp_name) = UPPER(?) LIMIT 1",
            (snp_name,),
        )
        row = cur.fetchone()
        if row and row[0]:
            conn.close()
            return row[0], candidate_snp, False

    # otherwise treat as haplogroup text
    token = h.split('-', 1)[0].strip()
    cur.execute("SELECT 1 FROM snp_reference WHERE haplogroup = ? LIMIT 1", (token,))
    has_ref = cur.fetchone() is not None
    if not has_ref:
        cur.execute("SELECT 1 FROM individuals WHERE y_haplogroup_clean = ? LIMIT 1", (token,))
        has_ref = cur.fetchone() is not None

    conn.close()
    if has_ref:
        return token, None, '-' in h
    else:
        return None, None, False
    return None, None, False


def match_by_genotype(db_path, user_haplo, user_genotypes, limit=1000, min_snps=3):
    # normalize user SNP names
    normalized_user_genotypes = {}
    for snp, val in user_genotypes.items():
        base = _strip_plink_suffix(snp)
        normalized_user_genotypes[base.upper()] = val
    user_genotypes = normalized_user_genotypes

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    # get all SNPs on this haplogroup branch + descendants
    cur.execute(
        """
        SELECT DISTINCT snp_name, position
        FROM snp_reference
        WHERE haplogroup = ? OR haplogroup LIKE ?
        """,
        (user_haplo, f"{user_haplo}%")
    )
    relevant_rows = [row for row in cur]
    if not relevant_rows:
        conn.close()
        return []

    relevant_positions = {}
    rsid_by_position = {}
    user_value_by_position = {}
    user_label_by_position = {}
    for row in relevant_rows:
        snp_name = row[0]
        position = row[1]
        if not snp_name:
            continue
        key = snp_name.upper()
        if position is not None and position not in relevant_positions:
            relevant_positions[position] = snp_name
        if position is not None and key.startswith('RS') and position not in rsid_by_position:
            rsid_by_position[position] = key
        if position is not None and key in user_genotypes:
            # keep the strongest user signal for this position (2 over 0).
            prev = user_value_by_position.get(position)
            cur_val = user_genotypes[key]
            if prev is None or cur_val > prev:
                user_value_by_position[position] = cur_val
                user_label_by_position[position] = key

    def _shared_mutation_label(pos):
        canonical = str(relevant_positions.get(pos, pos)).upper()
        rsid = rsid_by_position.get(pos)
        user_label = user_label_by_position.get(pos)

        # prefer what the user has when it is an rsID
        if user_label and user_label.startswith('RS'):
            return user_label if user_label == canonical else f"{user_label} ({canonical})"
        # if not, prefer a known rsID alias for display
        if rsid:
            return rsid if rsid == canonical else f"{rsid} ({canonical})"
        # fall back to user label or canonical ISOGG name
        if user_label:
            return user_label if user_label == canonical else f"{user_label} ({canonical})"
        return canonical

    # load metadata only for relevant haplogroup branch
    cur.execute("""
        SELECT id, y_haplogroup, y_haplogroup_clean, y_terminal,
               group_id, locality, country, lat, lon,
               date_mean, full_date, full_date_range, n_called
        FROM individuals
        WHERE y_haplogroup_clean = ?
           OR y_haplogroup_clean LIKE ?
    """, (user_haplo, f"{user_haplo}%"))
    ind_meta = {row['id']: dict(row) for row in cur}
    if not ind_meta:
        conn.close()
        return []

    # position-first targets, faster than snp name
    target_positions = sorted(user_value_by_position.keys())
    if not target_positions:
        conn.close()
        return []

    # identify user-derived SNPs
    user_derived_positions = {pos for pos, val in user_value_by_position.items() if val == 2}
    user_derived_total = len(user_derived_positions)
    if user_derived_total == 0:
        conn.close()
        return []

    # find best matches by position first
    placeholders = ','.join('?' for _ in target_positions)
    cur.execute(f"""
        SELECT g.individual_id, g.position, MAX(g.is_derived) AS is_derived
        FROM genotypes g
        JOIN individuals i ON i.id = g.individual_id
        WHERE (i.y_haplogroup_clean = ? OR i.y_haplogroup_clean LIKE ?)
          AND g.position IN ({placeholders})
        GROUP BY g.individual_id, g.position
    """, [user_haplo, f"{user_haplo}%", *target_positions])

    ind_genos = {}
    for row in cur:
        iid = row['individual_id']
        ind_genos.setdefault(iid, {})[row['position']] = row['is_derived']
    conn.close()


    # score each individual
    results = []
    for iid, genos in ind_genos.items():
        meta = ind_meta.get(iid)
        if not meta:
            continue


        snps_compared = 0
        snps_matched = 0
        shared_mutations = []

        for pos in target_positions:
            user_val = user_value_by_position[pos]
            sample_val = genos.get(pos)
            if sample_val not in (0, 1):
                continue

            snps_compared += 1
            sample_is_alt = sample_val == 1
            user_is_alt = user_val > 0

            # count only shared derived
            if pos in user_derived_positions and user_is_alt and sample_is_alt:
                snps_matched += 1
                shared_mutations.append(_shared_mutation_label(pos))

        # skip tiny overlaps (less than 3)
        if snps_compared < min_snps:
            continue

        match_score = snps_matched / max(user_derived_total,1)
        derived_agreement = snps_matched / snps_compared if snps_compared > 0 else 0
        low_overlap = snps_compared < 5

        results.append({
            **meta,
            'match_score': round(match_score,4),
            'derived_agreement': round(derived_agreement,4),
            'snps_compared': snps_compared,
            'snps_matched': snps_matched,
            'user_derived_total': user_derived_total,
            'low_overlap': low_overlap,
            'shared_mutations': shared_mutations
        })

    # best match_score, then most overlap
    results.sort(key=lambda r: (-r['match_score'], -r['snps_compared']))
    return results[:limit]


def match_with_broadening(db_path, user_haplo, user_genotypes, limit=1000, min_snps=3):
    # try matching with broader and broader haplogroup
    haplo = user_haplo
    for attempt in range(20):
        print("Trying haplogroup:", haplo)
        results = match_by_genotype(db_path, haplo, user_genotypes, limit=limit, min_snps=min_snps)
        if results or len(haplo) <= 1:
            return results, haplo, attempt
        haplo = haplo[:-1]  # shorten haplogroup
    return results, haplo, 20