import sqlite3
import logging

## Logging to matching_debug.log disabled
def _strip_plink_suffix(snp_name):
    """Normalize SNP names for matching."""
    if snp_name.startswith("snp_24_"):
        return "_".join(snp_name.split("_")[:-1])

    # remove _1 _2 etc
    if "_" in snp_name:
        snp_name = snp_name.rsplit("_", 1)[0]

    # remove .1 .2 etc
    if "." in snp_name:
        snp_name = snp_name.split(".", 1)[0]

    return snp_name


def resolve_haplogroup_simple(db_path, user_input):
    """Get haplogroup from user input (text or SNP name)."""

    h = (user_input or '').strip().rstrip('?*~')
    if not h:
        return None, None, False

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # Try SNP style first
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

    # Otherwise treat as haplogroup text
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
    # Normalize user SNP names
    normalized_user_genotypes = {}
    for snp, val in user_genotypes.items():
        base = _strip_plink_suffix(snp)
        normalized_user_genotypes[base] = val
    user_genotypes = normalized_user_genotypes

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    # get all SNPs on this haplogroup branch + descendants
    cur.execute(
        "SELECT DISTINCT snp_name FROM snp_reference WHERE haplogroup = ? OR haplogroup LIKE ?",
        (user_haplo, f"{user_haplo}%")
    )
    relevant_snps = {row[0] for row in cur}
    if not relevant_snps:
        conn.close()
        return []

    # map base SNP -> PLINK variants
    base_to_plinks = {}
    cur.execute("SELECT DISTINCT snp_name FROM genotypes")
    for (plink_name,) in cur:
        base = _strip_plink_suffix(plink_name)
        base_to_plinks.setdefault(base, []).append(plink_name)

    # keep only SNPs user has typed and exist in genotypes
    target_plink_names = set()
    base_targets = {}
    for base in relevant_snps:
        if base not in user_genotypes:
            continue
        plinks = base_to_plinks.get(base, [])
        if not plinks:
            continue
        base_targets[base] = plinks
        target_plink_names.update(plinks)

    if not base_targets:
        conn.close()
        return []

    # identify user-derived SNPs
    user_branch_bases = set(base_targets.keys())
    user_derived_bases = {b for b in user_branch_bases if user_genotypes[b] == 2}
    user_derived_total = len(user_derived_bases)

    # fetch genotype data for all target SNPs
    placeholders = ','.join('?' for _ in target_plink_names)
    cur.execute(f"""
        SELECT g.individual_id, g.snp_name, g.value
        FROM genotypes g
        WHERE g.snp_name IN ({placeholders})
    """, list(target_plink_names))

    ind_genos = {}
    for row in cur:
        iid = row['individual_id']
        ind_genos.setdefault(iid, {})[row['snp_name']] = row['value']

    # load metadata only for relevant haplogroup branch
    cur.execute("""
        SELECT id, y_haplogroup, y_haplogroup_clean, y_terminal,
               group_id, locality, country, lat, lon,
               date_mean, full_date, n_called
        FROM individuals
        WHERE y_haplogroup_clean = ?
           OR y_haplogroup_clean LIKE ?
    """, (user_haplo, f"{user_haplo}%"))
    ind_meta = {row['id']: dict(row) for row in cur}
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

        for base, plinks in base_targets.items():
            user_val = user_genotypes[base]
            sample_calls = [genos.get(p) for p in plinks if genos.get(p) in (0,1,2)]
            if base == "FGC2806":
                continue
            if not sample_calls:
                continue

            snps_compared += 1
            sample_is_alt = any(v > 0 for v in sample_calls)
            user_is_alt = user_val > 0

            # Count only shared derived
            if base in user_derived_bases and user_is_alt and sample_is_alt:
                snps_matched += 1
                shared_mutations.append(base)


        # Skip tiny overlaps
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

    # 8. Sort: best match_score, then most overlap
    results.sort(key=lambda r: (-r['match_score'], -r['snps_compared']))
    return results[:limit]


def match_with_broadening(db_path, user_haplo, user_genotypes, limit=1000, min_snps=3):
    """Try matching user haplogroup; if nothing found, shorten haplogroup step by step."""
    haplo = user_haplo
    for attempt in range(20):
        print("Trying haplogroup:", haplo)
        results = match_by_genotype(db_path, haplo, user_genotypes, limit=limit, min_snps=min_snps)
        if results or len(haplo) <= 1:
            return results, haplo, attempt
        haplo = haplo[:-1]  # shorten haplogroup
    return results, haplo, 20