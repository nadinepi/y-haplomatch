import logging
import re
import sqlite3


SNP24_POSITION_RE = re.compile(r'^snp_24_(\d+)(?:_|$)', re.IGNORECASE)


def _strip_plink_suffix(snp_name):
    # normalize SNP names for matching
    if snp_name.lower().startswith('snp_24_'):
        pos = _extract_marker_position(snp_name)
        if pos is not None:
            return f'snp_24_{pos}'

    if '_' in snp_name:
        snp_name = snp_name.rsplit('_', 1)[0]

    return snp_name


def _split_aliases(text):
    # some tree nodes use slash-separated aliases
    return [part.strip() for part in str(text).split('/') if part.strip() and part.strip() != '#']


def _extract_marker_position(snp_name):
    # pull the Y position from PLINK marker names
    m = SNP24_POSITION_RE.match(snp_name)
    if m:
        return int(m.group(1))
    return None


def _load_tree_data(cur):
    # load the 2016 tree and its aliases from the database
    cur.execute("SELECT child, parent FROM haplogroup_tree")
    tree_rows = cur.fetchall()

    parent_by_node = {}
    children_by_node = {}
    alias_to_node = {}

    for child_raw, parent_raw in tree_rows:
        child_aliases = _split_aliases(child_raw)
        parent_aliases = _split_aliases(parent_raw)
        if not child_aliases:
            continue

        child = child_aliases[0]
        parent = parent_aliases[0] if parent_aliases else None

        parent_by_node[child] = parent
        children_by_node.setdefault(child, set())
        if parent:
            children_by_node.setdefault(parent, set()).add(child)

        for alias in child_aliases:
            alias_to_node[alias.upper()] = child
        for alias in parent_aliases:
            alias_to_node[alias.upper()] = parent

    return alias_to_node, parent_by_node, children_by_node


def _resolve_tree_name(name, alias_to_node):
    # map any tree alias to its canonical 2016 node name
    if not name:
        return None
    return alias_to_node.get(name.strip().upper())


def _collect_descendants(start_node, children_by_node):
    # collect one node and all of its children
    wanted = set()
    stack = [start_node]

    while stack:
        node = stack.pop()
        if not node or node in wanted:
            continue
        wanted.add(node)
        for child in children_by_node.get(node, set()):
            if child not in wanted:
                stack.append(child)

    return wanted


def _make_placeholders(items):
    return ','.join('?' for _ in items)


def resolve_haplogroup_simple(db_path, user_input):
    # figure out which 2016 haplogroup branch to use
    h = (user_input or '').strip().rstrip('?*~')
    if not h:
        return None, None, False

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    alias_to_node, _, _ = _load_tree_data(cur)

    if '-' not in h:
        direct_node = _resolve_tree_name(h, alias_to_node)
        if direct_node:
            conn.close()
            return direct_node, None, False

    candidate_snp = h.split('-', 1)[1].strip() if '-' in h else h
    cur.execute(
        "SELECT haplogroup FROM snp_reference WHERE UPPER(snp_name) = UPPER(?) LIMIT 1",
        (candidate_snp,),
    )
    row = cur.fetchone()
    if row and row[0]:
        conn.close()
        return row[0], candidate_snp, False

    fallback_used = False
    fallback_node = None
    if '-' in h:
        prefix = h.split('-', 1)[0].strip()
        fallback_node = _resolve_tree_name(prefix, alias_to_node)
        fallback_used = fallback_node is not None
    else:
        fallback_node = _resolve_tree_name(h, alias_to_node)

    conn.close()
    if fallback_node:
        return fallback_node, None, fallback_used
    return None, None, False


def match_by_genotype(db_path, user_haplo, user_genotypes, user_labels=None, limit=1000, min_snps=1):
    logging.basicConfig(filename='match_debug.log', level=logging.DEBUG, format='%(message)s')
    user_labels = user_labels or {}

    # normalize the parsed user SNP names
    normalized_user_genotypes = {}
    normalized_user_labels = {}
    for snp, val in user_genotypes.items():
        base = _strip_plink_suffix(snp)
        normalized_user_genotypes[base.upper()] = val
        normalized_user_labels[base.upper()] = user_labels.get(snp, user_labels.get(base, snp))
    user_genotypes = normalized_user_genotypes
    user_labels = normalized_user_labels

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    alias_to_node, _, children_by_node = _load_tree_data(cur)
    branch_root = _resolve_tree_name(user_haplo, alias_to_node) or user_haplo
    branch_nodes = sorted(_collect_descendants(branch_root, children_by_node))
    if not branch_nodes:
        conn.close()
        return []

    branch_placeholders = _make_placeholders(branch_nodes)

    # load the canonical 2016 SNPs on this branch
    cur.execute(
        f"""
        SELECT DISTINCT snp_name, position
        FROM snp_reference
        WHERE haplogroup IN ({branch_placeholders})
        """,
        branch_nodes,
    )
    relevant_rows = cur.fetchall()
    if not relevant_rows:
        with open('match_debug.log', 'a') as logf:
            logf.write(f"No canonical 2016 SNPs found for haplogroup {branch_root}\n")
        conn.close()
        return []

    relevant_positions = {}
    user_value_by_position = {}
    user_label_by_position = {}
    for row in relevant_rows:
        snp_name = row['snp_name']
        position = row['position']
        if position is None:
            continue
        relevant_positions.setdefault(position, snp_name)

        key = str(snp_name).upper()
        if key in user_genotypes:
            prev = user_value_by_position.get(position)
            cur_val = user_genotypes[key]
            if prev is None or cur_val > prev:
                user_value_by_position[position] = cur_val
                user_label_by_position[position] = user_labels.get(key, key)

    if not user_value_by_position:
        with open('match_debug.log', 'a') as logf:
            logf.write(f"No user SNPs overlap the 2016 branch for {branch_root}\n")
        conn.close()
        return []

    target_positions = sorted(user_value_by_position.keys())
    user_derived_positions = {pos for pos, val in user_value_by_position.items() if val == 2}
    user_derived_total = len(user_derived_positions)
    if user_derived_total == 0:
        with open('match_debug.log', 'a') as logf:
            logf.write(f"No user derived SNPs for haplogroup {branch_root}\n")
        conn.close()
        return []

    # load individuals on the same branch
    cur.execute(
        f"""
        SELECT id, y_haplogroup, y_haplogroup_clean, y_terminal,
               group_id, locality, country, lat, lon,
               date_mean, full_date, full_date_range, n_called
        FROM individuals
        WHERE y_haplogroup_clean IN ({branch_placeholders})
        """,
        branch_nodes,
    )
    ind_meta = {row['id']: dict(row) for row in cur.fetchall()}
    if not ind_meta:
        with open('match_debug.log', 'a') as logf:
            logf.write(f"No individuals found for haplogroup {branch_root}\n")
        conn.close()
        return []

    geno_placeholders = _make_placeholders(target_positions)
    cur.execute(
        f"""
        SELECT g.individual_id, g.position, MAX(g.is_derived) AS is_derived
        FROM genotypes g
        JOIN individuals i ON i.id = g.individual_id
        WHERE i.y_haplogroup_clean IN ({branch_placeholders})
          AND g.position IN ({geno_placeholders})
        GROUP BY g.individual_id, g.position
        """,
        branch_nodes + target_positions,
    )

    ind_genos = {}
    for row in cur.fetchall():
        iid = row['individual_id']
        ind_genos.setdefault(iid, {})[row['position']] = row['is_derived']
    conn.close()

    def shared_mutation_label(pos):
        # show the user's SNP label when we have it
        canonical = str(relevant_positions.get(pos, pos)).upper()
        user_label = user_label_by_position.get(pos)

        if user_label and user_label.startswith('RS'):
            return user_label if user_label == canonical else f"{user_label} ({canonical})"
        if user_label:
            return user_label if user_label == canonical else f"{user_label} ({canonical})"
        return canonical

    results = []
    for iid, genos in ind_genos.items():
        meta = ind_meta.get(iid)
        if not meta:
            continue

        snps_compared = 0
        snps_matched = 0
        shared_mutations = []
        overlap_debug = []

        for pos in target_positions:
            user_val = user_value_by_position[pos]
            sample_val = genos.get(pos)
            if sample_val not in (0, 1):
                continue

            snps_compared += 1
            overlap_debug.append((pos, user_val, sample_val))
            if pos in user_derived_positions and sample_val == 1:
                snps_matched += 1
                shared_mutations.append(shared_mutation_label(pos))

        meta_id = meta.get('id', 'UNKNOWN')
        logging.debug(f"[DEBUG] Candidate {meta_id} SNP overlap:")
        for pos, user_val, sample_val in overlap_debug:
            logging.debug(f"    pos={pos}, user_val={user_val}, sample_val={sample_val}")

        with open('match_debug.log', 'a') as logf:
            logf.write(f"Candidate {meta_id}: snps_compared={snps_compared}, snps_matched={snps_matched}\n")

        if snps_compared < min_snps:
            continue

        match_score = snps_matched / max(user_derived_total, 1)
        derived_agreement = snps_matched / snps_compared if snps_compared > 0 else 0

        results.append({
            **meta,
            'match_score': round(match_score, 4),
            'derived_agreement': round(derived_agreement, 4),
            'snps_compared': snps_compared,
            'snps_matched': snps_matched,
            'user_derived_total': user_derived_total,
            'low_overlap': snps_compared < 5,
            'shared_mutations': shared_mutations,
        })

    results.sort(key=lambda r: (-r['match_score'], -r['snps_compared']))
    with open('match_debug.log', 'a') as logf:
        logf.write(f"Returning {len(results)} results: {[r.get('id', 'UNKNOWN') for r in results]}\n")
    return results[:limit]


def match_with_broadening(db_path, user_haplo, user_genotypes, user_labels=None, limit=1000, min_snps=1):
    # only broaden if the first exact branch gives no matches
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    alias_to_node, parent_by_node, _ = _load_tree_data(cur)
    conn.close()

    haplo = _resolve_tree_name(user_haplo, alias_to_node) or user_haplo
    for attempt in range(20):
        print(f"[DEBUG] Attempt {attempt}: Trying haplogroup '{haplo}'")
        results = match_by_genotype(db_path, haplo, user_genotypes, user_labels=user_labels, limit=limit, min_snps=min_snps)
        print(f"[DEBUG] Results found: {len(results)} for haplogroup '{haplo}'")
        if results:
            return results, haplo, attempt

        parent = parent_by_node.get(haplo)
        if not parent:
            return results, haplo, attempt
        haplo = parent

    return [], haplo, 20
