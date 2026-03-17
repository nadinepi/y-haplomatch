import logging
import os
import re
import sqlite3


SNP24_POSITION_RE = re.compile(r'^snp_24_(\d+)(?:_|$)', re.IGNORECASE)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TREE_2019_FILE = os.path.join(BASE_DIR, '..', 'data', 'AncientYDNA', 'chrY_hGrpTree_isogg2019.txt')

TREE_2019_CACHE = None
TREE_2019_PARENT_CACHE = None


def _split_aliases(text):
    # some tree nodes use slash-separated aliases
    return [part.strip() for part in str(text).split('/') if part.strip() and part.strip() != '#']


def _extract_marker_position(snp_name):
    # pull the Y position from PLINK marker names
    m = SNP24_POSITION_RE.match(str(snp_name))
    if m:
        return int(m.group(1))
    return None


def _normalize_position_key(value):
    # accept int positions, plain digit strings, or PLINK marker names
    if isinstance(value, int):
        return value
    text = str(value).strip()
    if text.isdigit():
        return int(text)
    return _extract_marker_position(text)


def _load_tree_data_2016(cur):
    # load the 2016 tree aliases from the database
    cur.execute("SELECT child, parent FROM haplogroup_tree")
    tree_rows = cur.fetchall()

    alias_to_node = {}
    for child_raw, parent_raw in tree_rows:
        child_aliases = _split_aliases(child_raw)
        if not child_aliases:
            continue

        child = child_aliases[0]
        for alias in child_aliases:
            alias_to_node[alias.upper()] = child

    return alias_to_node


def _load_tree_parents_2016(cur):
    # load canonical parent links for the 2016 tree
    cur.execute("SELECT child, parent FROM haplogroup_tree")
    parent_by_node = {}
    for child_raw, parent_raw in cur.fetchall():
        child_aliases = _split_aliases(child_raw)
        parent_aliases = _split_aliases(parent_raw)
        if not child_aliases:
            continue
        child = child_aliases[0]
        parent = parent_aliases[0] if parent_aliases else None
        parent_by_node[child] = parent
    return parent_by_node


def _load_tree_data_2019():
    # load the scraped 2019 tree from file
    global TREE_2019_CACHE
    if TREE_2019_CACHE is not None:
        return TREE_2019_CACHE

    alias_to_node = {}
    with open(TREE_2019_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            child_aliases = _split_aliases(parts[0])
            if not child_aliases:
                continue
            child = child_aliases[0]
            for alias in child_aliases:
                alias_to_node[alias.upper()] = child

    TREE_2019_CACHE = alias_to_node
    return alias_to_node


def _load_tree_parents_2019():
    # load canonical parent links for the 2019 tree
    global TREE_2019_PARENT_CACHE
    if TREE_2019_PARENT_CACHE is not None:
        return TREE_2019_PARENT_CACHE

    parent_by_node = {}
    with open(TREE_2019_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            child_aliases = _split_aliases(parts[0])
            parent_aliases = _split_aliases(parts[1])
            if not child_aliases:
                continue
            child = child_aliases[0]
            parent = parent_aliases[0] if parent_aliases else None
            parent_by_node[child] = parent

    TREE_2019_PARENT_CACHE = parent_by_node
    return parent_by_node


def _resolve_tree_name(name, alias_to_node):
    # map any tree alias to its canonical node name
    if not name:
        return None
    return alias_to_node.get(name.strip().upper())


def _make_placeholders(items):
    return ','.join('?' for _ in items)


def _major_haplogroup_cap(haplo):
    # stop broadening once we reach the main letter branch like R
    if not haplo:
        return None
    m = re.match(r'[A-Z]', haplo.upper())
    if m:
        return m.group(0)
    return haplo


def _broaden_one_level(haplo, cap):
    # make the haplogroup one step broader by trimming the last character
    if not haplo or haplo == cap:
        return None
    if len(haplo) <= len(cap):
        return None
    return haplo[:-1]


def _find_branch_or_ancestor(haplo, alias_to_node):
    # find the exact branch in a tree, or the nearest broader branch that exists there
    if not haplo:
        return None

    direct = _resolve_tree_name(haplo, alias_to_node)
    if direct:
        return direct

    current = haplo
    cap = _major_haplogroup_cap(haplo)
    while True:
        current = _broaden_one_level(current, cap)
        if not current:
            return None
        direct = _resolve_tree_name(current, alias_to_node)
        if direct:
            return direct


def _tree_depth(node, parent_by_node):
    # count how many parent steps we can walk from this node
    depth = 0
    current = node
    seen = set()
    while current and current not in seen and current != '#':
        seen.add(current)
        current = parent_by_node.get(current)
        depth += 1
    return depth


def _same_major_branch(a, b):
    # keep suggestions inside the same top-level branch like R, I, Q
    if not a or not b:
        return False
    return _major_haplogroup_cap(a) == _major_haplogroup_cap(b)


def _relation_between_branches(base_haplo, other_haplo):
    # explain how the suggested branch relates to the current branch
    if not base_haplo or not other_haplo:
        return 'different'
    if other_haplo == base_haplo:
        return 'same'
    if other_haplo.startswith(base_haplo):
        return 'downstream'
    if base_haplo.startswith(other_haplo):
        return 'broader'
    return 'different'


def _lineage_relation(user_haplo, sample_haplo):
    # classify how the ancient sample sits relative to the user branch
    if not sample_haplo:
        return None
    if sample_haplo == user_haplo:
        return 'exact'
    if user_haplo.startswith(sample_haplo):
        return 'broader'
    if sample_haplo.startswith(user_haplo):
        return 'downstream'
    return None


def _lineage_rank(relation):
    # exact first, then broader, then downstream, then unlabeled position matches
    if relation == 'exact':
        return 0
    if relation == 'broader':
        return 1
    if relation == 'downstream':
        return 2
    return 3


def _resolve_haplogroup_for_version(cur, user_input, tree_version):
    # try to resolve the input inside one tree/SNP system
    h = (user_input or '').strip().rstrip('?*~')
    if not h:
        return None, None, False

    if tree_version == '2019':
        alias_to_node = _load_tree_data_2019()
        cur.execute("SELECT snp_name, haplogroup FROM snp_reference WHERE source = '2019'")
        snp_to_haplo = {name.upper(): haplo for name, haplo in cur.fetchall() if haplo}
    else:
        alias_to_node = _load_tree_data_2016(cur)
        cur.execute("SELECT snp_name, haplogroup FROM snp_reference WHERE source = '2016'")
        snp_to_haplo = {name.upper(): haplo for name, haplo in cur.fetchall() if haplo}

    if '-' not in h:
        direct_node = _resolve_tree_name(h, alias_to_node)
        if direct_node:
            return direct_node, None, False, 3

    candidate_snp = h.split('-', 1)[1].strip() if '-' in h else h
    snp_haplo = snp_to_haplo.get(candidate_snp.upper())
    if snp_haplo:
        return snp_haplo, candidate_snp, False, 2

    fallback_used = False
    fallback_node = None
    if '-' in h:
        prefix = h.split('-', 1)[0].strip()
        fallback_node = _resolve_tree_name(prefix, alias_to_node)
        fallback_used = fallback_node is not None
    else:
        fallback_node = _resolve_tree_name(h, alias_to_node)

    if fallback_node:
        return fallback_node, None, fallback_used, 1
    return None, None, False, 0


def resolve_haplogroup_simple(db_path, user_input):
    # resolve the input against 2019 first, then 2016
    h = (user_input or '').strip().rstrip('?*~')
    if not h:
        return None, None, False, None

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    best = None
    for tree_version in ('2019', '2016'):
        resolved, snp, fallback, score = _resolve_haplogroup_for_version(cur, h, tree_version)
        if not resolved:
            continue

        candidate = (score, 1 if tree_version == '2019' else 0, resolved, snp, fallback, tree_version)
        if best is None or candidate > best:
            best = candidate

    conn.close()
    if best is None:
        return None, None, False, None
    _, _, resolved, snp, fallback, tree_version = best
    return resolved, snp, fallback, tree_version


def _load_relevant_rows_2016(cur, branch_root):
    # load 2016 SNP rows for one branch and its descendants
    cur.execute(
        """
        SELECT DISTINCT snp_name, position, ref_allele, alt_allele
        FROM snp_reference
        WHERE source = '2016'
          AND (haplogroup = ? OR haplogroup LIKE ?)
        """,
        (branch_root, f"{branch_root}%"),
    )
    return [dict(row) for row in cur.fetchall()]


def _load_relevant_rows_2019(cur, branch_root):
    # load 2019 SNP rows for one branch and its descendants
    cur.execute(
        """
        SELECT DISTINCT snp_name, position, ref_allele, alt_allele
        FROM snp_reference
        WHERE source = '2019'
          AND (haplogroup = ? OR haplogroup LIKE ?)
        """,
        (branch_root, f"{branch_root}%"),
    )
    return [dict(row) for row in cur.fetchall()]


def _load_relevant_rows_union(cur, branch_root, tree_version):
    # use the selected tree first, then add the other tree if the same branch exists there too
    rows = []
    seen = set()

    alias_2016 = _load_tree_data_2016(cur)
    alias_2019 = _load_tree_data_2019()

    root_2016 = _find_branch_or_ancestor(branch_root, alias_2016)
    root_2019 = _find_branch_or_ancestor(branch_root, alias_2019)

    source_order = ['2019', '2016'] if tree_version == '2019' else ['2016', '2019']
    for source in source_order:
        if source == '2016' and root_2016:
            source_rows = _load_relevant_rows_2016(cur, root_2016)
        elif source == '2019' and root_2019:
            source_rows = _load_relevant_rows_2019(cur, root_2019)
        else:
            continue

        for row in source_rows:
            key = (
                row['snp_name'],
                row['position'],
                row.get('ref_allele'),
                row.get('alt_allele'),
            )
            if key in seen:
                continue
            seen.add(key)
            rows.append(row)

    return rows


def suggest_haplogroup_from_user_data(db_path, resolved_haplo, user_alleles, tree_version):
    # find a better-supported branch in the selected tree using the user's SNP overlap
    if not resolved_haplo:
        return None

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    if tree_version == '2019':
        alias_to_node = _load_tree_data_2019()
        parent_by_node = _load_tree_parents_2019()
    else:
        alias_to_node = _load_tree_data_2016(cur)
        parent_by_node = _load_tree_parents_2016(cur)

    branch_root = _find_branch_or_ancestor(resolved_haplo, alias_to_node) or resolved_haplo

    normalized_alleles = {}
    for key, allele in user_alleles.items():
        pos = _normalize_position_key(key)
        if pos is None:
            continue
        normalized_alleles[pos] = str(allele).upper()
    if not normalized_alleles:
        conn.close()
        return None

    pos_placeholders = _make_placeholders(normalized_alleles)
    cur.execute(
        f"""
        SELECT snp_name, position, ref_allele, haplogroup, alt_allele
        FROM snp_reference
        WHERE source = ?
          AND position IN ({pos_placeholders})
        """,
        [tree_version, *sorted(normalized_alleles.keys())],
    )

    derived_support = {}
    ancestral_support = {}
    for row in cur.fetchall():
        haplo = row['haplogroup']
        pos = row['position']
        allele = normalized_alleles.get(pos)
        ref = (row['ref_allele'] or '').upper()
        alt = (row['alt_allele'] or '').upper()
        if not haplo or not allele or not ref or not alt:
            continue
        if not _same_major_branch(branch_root, haplo):
            continue

        if allele == alt:
            derived_support[haplo] = derived_support.get(haplo, 0) + 1
        elif allele == ref:
            ancestral_support[haplo] = ancestral_support.get(haplo, 0) + 1

    conn.close()

    current_derived = derived_support.get(branch_root, 0)
    current_ancestral = ancestral_support.get(branch_root, 0)
    current_score = current_derived - current_ancestral

    best = None
    candidate_haplos = set(derived_support) | set(ancestral_support)
    for haplo in candidate_haplos:
        derived = derived_support.get(haplo, 0)
        ancestral = ancestral_support.get(haplo, 0)
        score = derived - ancestral
        depth = _tree_depth(haplo, parent_by_node)

        if haplo == branch_root:
            continue
        if derived < 2:
            continue
        if score <= 0:
            continue

        candidate = (score, derived, depth, haplo, ancestral)
        if best is None or candidate > best:
            best = candidate

    if best is None:
        return None

    score, derived, depth, haplo, ancestral = best
    if score < current_score:
        return None
    if score == current_score and depth <= _tree_depth(branch_root, parent_by_node):
        return None

    relation = _relation_between_branches(branch_root, haplo)
    # Only suggest if the relation is 'exact' (i.e., exact branch match)
    if relation != 'same':
        return None

    return {
        'haplogroup': haplo,
        'tree_version': tree_version,
        'derived_support': derived,
        'ancestral_support': ancestral,
        'support_score': score,
        'relation_to_resolved': relation,
    }


def match_by_genotype(db_path, user_haplo, user_alleles, user_labels=None, limit=1000, min_snps=1, tree_version='2016'):
    logging.basicConfig(filename='match_debug.log', level=logging.DEBUG, format='%(message)s')
    user_labels = user_labels or {}

    normalized_alleles = {}
    normalized_labels = {}
    for key, allele in user_alleles.items():
        pos = _normalize_position_key(key)
        if pos is None:
            continue
        normalized_alleles[pos] = str(allele).upper()
        normalized_labels[pos] = user_labels.get(key, user_labels.get(str(key), str(key)))

    user_alleles = normalized_alleles
    user_labels = normalized_labels

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    if tree_version == '2019':
        alias_to_node = _load_tree_data_2019()
    else:
        alias_to_node = _load_tree_data_2016(cur)
    branch_root = _resolve_tree_name(user_haplo, alias_to_node) or user_haplo

    # load relevant branch SNPs from both 2016 and 2019 when the branch exists in both
    relevant_rows = _load_relevant_rows_union(cur, branch_root, tree_version)

    if not relevant_rows:
        with open('match_debug.log', 'a') as logf:
            logf.write(f"No SNPs found for haplogroup {branch_root} in {tree_version}\n")
        conn.close()
        return []

    relevant_positions = {}
    user_value_by_position = {}
    user_label_by_position = {}
    for row in relevant_rows:
        pos = row['position']
        if pos is None or pos not in user_alleles:
            continue

        allele = user_alleles[pos]
        ref = (row.get('ref_allele') or '').upper()
        alt = (row.get('alt_allele') or '').upper()
        if not ref or not alt:
            continue

        relevant_positions.setdefault(pos, row['snp_name'])

        if allele == alt:
            user_value_by_position[pos] = 2
            user_label_by_position[pos] = user_labels.get(pos, row['snp_name'])
        elif allele == ref and pos not in user_value_by_position:
            user_value_by_position[pos] = 0
            user_label_by_position[pos] = user_labels.get(pos, row['snp_name'])

    if not user_value_by_position:
        with open('match_debug.log', 'a') as logf:
            logf.write(f"No user SNPs overlap the {tree_version} branch for {branch_root}\n")
        conn.close()
        return []

    target_positions = sorted(user_value_by_position.keys())
    user_derived_positions = {pos for pos, val in user_value_by_position.items() if val == 2}
    user_derived_total = len(user_derived_positions)
    if user_derived_total == 0:
        with open('match_debug.log', 'a') as logf:
            logf.write(f"No user derived SNPs for haplogroup {branch_root} in {tree_version}\n")
        conn.close()
        return []

    # only compare against ancient individuals on the same lineage path
    cur.execute(
        """
        SELECT id, y_haplogroup, y_haplogroup_clean, y_terminal,
               group_id, locality, country, lat, lon,
               date_mean, full_date, full_date_range, n_called
        FROM individuals
        WHERE y_haplogroup_clean IS NOT NULL
          AND (
            y_haplogroup_clean = ?
            OR y_haplogroup_clean LIKE ?
            OR ? LIKE y_haplogroup_clean || '%'
          )
        """,
        (branch_root, f"{branch_root}%", branch_root),
    )
    ind_meta = {}
    for row in cur.fetchall():
        meta = dict(row)
        relation = _lineage_relation(branch_root, meta.get('y_haplogroup_clean'))
        if not relation:
            continue
        meta['lineage_relation'] = relation
        ind_meta[meta['id']] = meta

    if not ind_meta:
        with open('match_debug.log', 'a') as logf:
            logf.write(f"No lineage-relevant individuals found for haplogroup {branch_root}\n")
        conn.close()
        return []

    geno_placeholders = _make_placeholders(target_positions)
    cur.execute(
        f"""
        SELECT g.individual_id, g.position, MAX(g.is_derived) AS is_derived
        FROM genotypes g
        WHERE g.position IN ({geno_placeholders})
        GROUP BY g.individual_id, g.position
        """,
        target_positions,
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
        if user_label:
            return user_label if str(user_label).upper() == canonical else f"{user_label} ({canonical})"
        return canonical

    results = []
    for iid, genos in ind_genos.items():
        meta = ind_meta.get(iid)
        if not meta:
            continue

        snps_compared = 0
        snps_matched = 0
        shared_mutations = []

        for pos in target_positions:
            sample_val = genos.get(pos)
            if sample_val not in (0, 1):
                continue

            snps_compared += 1
            if pos in user_derived_positions and sample_val == 1:
                snps_matched += 1
                shared_mutations.append(shared_mutation_label(pos))

        if snps_compared < min_snps or snps_matched == 0:
            continue

        match_score = snps_matched / max(user_derived_total, 1)
        derived_agreement = snps_matched / snps_compared if snps_compared > 0 else 0

        results.append({
            **meta,
            'lineage_rank': _lineage_rank(meta.get('lineage_relation')),
            'match_score': round(match_score, 4),
            'derived_agreement': round(derived_agreement, 4),
            'snps_compared': snps_compared,
            'snps_matched': snps_matched,
            'user_derived_total': user_derived_total,
            'shared_mutations': shared_mutations,
        })

    results.sort(key=lambda r: (r['lineage_rank'], -r['match_score'], -r['snps_compared']))
    with open('match_debug.log', 'a') as logf:
        logf.write(f"Returning {len(results)} results for {branch_root} in {tree_version}\n")
    return results[:limit]


def match_with_broadening(db_path, user_haplo, user_alleles, user_labels=None, limit=1000, min_snps=1, tree_version='2016'):
    # only broaden if the first exact branch gives no matches
    if tree_version == '2019':
        alias_to_node = _load_tree_data_2019()
    else:
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        alias_to_node = _load_tree_data_2016(cur)
        conn.close()

    haplo = _resolve_tree_name(user_haplo, alias_to_node) or user_haplo
    cap = _major_haplogroup_cap(haplo)
    for attempt in range(20):
        print(f"[DEBUG] Attempt {attempt}: Trying haplogroup '{haplo}' in {tree_version}")
        results = match_by_genotype(
            db_path,
            haplo,
            user_alleles,
            user_labels=user_labels,
            limit=limit,
            min_snps=min_snps,
            tree_version=tree_version,
        )
        print(f"[DEBUG] Results found: {len(results)} for haplogroup '{haplo}'")
        if results:
            return results, haplo, attempt

        next_haplo = _broaden_one_level(haplo, cap)
        if not next_haplo:
            return results, haplo, attempt
        haplo = next_haplo

    return [], haplo, 20
