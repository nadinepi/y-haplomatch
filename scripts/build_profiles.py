import csv
import os
import re
import sqlite3
import openpyxl


def is_invalid_haplo(h):
    # check if haplogroup value is missing or uninformative
    if h is None:
        return True
    h = re.sub(r'\s+', ' ', str(h).strip().lower())
    if h in ('', 'na', 'n/a', 'nan', '?', '*', '~', 'none', '..', '.', '0'):
        return True
    if h.startswith('n/a'):
        return True
    if 'not published' in h:
        return True
    if h in ('sex unknown', 'unknown'):
        return True
    return False


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'plink_out', 'aadr_chrY.raw')
ANNO_CSV_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'AADR_54.1', 'AADR Annotations 2025.csv')
ANNO_XLSX_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'AADR_54.1', 'AADR Annotations 2025.xlsx')
DB_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'ydna.db')
TREE_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'AncientYDNA', 'chrY_hGrpTree_isogg2016.txt')
LOCUS_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'AncientYDNA', 'chrY_locusFile_b37_isogg2016.txt')
TREE_2019_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'AncientYDNA', 'chrY_hGrpTree_isogg2019.txt')
SNP_2019_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'AncientYDNA', 'snpFile_b37_isogg2019.txt')
ANCIENT_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'AADR_54.1', 'Ancient_samples.txt')

SNP24_POSITION_RE = re.compile(r'^snp_24_(\d+)(?:_|$)', re.IGNORECASE)


def require_file(path):
    # stop early if an important input file is missing
    if not os.path.exists(path):
        raise FileNotFoundError(f'Missing required file: {path}')


def find_annotation_file():
    # accept either the csv or the original xlsx
    if os.path.exists(ANNO_CSV_FILE):
        return ANNO_CSV_FILE
    if os.path.exists(ANNO_XLSX_FILE):
        return ANNO_XLSX_FILE
    raise FileNotFoundError(f'Missing required file: {ANNO_CSV_FILE}')


def normalize_plink_marker_name(snp_name):
    # strip the extra PLINK assay suffix
    if snp_name.lower().startswith('snp_24_'):
        pos = extract_position_from_marker(snp_name)
        if pos is not None:
            return f'snp_24_{pos}'

    if '_' in snp_name:
        head, tail = snp_name.rsplit('_', 1)
        if tail.isdigit():
            return head

    return snp_name


def extract_position_from_marker(marker_name):
    # pull the Y position out of a marker name
    m = SNP24_POSITION_RE.match(marker_name)
    if m:
        return int(m.group(1))
    return None


def split_aliases(text):
    # some ISOGG nodes are written like LT/K1
    return [part.strip() for part in str(text).split('/') if part.strip() and part.strip() != '#']


def load_haplogroup_tree_from_file(tree_path):
    # load the tree file exactly as it is written
    tree = []
    with open(tree_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            tree.append((parts[0].strip(), parts[1].strip()))
    return tree


def build_tree_maps(tree_rows):
    # keep one canonical name per node, but also remember aliases
    alias_to_canonical = {}
    canonical_tree_rows = []
    haplo_alias_rows = []
    seen_tree = set()
    seen_alias = set()

    for child_raw, parent_raw in tree_rows:
        child_aliases = split_aliases(child_raw)
        parent_aliases = split_aliases(parent_raw)
        if not child_aliases:
            continue

        child_name = child_aliases[0]
        parent_name = parent_aliases[0] if parent_aliases else None

        if (child_name, parent_name) not in seen_tree:
            canonical_tree_rows.append((child_name, parent_name))
            seen_tree.add((child_name, parent_name))

        for alias in child_aliases:
            alias_to_canonical[alias.upper()] = child_name
            if alias.upper() not in seen_alias:
                haplo_alias_rows.append((alias, child_name))
                seen_alias.add(alias.upper())

        for alias in parent_aliases:
            alias_to_canonical[alias.upper()] = parent_name
            if alias.upper() not in seen_alias:
                haplo_alias_rows.append((alias, parent_name))
                seen_alias.add(alias.upper())

    return alias_to_canonical, canonical_tree_rows, haplo_alias_rows


def normalize_haplogroup(raw, terminal_raw, alias_to_canonical):
    # try to turn noisy text into one 2016 tree node
    for source in (raw, terminal_raw):
        if is_invalid_haplo(source):
            continue

        text = str(source).strip()
        text = text.replace('–', '-').replace('—', '-').replace('−', '-')
        text = re.sub(r'\([^)]*\)', ' ', text)
        text = text.replace(';', ' ').replace(',', ' ').replace('/', ' ')

        for token in text.split():
            token = token.strip().strip("()[]{}")
            token = token.rstrip('?*~')
            if not token:
                continue

            if token.lower().startswith('x') and len(token) > 1:
                token = token[1:]

            canonical = alias_to_canonical.get(token.upper())
            if canonical:
                return canonical

    return None


def load_canonical_snp_reference(alias_to_canonical):
    # use the 2016 locus file as the SNP list
    snp_to_haplo = {}
    snp_rows = []
    rows_by_pos = {}

    with open(LOCUS_FILE) as f:
        for line in f:
            if not line.strip() or line.startswith('#') or line.startswith('locName'):
                continue
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            snp_name = row[0].strip()
            position = row[1].strip()
            ref_allele = row[2].strip().upper()
            raw_haplo = row[3].strip()
            alt_allele = row[4].strip().upper()

            if not snp_name or not position.isdigit():
                continue

            cleaned_haplo = normalize_haplogroup(raw_haplo, None, alias_to_canonical)
            if not cleaned_haplo:
                continue

            snp_row = (snp_name, int(position), ref_allele, cleaned_haplo, alt_allele, '2016')
            snp_rows.append(snp_row)
            snp_to_haplo[snp_name.upper()] = cleaned_haplo
            rows_by_pos.setdefault(int(position), []).append(snp_row)

    return snp_to_haplo, snp_rows, rows_by_pos


def load_2019_snp_reference(alias_to_canonical):
    # use the 2019 SNP file as an extra reference set
    snp_to_haplo = {}
    snp_rows = []
    rows_by_pos = {}

    with open(SNP_2019_FILE) as f:
        for line in f:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            snp_name = row[0].strip()
            position = row[1].strip()
            ref_allele = row[2].strip().upper()
            raw_haplo = row[3].strip()
            alt_allele = row[4].strip().upper()

            if not snp_name or not position.isdigit():
                continue

            cleaned_haplo = normalize_haplogroup(raw_haplo, None, alias_to_canonical)
            if not cleaned_haplo:
                continue

            snp_row = (snp_name, int(position), ref_allele, cleaned_haplo, alt_allele, '2019')
            snp_rows.append(snp_row)
            snp_to_haplo[snp_name.upper()] = cleaned_haplo
            rows_by_pos.setdefault(int(position), []).append(snp_row)

    return snp_to_haplo, snp_rows, rows_by_pos


def load_ancient_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                ids.add(parts[1])
            else:
                ids.add(parts[0])
    return ids


def load_annotations(annotation_file):
    if annotation_file.endswith('.csv'):
        with open(annotation_file, newline='', encoding='utf-8') as f:
            rows = list(csv.reader(f))
    else:
        wb = openpyxl.load_workbook(annotation_file, read_only=True, data_only=True)
        ws = wb.active
        rows = [[cell for cell in row] for row in ws.iter_rows(values_only=True)]

    headers = [str(h).lower().strip() if h else '' for h in rows[0]]

    # AADR column names are long so match on keywords
    col_gid = col_group = col_locality = col_country = None
    col_lat = col_lon = col_date = col_fulldate = None
    col_yhaplo = col_yterminal = None

    for i, h in enumerate(headers):
        if 'genetic id' in h:
            col_gid = i
        elif 'group id' in h:
            col_group = i
        elif 'locality' in h:
            col_locality = i
        elif 'political entity' in h:
            col_country = i
        elif 'lat' in h and 'long' not in h:
            col_lat = i
        elif 'long' in h:
            col_lon = i
        elif 'date mean' in h:
            col_date = i
        elif 'full date' in h:
            col_fulldate = i
        elif 'y haplogroup' in h and 'isogg' in h:
            col_yhaplo = i
        elif 'y haplogroup' in h and 'terminal' in h:
            col_yterminal = i

    annotations = {}
    for row in rows[1:]:
        if col_gid is None or len(row) <= col_gid:
            continue
        gid = str(row[col_gid]).strip() if row[col_gid] else ''
        if not gid:
            continue
        annotations[gid] = {
            'group': str(row[col_group]).strip() if col_group is not None and len(row) > col_group and row[col_group] else None,
            'locality': str(row[col_locality]).strip() if col_locality is not None and len(row) > col_locality and row[col_locality] else None,
            'country': str(row[col_country]).strip() if col_country is not None and len(row) > col_country and row[col_country] else None,
            'lat': str(row[col_lat]).strip() if col_lat is not None and len(row) > col_lat and row[col_lat] else None,
            'lon': str(row[col_lon]).strip() if col_lon is not None and len(row) > col_lon and row[col_lon] else None,
            'date_mean': str(row[col_date]).strip() if col_date is not None and len(row) > col_date and row[col_date] else None,
            'full_date': str(row[col_fulldate]).strip() if col_fulldate is not None and len(row) > col_fulldate and row[col_fulldate] else None,
            'y_haplo': str(row[col_yhaplo]).strip() if col_yhaplo is not None and len(row) > col_yhaplo and row[col_yhaplo] else None,
            'y_terminal': str(row[col_yterminal]).strip() if col_yterminal is not None and len(row) > col_yterminal and row[col_yterminal] else None,
        }
    return annotations


def extract_date_range(full_date):
    if not full_date:
        return None
    if 'present' in full_date.lower():
        return '1950 CE'

    m = re.search(r'(\d{3,5}-\d{3,5} (?:cal)?(?:BCE|CE))', full_date)
    if m:
        return m.group(1).replace('calBCE', 'BCE').replace('calCE', 'CE').replace('cal BCE', 'BCE').replace('cal CE', 'CE')

    m = re.search(r'(\d{3,5} (?:cal)?(?:BCE|CE))\s*-\s*(\d{3,5} (?:cal)?(?:BCE|CE))', full_date)
    if m:
        left = m.group(1).replace('calBCE', 'BCE').replace('calCE', 'CE').replace('cal BCE', 'BCE').replace('cal CE', 'CE')
        right = m.group(2).replace('calBCE', 'BCE').replace('calCE', 'CE').replace('cal BCE', 'BCE').replace('cal CE', 'CE')
        return f"{left} - {right}"

    return full_date


def create_db(db_path):
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute('''CREATE TABLE individuals (
        id TEXT PRIMARY KEY,
        y_haplogroup TEXT,
        y_haplogroup_clean TEXT,
        y_terminal TEXT,
        group_id TEXT,
        locality TEXT,
        country TEXT,
        lat REAL,
        lon REAL,
        date_mean REAL,
        full_date TEXT,
        full_date_range TEXT,
        n_called INTEGER
    )''')

    cur.execute('''CREATE TABLE genotypes (
        individual_id TEXT,
        snp_name TEXT,
        position INTEGER,
        value INTEGER,
        is_derived INTEGER,
        PRIMARY KEY (individual_id, snp_name)
    )''')

    cur.execute('''CREATE TABLE haplogroup_tree (
        child TEXT PRIMARY KEY,
        parent TEXT
    )''')

    cur.execute('''CREATE TABLE haplogroup_aliases (
        alias_name TEXT PRIMARY KEY,
        canonical_name TEXT
    )''')

    cur.execute('''CREATE TABLE snp_reference (
        snp_name TEXT,
        position INTEGER,
        ref_allele TEXT,
        haplogroup TEXT,
        alt_allele TEXT,
        source TEXT
    )''')

    # add a few indexes for lookups
    cur.execute('CREATE INDEX idx_geno_ind ON genotypes(individual_id)')
    cur.execute('CREATE INDEX idx_geno_pos ON genotypes(position)')
    cur.execute('CREATE INDEX idx_ind_haplo ON individuals(y_haplogroup_clean)')
    cur.execute('CREATE INDEX idx_tree_parent ON haplogroup_tree(parent)')
    cur.execute('CREATE INDEX idx_haplo_alias_canon ON haplogroup_aliases(canonical_name)')
    cur.execute('CREATE INDEX idx_snpref_haplo ON snp_reference(haplogroup)')
    cur.execute('CREATE INDEX idx_snpref_pos ON snp_reference(position)')
    cur.execute('CREATE INDEX idx_snpref_source ON snp_reference(source)')
    cur.execute('CREATE INDEX idx_snpref_source_haplo ON snp_reference(source, haplogroup)')

    conn.commit()
    return conn


def main():
    # stop early if the required files are not in place
    annotation_file = find_annotation_file()
    for path in (RAW_FILE, annotation_file, ANCIENT_FILE, TREE_FILE, LOCUS_FILE, TREE_2019_FILE, SNP_2019_FILE):
        require_file(path)

    os.makedirs(os.path.dirname(DB_FILE), exist_ok=True)

    print('Loading 2016 haplogroup tree...')
    raw_tree_rows = load_haplogroup_tree_from_file(TREE_FILE)
    alias_to_canonical, canonical_tree_rows, haplo_alias_rows = build_tree_maps(raw_tree_rows)
    print(f'  {len(canonical_tree_rows)} canonical tree edges')
    print(f'  {len(haplo_alias_rows)} haplogroup aliases')

    print('Loading canonical SNP reference from the 2016 locus file...')
    snp_to_haplo_2016, snp_rows_2016, rows_by_pos_2016 = load_canonical_snp_reference(alias_to_canonical)
    print(f'  {len(snp_rows_2016)} canonical SNP rows')

    print('Loading 2019 haplogroup tree...')
    raw_tree_rows_2019 = load_haplogroup_tree_from_file(TREE_2019_FILE)
    alias_to_canonical_2019, _, _ = build_tree_maps(raw_tree_rows_2019)

    print('Loading extra SNP reference from the 2019 SNP file...')
    snp_to_haplo_2019, snp_rows_2019, rows_by_pos_2019 = load_2019_snp_reference(alias_to_canonical_2019)
    print(f'  {len(snp_rows_2019)} extra 2019 SNP rows')

    snp_rows = snp_rows_2016 + snp_rows_2019
    ref_pos_by_name = {}
    for row in snp_rows:
        ref_pos_by_name.setdefault(row[0].upper(), row[1])
    ref_positions = {row[1] for row in snp_rows}
    ref_names = set(snp_to_haplo_2016.keys()) | set(snp_to_haplo_2019.keys())

    print('Loading annotations...')
    annotations = load_annotations(annotation_file)
    print(f'  {len(annotations)} individuals annotated')

    print('Loading ancient sample IDs...')
    ancient_ids = load_ancient_ids(ANCIENT_FILE)
    print(f'  {len(ancient_ids)} ancient IDs')

    print('Creating database...')
    conn = create_db(DB_FILE)
    cur = conn.cursor()

    print('Inserting tree and SNP tables...')
    cur.executemany('INSERT INTO haplogroup_tree VALUES (?,?)', canonical_tree_rows)
    cur.executemany('INSERT INTO haplogroup_aliases VALUES (?,?)', haplo_alias_rows)
    cur.executemany('INSERT INTO snp_reference VALUES (?,?,?,?,?,?)', snp_rows)
    conn.commit()

    print('Loading plink .raw file...')
    with open(RAW_FILE) as f:
        header = f.readline().strip().split()
        snp_names = header[6:]

        # only keep Y SNPs that connect to the 2016 or 2019 SNP list
        snp_is_relevant = []
        for snp_name in snp_names:
            base = normalize_plink_marker_name(snp_name)
            pos = extract_position_from_marker(base)
            snp_is_relevant.append(base.upper() in ref_names or pos in ref_positions)

        n_relevant = sum(snp_is_relevant)
        print(f'  {len(snp_names)} Y-SNPs total, {n_relevant} relevant to the 2016/2019 SNP list')

        inserted = 0
        skipped = 0
        skipped_non_ancient = 0
        norm_matched = 0
        malformed_rows = 0

        for line in f:
            if not line.strip():
                continue

            row = line.split()
            if len(row) != len(header):
                malformed_rows += 1
                continue

            iid = row[1]
            if iid not in ancient_ids:
                skipped += 1
                skipped_non_ancient += 1
                continue

            info = annotations.get(iid, {})
            y_haplo = info.get('y_haplo')
            if is_invalid_haplo(y_haplo):
                skipped += 1
                continue

            geno_rows = []
            for snp_name, val, relevant in zip(snp_names, row[6:], snp_is_relevant):
                if not relevant or val == 'NA':
                    continue

                geno_val = int(val)
                norm_snp_name = normalize_plink_marker_name(snp_name)
                position = extract_position_from_marker(norm_snp_name)
                if position is None:
                    position = ref_pos_by_name.get(norm_snp_name.upper())

                is_derived = 1 if geno_val == 2 else 0
                geno_rows.append((iid, norm_snp_name, position, geno_val, is_derived))

            if not geno_rows:
                skipped += 1
                continue

            lat = None
            lon = None
            try:
                lat = float(info['lat']) if info.get('lat') else None
                lon = float(info['lon']) if info.get('lon') else None
            except (ValueError, TypeError):
                pass

            date_mean = None
            try:
                date_mean = float(info['date_mean']) if info.get('date_mean') else None
            except (ValueError, TypeError):
                pass

            y_clean = normalize_haplogroup(y_haplo, info.get('y_terminal'), alias_to_canonical)
            if y_clean:
                norm_matched += 1

            full_date = info.get('full_date')
            full_date_range = extract_date_range(full_date)

            cur.execute(
                'INSERT INTO individuals VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',
                (
                    iid, y_haplo, y_clean, info.get('y_terminal'), info.get('group'),
                    info.get('locality'), info.get('country'),
                    lat, lon, date_mean, full_date, full_date_range, len(geno_rows)
                )
            )
            cur.executemany('INSERT INTO genotypes VALUES (?,?,?,?,?)', geno_rows)

            inserted += 1
            if inserted % 200 == 0:
                print(f'  ... {inserted} individuals inserted')
                conn.commit()

    conn.commit()
    conn.close()

    print(f'\nDone! {inserted} individuals inserted, {skipped} skipped')
    print(f'  {skipped_non_ancient} skipped because they are not in Ancient_samples.txt')
    print(f'  {norm_matched}/{inserted} haplogroups resolved to 2016 tree nodes')
    if malformed_rows:
        print(f'Ignored {malformed_rows} malformed .raw rows')
    print(f'Database saved to {DB_FILE}')


if __name__ == '__main__':
    main()
