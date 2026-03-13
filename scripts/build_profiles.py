#!/usr/bin/env python3
#
# Merge plink-extracted Y chromosome genotypes with AADR annotations
# and store everything in a SQLite database.
# Run extract_y.sh first to get the .raw file.
#
import sqlite3
import os
import re
import csv
import openpyxl
def is_invalid_haplo(h):
    """Return True if haplogroup value is missing or uninformative (AADR style)."""
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
ANNO_FILE = os.path.join(SCRIPT_DIR, '..', 'ProjetcFiles', 'AADR_54.1', 'AADR Annotations 2025.csv')
DB_FILE = os.path.join(SCRIPT_DIR, '..', 'data', 'ydna.db')
LOCUS_FILE = os.path.join(SCRIPT_DIR, '..', 'ProjetcFiles', 'AncientYDNA', 'chrY_locusFile_b37_isogg2016.txt')
ANCIENT_FILE = os.path.join(SCRIPT_DIR, '..', 'ProjetcFiles', 'AADR_54.1', 'Ancient_samples.txt')



def load_snp_reference():
    """Load SNP reference from snpFile_b37_isogg2019.txt.
    Returns:
        snp_to_haplo: dict  snp_name -> ISOGG haplogroup
        snp_rows: list of (snp_name, position, ref_allele, haplogroup, alt_allele)
    """
    snp_to_haplo = {}
    snp_rows = []
    with open(LOCUS_FILE) as f:
        for line in f:
            if not line.strip() or line.startswith('#') or line.startswith('locName'):
                continue
            row = line.strip().split('\t')
            if len(row) < 5:
                continue
            snp_name = row[0].strip()
            position = row[1].strip()
            ref_allele = row[2].strip()
            isogg_haplo = row[3].strip()
            alt_allele = row[4].strip()

            if not position.isdigit():
                continue
            if not isogg_haplo:
                continue
            if snp_name not in snp_to_haplo:
                snp_to_haplo[snp_name] = isogg_haplo
                snp_rows.append((snp_name, int(position), ref_allele, isogg_haplo, alt_allele))

    return snp_to_haplo, snp_rows



def load_haplogroup_tree_from_file(tree_path):
    """Load haplogroup tree from a tab-separated file (child<TAB>parent)."""
    tree = []
    with open(tree_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            child, parent = parts
            tree.append((child.strip(), parent.strip()))
    return tree


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


def normalize_haplogroup(raw, terminal_raw, snp_to_haplo, tree_by_name):
    """Try to convert noisy AADR Y labels to one tree node."""
    snp_lookup = {k.upper(): v for k, v in snp_to_haplo.items()}

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

            candidates = [token]

            if '-' in token:
                candidates.extend([p for p in token.split('-') if p])

            if "'" in token:
                candidates.append(token.split("'", 1)[0])

            for cand in candidates:
                cand_up = cand.strip().upper()
                if not cand_up:
                    continue

                if cand_up in tree_by_name:
                    return tree_by_name[cand_up]

                if cand_up in snp_lookup:
                    mapped_haplo = snp_lookup[cand_up].upper()
                    if mapped_haplo in tree_by_name:
                        return tree_by_name[mapped_haplo]

    return None


import csv
def load_annotations():
    with open(ANNO_FILE, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        headers = [str(h).lower().strip() if h else '' for h in next(reader)]

        # AADR column names are super long so match on keywords
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
        for row in reader:
            vals = row
            if len(vals) <= col_gid:
                continue
            gid = vals[col_gid]
            if not gid:
                continue
            gid = str(gid).strip()
            annotations[gid] = {
                'group':     str(vals[col_group]).strip() if col_group is not None and len(vals) > col_group and vals[col_group] else None,
                'locality':  str(vals[col_locality]).strip() if col_locality is not None and len(vals) > col_locality and vals[col_locality] else None,
                'country':   str(vals[col_country]).strip() if col_country is not None and len(vals) > col_country and vals[col_country] else None,
                'lat':       str(vals[col_lat]).strip() if col_lat is not None and len(vals) > col_lat and vals[col_lat] else None,
                'lon':       str(vals[col_lon]).strip() if col_lon is not None and len(vals) > col_lon and vals[col_lon] else None,
                'date_mean': str(vals[col_date]).strip() if col_date is not None and len(vals) > col_date and vals[col_date] else None,
                'full_date': str(vals[col_fulldate]).strip() if col_fulldate is not None and len(vals) > col_fulldate and vals[col_fulldate] else None,
                'y_haplo':   str(vals[col_yhaplo]).strip() if col_yhaplo is not None and len(vals) > col_yhaplo and vals[col_yhaplo] else None,
                'y_terminal': str(vals[col_yterminal]).strip() if col_yterminal is not None and len(vals) > col_yterminal and vals[col_yterminal] else None,
            }
    return annotations


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
        n_called INTEGER
    )''')

    cur.execute('''CREATE TABLE genotypes (
        individual_id TEXT,
        snp_name TEXT,
        value INTEGER,
        PRIMARY KEY (individual_id, snp_name)
    )''')

    cur.execute('''CREATE TABLE haplogroup_tree (
        child TEXT PRIMARY KEY,
        parent TEXT
    )''')

    cur.execute('''CREATE TABLE snp_reference (
        snp_name TEXT PRIMARY KEY,
        position INTEGER,
        ref_allele TEXT,
        haplogroup TEXT,
        alt_allele TEXT
    )''')

    # for fast lookups
    cur.execute('CREATE INDEX idx_geno_ind ON genotypes(individual_id)')
    cur.execute('CREATE INDEX idx_geno_snp ON genotypes(snp_name)')
    cur.execute('CREATE INDEX idx_ind_haplo ON individuals(y_haplogroup_clean)')
    cur.execute('CREATE INDEX idx_snpref_haplo ON snp_reference(haplogroup)')

    conn.commit()
    return conn

def main():
    os.makedirs(os.path.dirname(DB_FILE), exist_ok=True)

    print("Loading SNP reference from locusFile...")
    snp_to_haplo, snp_rows = load_snp_reference()
    print(f"  {len(snp_to_haplo)} SNP->haplogroup mappings")

    print("Loading haplogroup tree from file...")
    TREE_FILE = os.path.join(SCRIPT_DIR, '..', 'ProjetcFiles', 'AncientYDNA', 'chrY_hGrpTree_isogg2016.txt')
    tree = load_haplogroup_tree_from_file(TREE_FILE)
    print(f"  {len(tree)} tree edges loaded")
    # Build tree_nodes set from file-based tree, including all alternative names
    tree_nodes = set()
    for child, parent in tree:
        for name in child.split('/'):
            tree_nodes.add(name)
        for name in parent.split('/'):
            tree_nodes.add(name)
    tree_by_name = {node.upper(): node for node in tree_nodes}

    print("Loading annotations...")
    annotations = load_annotations()
    print(f"  {len(annotations)} individuals annotated")

    print("Loading ancient sample IDs...")
    ancient_ids = load_ancient_ids(ANCIENT_FILE)
    print(f"  {len(ancient_ids)} ancient IDs")

    print("Loading plink .raw file...")
    with open(RAW_FILE) as f:
        header = f.readline().strip().split()
        snp_names = header[6:]

        # Build set of ISOGG-relevant SNP base names for filtering
        isogg_snp_names = set(snp_to_haplo.keys())
        # Map column index -> True if this SNP is ISOGG-relevant
        snp_is_relevant = []
        for sn in snp_names:
            # Strip PLINK _N suffix: 'CTS1595_4' -> 'CTS1595'
            if sn.startswith('snp_24_'):
                base = '_'.join(sn.split('_')[:-1])
            else:
                base = sn.rsplit('_', 1)[0]
            snp_is_relevant.append(base in isogg_snp_names)

        n_relevant = sum(snp_is_relevant)
        print(f"  {len(snp_names)} Y-SNPs total, {n_relevant} ISOGG-relevant (kept)")
        print("Creating database...")
        conn = create_db(DB_FILE)
        cur = conn.cursor()

        # Insert tree and SNP reference 
        print("Inserting haplogroup tree...")
        cur.executemany('INSERT INTO haplogroup_tree VALUES (?,?)', tree)

        print("Inserting SNP reference...")
        cur.executemany('INSERT OR IGNORE INTO snp_reference VALUES (?,?,?,?,?)',
                        snp_rows)
        print(f"  {len(snp_rows)} SNP reference rows")
        conn.commit()

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


            # skip individuals with missing/uninformative haplogroup
            if is_invalid_haplo(y_haplo):
                skipped += 1
                continue

            # count non-missing genotypes (only ISOGG-relevant SNPs)
            geno_rows = []
            for snp_name, val, relevant in zip(snp_names, row[6:], snp_is_relevant):
                if relevant and val != 'NA':
                    geno_rows.append((iid, snp_name, int(val)))

            if not geno_rows:
                skipped += 1
                continue

            # parse lat/lon as floats
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

            y_clean = normalize_haplogroup(
                y_haplo,
                info.get('y_terminal'),
                snp_to_haplo,
                tree_by_name
            )
            if y_clean and y_clean in tree_nodes:
                norm_matched += 1

            cur.execute(
                'INSERT INTO individuals VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                (iid, y_haplo, y_clean, info.get('y_terminal'), info.get('group'),
                 info.get('locality'), info.get('country'),
                 lat, lon, date_mean, info.get('full_date'), len(geno_rows))
            )

            cur.executemany(
                'INSERT INTO genotypes VALUES (?,?,?)',
                geno_rows
            )

            inserted += 1
            if inserted % 200 == 0:
                print(f"  ... {inserted} individuals inserted")
                conn.commit()

    conn.commit()
    conn.close()
    print(f"\nDone! {inserted} individuals inserted, {skipped} skipped")
    print(f"  {skipped_non_ancient} skipped because they are not in Ancient_samples.txt")
    print(f"  {norm_matched}/{inserted} haplogroups resolved to tree nodes")
    if malformed_rows:
        print(f"Ignored {malformed_rows} malformed .raw rows")
    print(f"Database saved to {DB_FILE}")


if __name__ == '__main__':
    main()
