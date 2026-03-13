# Parse user-uploaded genotype files for Y chromosome SNPs

import re
import sqlite3


def _split_row(raw_line):
    # Split a line by comma, tab, or whitespace
    raw = raw_line.strip()
    if not raw:
        return []
    if raw.startswith('#'):
        raw = raw[1:].strip()
    if ',' in raw:
        return [p.strip() for p in raw.split(',')]
    if '\t' in raw:
        return [p.strip() for p in raw.split('\t')]
    return re.split(r'\s+', raw)


def parse_user_file(file_content, db_path):
    # Parse a user genotype file and return a dict of SNPs
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT snp_name, position, ref_allele, alt_allele FROM snp_reference")

    ref_by_pos = {}
    ref_by_name = {}
    for name, pos, ref, alt in cur:
        ref_by_pos[str(pos)] = (name, ref)
        ref_by_name[name.upper()] = (name, ref, alt)
    conn.close()

    # Acceptable Y chromosome labels
    y_labels = {'24', 'Y', 'CHRY', 'CHR Y', 'CHR Y', 'Y_CHR', 'YCHR', 'Y-CHR', 'CHR Y', 'CHRY', 'Y', 'chrY', 'y', 'chrY'.lower(), 'chrY'.upper()}

    # Try to auto-detect header and column order
    header = None
    for line in file_content.splitlines():
        if not line.strip() or line.strip().startswith('#'):
            continue
        parts = _split_row(line)
        if len(parts) < 4:
            continue
        # Heuristic: header if any field contains 'chrom' or 'allele' or 'geno' or 'pos'
        if any(any(x in p.lower() for x in ['chrom', 'allele', 'geno', 'pos', 'snp', 'rsid', 'name']) for p in parts):
            header = [p.strip().lower() for p in parts]
            break

    # Default column indices
    idx_rsid = 0
    idx_chrom = 1
    idx_pos = 2
    idx_allele1 = 3
    idx_allele2 = 4
    idx_genotype = 3

    if header:
        # Try to find columns by name
        for i, col in enumerate(header):
            if 'chrom' in col:
                idx_chrom = i
            if 'pos' in col:
                idx_pos = i
            if 'allele1' in col or (col == 'allele' and idx_allele1 == 3):
                idx_allele1 = i
            if 'allele2' in col:
                idx_allele2 = i
            if 'geno' in col or 'result' in col:
                idx_genotype = i
            if 'rsid' in col or 'snp' in col or 'name' in col:
                idx_rsid = i

    result = {}
    for line in file_content.splitlines():
        if not line.strip() or line.strip().startswith('#'):
            continue
        parts = _split_row(line)
        if header and [p.strip().lower() for p in parts] == header:
            continue
        if len(parts) < 4:
            continue

        # Get chromosome label
        chrom = parts[idx_chrom].strip().upper() if len(parts) > idx_chrom else ''
        # Accept more Y chromosome labels
        if chrom not in y_labels:
            continue

        # Get position
        pos = parts[idx_pos].strip() if len(parts) > idx_pos else ''
        try:
            int(pos)
        except Exception:
            continue

        # Try to get alleles or genotype
        allele1 = allele2 = None
        if len(parts) > idx_allele2 and parts[idx_allele1] and parts[idx_allele2]:
            allele1 = parts[idx_allele1].strip().upper()
            allele2 = parts[idx_allele2].strip().upper()
        elif len(parts) > idx_genotype and parts[idx_genotype]:
            genotype = parts[idx_genotype].strip().upper()
            if len(genotype) >= 2:
                allele1, allele2 = genotype[0], genotype[1]
            elif len(genotype) == 1:
                allele1 = allele2 = genotype
        if not allele1:
            continue

        # Try to match by position or rsid
        if pos not in ref_by_pos:
            rsid = parts[idx_rsid].strip().upper() if len(parts) > idx_rsid else ''
            if rsid not in ref_by_name:
                continue
            snp_name, ref_allele, alt_allele = ref_by_name[rsid]
        else:
            snp_name, ref_allele, alt_allele = ref_by_name.get(
                ref_by_pos[pos][0].upper(),
                (ref_by_pos[pos][0], ref_by_pos[pos][1], None),
            )

        ref_allele = (ref_allele or '').upper()
        alt_allele = (alt_allele or '').upper()

        # Y is haploid, just use allele1
        allele = allele1
        if allele in ('', '-', '--', '0', 'N'):
            continue
        if allele == ref_allele:
            result[snp_name] = 0
        elif allele == alt_allele:
            result[snp_name] = 2
        else:
            pass

    return result
