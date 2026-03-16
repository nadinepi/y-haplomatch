# Parse user-uploaded genotype files for Y chromosome SNPs

import re
import sqlite3


def _split_row(raw_line):
    # split a line by comma, tab, or whitespace
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
    # parse a user genotype file and return 2016 SNP names
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute(
        "SELECT snp_name, position, ref_allele, alt_allele FROM snp_reference ORDER BY position, snp_name"
    )
    snp_by_name = {}
    snp_by_pos = {}
    for name, pos, ref, alt in cur.fetchall():
        snp_by_name[name.upper()] = (name, pos, ref, alt)
        snp_by_pos.setdefault(str(pos), (name, pos, ref, alt))
    conn.close()

    # acceptable labels for y chromosome in the chromosome column
    y_labels = {
        '24', 'Y', 'CHRY', 'CHR Y', 'Y_CHR', 'YCHR', 'Y-CHR', 'CHR24', 'CHROMOSOMEY'
    }

    # try to auto-detect the header
    header = None
    for line in file_content.splitlines():
        if not line.strip() or line.strip().startswith('#'):
            continue
        parts = _split_row(line)
        if len(parts) < 4:
            continue
        if any(any(x in p.lower() for x in ['chrom', 'allele', 'geno', 'pos', 'snp', 'rsid', 'name']) for p in parts):
            header = [p.strip().lower() for p in parts]
            break

    # default column indices
    idx_rsid = 0
    idx_chrom = 1
    idx_pos = 2
    idx_allele1 = 3
    idx_allele2 = 4
    idx_genotype = 3

    if header:
        # try to find columns by name
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
    display_labels = {}
    log_lines = []
    for line in file_content.splitlines():
        if not line.strip() or line.strip().startswith('#'):
            continue
        parts = _split_row(line)
        if header and [p.strip().lower() for p in parts] == header:
            continue
        if len(parts) < 4:
            continue

        # get chromosome label
        chrom = parts[idx_chrom].strip().upper() if len(parts) > idx_chrom else ''
        if chrom not in y_labels:
            continue

        # get position
        pos = parts[idx_pos].strip() if len(parts) > idx_pos else ''
        try:
            int(pos)
        except Exception:
            continue

        # get alleles or genotype
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

        # match by 2016 SNP name first, then by position
        rsid = parts[idx_rsid].strip().upper() if len(parts) > idx_rsid else ''
        match_kind = None
        snp_info = None

        if rsid and rsid in snp_by_name:
            snp_info = snp_by_name[rsid]
            match_kind = 'DIRECT_2016_NAME'
        elif pos in snp_by_pos:
            snp_info = snp_by_pos[pos]
            match_kind = 'POS_MATCH_2016'
        else:
            log_lines.append(f"User SNP {rsid or pos}: no 2016 match in SNP table")
            continue

        snp_name, _, ref_allele, alt_allele = snp_info
        ref_allele = (ref_allele or '').upper()
        alt_allele = (alt_allele or '').upper()

        # Y is haploid, so just use the first allele
        allele = allele1
        if allele in ('', '-', '--', '0', 'N'):
            continue

        if allele == ref_allele:
            result[snp_name] = 0
            display_labels[snp_name] = rsid if rsid else snp_name
            log_lines.append(
                f"User SNP {snp_name}: match={match_kind}, allele={allele}, ref={ref_allele}, alt={alt_allele}, result=0 (ancestral)"
            )
        elif allele == alt_allele:
            result[snp_name] = 2
            display_labels[snp_name] = rsid if rsid else snp_name
            log_lines.append(
                f"User SNP {snp_name}: match={match_kind}, allele={allele}, ref={ref_allele}, alt={alt_allele}, result=2 (derived)"
            )
        else:
            log_lines.append(
                f"User SNP {snp_name}: match={match_kind}, allele={allele}, ref={ref_allele}, alt={alt_allele}, result=NA (no match)"
            )

    # write a small debug log for review
    with open("user_snp_debug.log", "w") as logf:
        for line in log_lines:
            logf.write(line + "\n")

    return result, display_labels
