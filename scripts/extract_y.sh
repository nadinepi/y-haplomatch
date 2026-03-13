#!/bin/bash
# Extract Y chromosome (chr 24) data from AADR PLINK files using plink,
# then convert to a recoded text format we can parse in Python.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

AADR_DIR="$PROJECT_DIR/ProjetcFiles/AADR_54.1"
OUT_DIR="$PROJECT_DIR/data/plink_out"
mkdir -p "$OUT_DIR"

echo "Extracting chr 24 (Y) from AADR..."
plink --bfile "$AADR_DIR/v54.1_1240K_public" \
      --chr 24 \
      --filter-males \
      --recodeA \
      --out "$OUT_DIR/aadr_chrY"

echo ""
echo "Done. Output files in $OUT_DIR/"
ls -lh "$OUT_DIR/"
