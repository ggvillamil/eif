#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------
# Generate RNA-seq–filtered GTF and FASTA files
# ------------------------------------------------------------------------
# Usage:
#   ./filter_and_combine_transcriptome.sh REMOVE_TX GTF_IN FA_IN YEAST_GTF YEAST_FA OUT_PREFIX
# Example:
#   ./filter_and_combine_transcriptome.sh \
#     results/post/remove_tx.txt \
#     resources/human.gtf \
#     resources/human.fa \
#     resources/yeast.gtf \
#     resources/yeast.fa \
#     resources/filtered_output
# ------------------------------------------------------------------------

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 REMOVE_TX GTF_IN FA_IN YEAST_GTF YEAST_FA OUT_PREFIX"
  exit 1
fi

# Assign arguments
REMOVE_TX="$1"
GTF_IN="$2"
FA_IN="$3"
YEAST_GTF="$4"
YEAST_FA="$5"
OUT_PREFIX="$6"

# Derived output paths
GTF_FILTERED="${OUT_PREFIX}.human.filtered.gtf"
FA_FILTERED="${OUT_PREFIX}.human.filtered.fa"
GTF_COMBINED="${OUT_PREFIX}.combined.gtf"
FA_COMBINED="${OUT_PREFIX}.combined.fa"

echo "Filtering GTF..."
grep -vF -f "$REMOVE_TX" "$GTF_IN" > "$GTF_FILTERED"
echo "Filtered GTF written to $GTF_FILTERED"

echo "Filtering FASTA..."
awk 'NR==FNR {remove[$1]; next} /^>/ {keep=!($2 in remove)} keep' "$REMOVE_TX" "$FA_IN" > "$FA_FILTERED"
echo "Filtered FASTA written to $FA_FILTERED"

echo "Combining with yeast GTF and FASTA..."
cat "$GTF_FILTERED" "$YEAST_GTF" > "$GTF_COMBINED"
cat "$FA_FILTERED" "$YEAST_FA" > "$FA_COMBINED"
echo "Combined GTF written to $GTF_COMBINED"
echo "Combined FASTA written to $FA_COMBINED"

echo "Done."
