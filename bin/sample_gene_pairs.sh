#!/usr/bin/env bash
# sample_gene_pairs.sh
# Sample 1000 unique unordered gene pairs from a list of unique gene names.
# Usage:
#   ./sample_gene_pairs.sh [-i genes.txt] [-n 1000] [-o output.tsv] [--seed 123]
#
# Notes:
# - Input is assumed to be one gene name per line (blank lines are ignored).
# - Pairs are printed as two columns separated by a TAB.
# - If not enough unique pairs exist (nC2 < target), exits with error.

set -euo pipefail

# Defaults
INPUT_FILE="genes.txt"
TARGET=1000
OUTPUT_FILE=""
SEED=""

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   INPUT_FILE="$2"; shift 2 ;;
    -n|--num)     TARGET="$2"; shift 2 ;;
    -o|--out)     OUTPUT_FILE="$2"; shift 2 ;;
    --seed)       SEED="$2"; shift 2 ;;
    -h|--help)
      sed -n '1,200p' "$0" | sed -n '1,40p'
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# Choose output destination
if [[ -n "$OUTPUT_FILE" ]]; then
  exec >"$OUTPUT_FILE"
fi

# AWK does the heavy lifting:
#  - Reads non-empty lines into an array genes[1..N]
#  - Checks N choose 2 >= TARGET, else errors
#  - Uniformly samples unordered pairs without replacement using a seen[] set
awk -v target="$TARGET" -v seed="$SEED" '
BEGIN {
  # Seed RNG
  if (length(seed) > 0) srand(seed + 0); else srand();
}
NF > 0 { genes[++N] = $0 }  # skip blank lines
END {
  if (N < 2) {
    print "ERROR: Need at least 2 genes (after ignoring blank lines)." > "/dev/stderr";
    exit 1;
  }
  total = (N * (N - 1)) / 2;
  if (total < target) {
    printf "ERROR: Only %d unique unordered pairs possible; need %d.\n", total, target > "/dev/stderr";
    exit 1;
  }

  count = 0;
  while (count < target) {
    i = int(rand() * N) + 1;
    j = int(rand() * N) + 1;
    if (i == j) continue;

    # Make unordered (i < j)
    if (i > j) { tmp = i; i = j; j = tmp; }

    key = i SUBSEP j;
    if (!(key in seen)) {
      seen[key] = 1;
      # Output two columns (TAB-separated)
      print genes[i] "\t" genes[j];
      count++;
    }
  }
}
' "$INPUT_FILE"
