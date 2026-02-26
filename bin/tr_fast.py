#!/usr/bin/env python3

# Usage: in the following example the words found in the second field
# (sys.argv[5]) of file.tsv (sys.argv[2]) that are specified in the first
# column (sys.argv[3]) of dict.tsv (sys.argv[1]) are translated with the
# corresponding word found in the third column (sys.argv[4]) of dict.tsv.
# Rows of file.tsv where no translation could occurr won't be discarded 
# (sys.argv[6]).
#
#   translator.py \
#       dict.tsv \
#       file.tsv \
#       1 \
#       3 \
#       2 \
#       1 \
#       > translated_file.tsv
#
# Specifying multiple fields on which to perform the translation is also
# possible, e.g. to translate on the first and the second field of file.tsv:
#
#   translator.py \
#       dict.tsv \
#       file.tsv \
#       1 \
#       3 \
#       1,2 \
#       1 \
#       > translated_file.tsv
#
# or you can also use 'all' to translate everywhere:
#
#   translator.py \
#       dict.tsv \
#       file.tsv \
#       1 \
#       3 \
#       all \
#       1 \
#       > translated_file.tsv
#

import sys
import csv


def parse_translate_on(arg: str, n_fields: int | None = None) -> list[int] | None:
    """
    Returns:
      None  -> translate on all fields
      list of 0-based column indices -> translate only those
    """
    if arg == "all":
        return None
    cols = [int(x) - 1 for x in arg.split(",") if x.strip() != ""]
    if n_fields is not None:
        cols = [c for c in cols if 0 <= c < n_fields]
    return cols


def load_mapping(dict_file: str, from_col_1based: int, to_col_1based: int) -> dict[str, str]:
    from_idx = from_col_1based - 1
    to_idx = to_col_1based - 1
    mapping: dict[str, str] = {}

    with open(dict_file, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            # skip empty / short rows
            if not row or max(from_idx, to_idx) >= len(row):
                continue
            src = row[from_idx]
            dst = row[to_idx]
            if src == "" or dst == "":
                continue
            # first match wins (matches your old behavior of taking values[0])
            if src not in mapping:
                mapping[src] = dst
    return mapping


def main() -> int:
    dict_file = sys.argv[1]
    target_file = sys.argv[2]
    from_col = int(sys.argv[3])
    to_col = int(sys.argv[4])
    translate_on_arg = sys.argv[5]
    keep_untranslated = bool(int(sys.argv[6]))

    mapping = load_mapping(dict_file, from_col, to_col)

    with open(target_file, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        out = sys.stdout

        for row in reader:
            if translate_on_arg == "all":
                indices = range(len(row))
            else:
                # translate only requested columns (1-based -> 0-based)
                indices = [int(x) - 1 for x in translate_on_arg.split(",")]

            translated_all = True
            for idx in indices:
                if idx < 0 or idx >= len(row):
                    translated_all = False
                    continue
                tr = mapping.get(row[idx])
                if tr is None:
                    translated_all = False
                else:
                    row[idx] = tr

            if translated_all or keep_untranslated:
                out.write("\t".join(row) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
