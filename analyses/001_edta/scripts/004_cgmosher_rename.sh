#!/bin/bash
# Make a renamed COPY of CgMosher (results/00_prep/CgMosher.fa is itself a
# header-trimmed copy; the source in /rhome/.../GENOMES is never touched).
# Why: RepeatMasker rejects IDs >50 chars; EDTA appends ":start..end#class",
# so 40-char scaffold names overflowed in `EDTA.pl --step final`.
# Rename: lcl_Scaffold_<N>_1_contigs_length_<L> -> lcl_Scaffold_<N>
# (<N> unique; "_1_" constant; length is redundant). lcl_Cg*.2/_ch/_mt unchanged.
# Mapping recorded in chrmap.tsv. Run directly (light): bash 004_cgmosher_rename.sh
set -euo pipefail

SRC=$(realpath ../results/00_prep/CgMosher.fa)
OUT=../results/03_cgmosher_renamed
mkdir -pv "$OUT"
OUT=$(realpath "$OUT")

awk -F'\t' -v OFS='\t' -v map="$OUT/chrmap.tsv" '
    /^>/ {
        old = substr($0, 2); new = old
        if (match(old, /^lcl_Scaffold_[0-9]+_1_contigs_length_[0-9]+$/)) {
            split(old, a, "_"); new = "lcl_Scaffold_" a[3]
        }
        if (new in seen) { print "ERROR: duplicate new id " new > "/dev/stderr"; exit 1 }
        seen[new] = 1
        print old, new > map
        print ">" new; next
    }
    { print }
' "$SRC" > "$OUT/CgMosher.fa"
sed -i '1i old_id\tnew_id' "$OUT/chrmap.tsv"

# checks: same count, same sequences, max ID length
N_OLD=$(grep -c '^>' "$SRC"); N_NEW=$(grep -c '^>' "$OUT/CgMosher.fa")
[ "$N_OLD" = "$N_NEW" ] || { echo "ERROR: count $N_OLD vs $N_NEW" >&2; exit 1; }
diff <(grep -v '^>' "$SRC" | md5sum) <(grep -v '^>' "$OUT/CgMosher.fa" | md5sum) \
    && echo "sequence bytes identical"
echo "sequences: $N_NEW; max ID length: $(grep '^>' "$OUT/CgMosher.fa" | awk '{l=length($0)-1; if(l>m)m=l} END{print m}')"
