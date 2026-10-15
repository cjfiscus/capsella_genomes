#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=std/norm_%a_%A.stdout
#SBATCH --error=std/norm_%a_%A.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=0-04:00:00
#SBATCH --job-name="aln_norm"
#SBATCH -p koeniglab
#SBATCH --array=1-9

# Assign each chromosome to a subgenome (Cbp genomes) and record its homologous
# chromosome number + orientation in results/01_norm/<name>.chrmap.tsv, then write
# the per-subgenome FASTAs that syri consumes. Sequence NAMES ARE NEVER CHANGED;
# the only edit to sequence is reverse-complementing chromosomes that are
# '-' relative to the reference (recorded in chrmap.tsv), so syri doesn't call
# spurious whole-chromosome inversions.

# software dependencies
## minimap2 2.31, seqkit, samtools (conda env syri_1.8.2)

source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate syri_1.8.2
set -euo pipefail
export PYTHONUNBUFFERED=1

# SET VARIABLES
LIST=../data/genomes.tsv
FIXED=../data/subgenome_fixed.tsv
PREP=../results/00_prep
OUTDIR=${OUTDIR:-../results/01_norm}     # overridable for dry runs
SUBFILES=../data/subgenome_files.tsv
MIN_FRAC=0.7            # warn if < this fraction of a chromosome's aligned bases hit its best reference chromosome
NCHR_DIPLOID=8
TEMP_DIR="$TMPDIR"/norm_${SLURM_ARRAY_TASK_ID}   # $TMPDIR is the per-job scratch dir (see CLAUDE.md)

#### PIPELINE #####
mkdir -pv "$OUTDIR" "$TEMP_DIR"
trap 'rm -rf "$TEMP_DIR"' EXIT

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
NAME=$(echo "$LINE" | cut -f1)
N_CHR=$(echo "$LINE" | cut -f3)
GROUP=$(echo "$LINE" | cut -f4)     # Cbp | Co | Cr
MODE=$(echo "$LINE" | cut -f5)      # diploid | fixed | cbp  (see README)
IN_FA="$PREP"/"$NAME".fa

echo "[$NAME] group=$GROUP n_chrom=$N_CHR mode=$MODE which minimap2: $(which minimap2)"

# reference chromosome sets, named <Sub>_chrN
prefixed() {   # prefixed <fasta> <old-prefix-regex> <new-prefix>
    seqkit replace -p "$2" -r "$3" "$1"
}
check_ref() {  # check_ref <fasta>  -> must contain exactly Co/Cr_chr1..8 (whichever subs are present)
    local subs; subs=$(seqkit seq -n -i "$1" | sed 's/_chr.*//' | sort -u | paste -sd' ')
    local exp obs
    exp=$(for s in $subs; do for i in $(seq 1 $NCHR_DIPLOID); do echo "${s}_chr${i}"; done; done | sort | paste -sd,)
    obs=$(seqkit seq -n -i "$1" | sort | paste -sd,)
    if [ "$exp" != "$obs" ]; then echo "[$NAME] ERROR: reference chromosome names unexpected: $obs" >&2; exit 1; fi
}

# Cbp2-2 subgenome reference with original names: for the same-species passes.
# Renaming (to <Sub>_chrN, using Cbp2-2.chrmap.tsv) only happens in this
# temporary reference; it is not written to any output.
cbp2_ref() {   # cbp2_ref <sub> <out.fa>
    local sub=$1 out=$2 m="$OUTDIR"/Cbp2-2.chrmap.tsv
    [ -s "$m" ] || { echo "[$NAME] ERROR: $m missing — Cbp2-2 must be processed first" >&2; exit 1; }
    awk -F'\t' -v s="$sub" 'NR>1 && $2==s { print $1 }' "$m" > "$TEMP_DIR"/cbp2_$sub.ids
    awk -F'\t' -v OFS='\t' -v s="$sub" 'NR>1 && $2==s { print $1, s"_chr"$3 }' "$m" > "$TEMP_DIR"/cbp2_$sub.kv
    seqkit grep -f "$TEMP_DIR"/cbp2_$sub.ids "$PREP"/Cbp2-2.fa \
        | seqkit replace -p '^(\S+)$' -r '{kv}' -k "$TEMP_DIR"/cbp2_$sub.kv > "$out"
}

# vote: per query chromosome -> best reference chromosome (most matching bases),
# majority strand on it, and the fraction of its matching bases that went there.
# Output lines: qchrom  reference_chrom  strand  support_frac  aligned_bases
vote() {       # vote <paf>
    awk -F'\t' -v OFS='\t' '
        { m[$1,$6,$5]+=$10; t[$1,$6]+=$10; tot[$1]+=$10; q[$1]=1 }
        END {
            for (k in t) {
                split(k, a, SUBSEP)
                if (t[k] > best[a[1]]) { best[a[1]]=t[k]; bt[a[1]]=a[2] }
            }
            for (c in q) {
                strand = (m[c,bt[c],"+"] >= m[c,bt[c],"-"]) ? "+" : "-"
                printf "%s\t%s\t%s\t%.3f\t%d\n", c, bt[c], strand, best[c]/tot[c], tot[c]
            }
        }' "$1"
}

RAW="$TEMP_DIR"/chrmap.raw
: > "$RAW"

case "$MODE" in
    diploid|cbp)
        # one pass, all chromosomes vs a combined Co+Cr reference
        REF="$TEMP_DIR"/ref.fa
        if [ "$MODE" = "diploid" ]; then
            # Co39 / Cr145 name their chromosomes SCF_1..SCF_8
            prefixed "$PREP"/Co39.fa  '^SCF_' 'Co_chr'  > "$REF"
            prefixed "$PREP"/Cr145.fa '^SCF_' 'Cr_chr' >> "$REF"
        else
            # same-species reference: our Cbp2-2 subgenomes
            cbp2_ref Co "$TEMP_DIR"/ref_Co.fa; cbp2_ref Cr "$TEMP_DIR"/ref_Cr.fa
            cat "$TEMP_DIR"/ref_Co.fa "$TEMP_DIR"/ref_Cr.fa > "$REF"
        fi
        check_ref "$REF"
        minimap2 -x asm20 --secondary=no -t "$SLURM_NTASKS" "$REF" "$IN_FA" > "$TEMP_DIR"/map.paf
        vote "$TEMP_DIR"/map.paf >> "$RAW"
        ;;
    given)
        # Cbp2-2: subgenome membership is GIVEN by the existing separate
        # subgenome FASTAs (data/subgenome_files.tsv, read-only). No alignment,
        # nothing inferred. chr = rank of the sequence within its subgenome file
        # (natural order of names), strand '+', informational only.
        for s in Co Cr; do
            SRC=$(awk -F'\t' -v g="$NAME" -v s="$s" '$1==g && $2==s { print $3 }' "$SUBFILES")
            [ -s "$SRC" ] || { echo "[$NAME] ERROR: no $s subgenome file for $NAME in $SUBFILES" >&2; exit 1; }
            seqkit seq -n -i "$SRC" | sort -V | awk -F'\t' -v OFS='\t' -v s="$s" '{ printf "%s\t%s_chr%d\t+\tNA\tNA\n", $1, s, NR }' >> "$RAW"
        done
        # every chromosome of the genome must be assigned exactly once
        if [ "$(cut -f1 "$RAW" | sort | uniq -d | wc -l)" -ne 0 ]; then
            echo "[$NAME] ERROR: a sequence appears in both subgenome files" >&2; exit 1
        fi
        ;;
    fixed|fixed_cbp)
        # subgenome of each chromosome is GIVEN (data/subgenome_fixed.tsv) and
        # never inferred. Each chromosome is only compared against the reference
        # of its own subgenome, for its homolog number + orientation
        # (a combined-ref vote is unreliable for homeologous-exchange-affected
        # chromosomes, e.g. Cbp2-2 SCF_2 splits 50/50 Co/Cr).
        #   fixed     : reference = diploid Co39 (Co) / Cr145 (Cr)   [Cbp2-2]
        #   fixed_cbp : reference = Cbp2-2 subgenomes                [CbpChinese]
        for s in Co Cr; do
            awk -F'\t' -v g="$NAME" -v s="$s" '$1==g && $3==s { print $2 }' "$FIXED" > "$TEMP_DIR"/q_$s.txt
            [ -s "$TEMP_DIR"/q_$s.txt ] || { echo "[$NAME] ERROR: no $s chromosomes for $NAME in $FIXED" >&2; exit 1; }
            seqkit grep -f "$TEMP_DIR"/q_$s.txt "$IN_FA" > "$TEMP_DIR"/q_$s.fa
            if [ "$MODE" = "fixed" ]; then
                if [ "$s" = "Co" ]; then DIP=Co39; else DIP=Cr145; fi
                prefixed "$PREP"/"$DIP".fa '^SCF_' "${s}_chr" > "$TEMP_DIR"/ref_$s.fa
            else
                cbp2_ref "$s" "$TEMP_DIR"/ref_$s.fa
            fi
            check_ref "$TEMP_DIR"/ref_$s.fa
            minimap2 -x asm20 --secondary=no -t "$SLURM_NTASKS" "$TEMP_DIR"/ref_$s.fa "$TEMP_DIR"/q_$s.fa > "$TEMP_DIR"/map_$s.paf
            vote "$TEMP_DIR"/map_$s.paf >> "$RAW"
        done
        ;;
    *) echo "[$NAME] ERROR: unknown mode '$MODE'" >&2; exit 1 ;;
esac

N_MAPPED=$(wc -l < "$RAW")
if [ "$N_MAPPED" -ne "$N_CHR" ]; then
    echo "[$NAME] ERROR: only $N_MAPPED of $N_CHR chromosomes had alignments to the reference" >&2
    exit 1
fi

# chrmap.tsv: qchrom  subgenome  chr  strand  support_frac  aligned_bases
# (chr = homologous reference chromosome number, informational only; qchrom names are never changed)
MAP="$OUTDIR"/"$NAME".chrmap.tsv
printf 'qchrom\tsubgenome\tchr\tstrand\tsupport_frac\taligned_bases\n' > "$MAP"
awk -F'\t' -v OFS='\t' '{
    split($2, r, "_chr")          # Co_chr3 -> Co, 3
    print $1, r[1], r[2], $3, $4, $5
}' "$RAW" | sort -k2,2 -k3,3n >> "$MAP"
cat "$MAP"

# 4) require a clean 1:1 assignment
if [ "$GROUP" = "Cbp" ]; then SUBS="Co Cr"; else SUBS="$GROUP"; fi
EXPECT="$TEMP_DIR"/expected.txt
for s in $SUBS; do for i in $(seq 1 $NCHR_DIPLOID); do printf '%s\t%s\n' "$s" "$i"; done; done | sort > "$EXPECT"
tail -n +2 "$MAP" | cut -f2,3 | sort > "$TEMP_DIR"/observed.txt
if ! diff "$EXPECT" "$TEMP_DIR"/observed.txt; then
    echo "[$NAME] ERROR: chromosome assignment is not 1:1 onto expected {$SUBS} x chr1-$NCHR_DIPLOID (diff above: < expected, > observed)" >&2
    exit 1
fi
tail -n +2 "$MAP" | awk -F'\t' -v f="$MIN_FRAC" -v n="$NAME" \
    '$5 < f { printf "[%s] WARN: %s -> %s_chr%s support_frac=%s < %s\n", n, $1, $2, $3, $5, f > "/dev/stderr" }'

# 5) write per-subgenome fasta(s) (Cbp split into Co / Cr). Sequence names are
#    left exactly as in the source; '-' strand chromosomes are reverse-complemented.
for s in $SUBS; do
    if [ "$GROUP" = "Cbp" ]; then OUT_NAME="${NAME}${s}"; else OUT_NAME="$NAME"; fi
    OUT_FA="$OUTDIR"/"$OUT_NAME".fa

    tail -n +2 "$MAP" | awk -F'\t' -v s="$s" '$2==s { print $1 }'            > "$TEMP_DIR"/keep.txt
    tail -n +2 "$MAP" | awk -F'\t' -v s="$s" '$2==s && $4=="-" { print $1 }' > "$TEMP_DIR"/minus.txt

    seqkit grep -f "$TEMP_DIR"/keep.txt "$IN_FA" > "$TEMP_DIR"/sel.fa
    {
        if [ -s "$TEMP_DIR"/minus.txt ]; then
            seqkit grep -f "$TEMP_DIR"/minus.txt "$TEMP_DIR"/sel.fa | seqkit seq -r -p -t dna
            seqkit grep -v -f "$TEMP_DIR"/minus.txt "$TEMP_DIR"/sel.fa
        else
            cat "$TEMP_DIR"/sel.fa
        fi
    } | seqkit sort -n -N | seqkit seq -w 80 > "$OUT_FA"

    # sanity: exactly the selected sequence names (unchanged), same total length
    OBS_NAMES=$(seqkit seq -n -i "$OUT_FA" | sort | paste -sd,)
    EXP_NAMES=$(sort "$TEMP_DIR"/keep.txt | paste -sd,)
    if [ "$OBS_NAMES" != "$EXP_NAMES" ]; then
        echo "[$NAME] ERROR: $OUT_FA sequence names '$OBS_NAMES' != expected '$EXP_NAMES'" >&2
        exit 1
    fi
    LEN_SEL=$(seqkit stats -T "$TEMP_DIR"/sel.fa | tail -n1 | cut -f5)
    LEN_OUT=$(seqkit stats -T "$OUT_FA" | tail -n1 | cut -f5)
    if [ "$LEN_SEL" != "$LEN_OUT" ]; then
        echo "[$NAME] ERROR: total length changed ($LEN_SEL -> $LEN_OUT)" >&2
        exit 1
    fi
    samtools faidx "$OUT_FA"
    echo "[$NAME] wrote $OUT_FA ($LEN_OUT bp, minus-strand chromosomes flipped: $(wc -l < "$TEMP_DIR"/minus.txt))"
done
