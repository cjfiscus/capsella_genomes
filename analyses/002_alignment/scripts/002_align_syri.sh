#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=std/syri_%a_%A.stdout
#SBATCH --error=std/syri_%a_%A.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=2-00:00:00
#SBATCH --job-name="aln_syri"
#SBATCH -p koeniglab
#SBATCH --array=1-30

# software dependencies
## minimap2 2.31, samtools, syri 1.8.2, plotsr 1.2.0 (conda env syri_1.8.2)

source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate syri_1.8.2
set -euo pipefail
export PYTHONUNBUFFERED=1

# SET VARIABLES
LIST=../data/pairs.tsv
NORM=../results/01_norm
OUTROOT=${OUTROOT:-../results/02_syri}     # overridable, e.g. for preset tests
PRESET=${PRESET:-asm10}            # one preset for all 30 pairs so SV counts stay comparable. asm10 chosen over asm20 from pilots on Co39_vs_Cbp2-2Co and Cr145_vs_CgMosher: more aligned, no loss on the most divergent pair (~0.9% syntenic SNP density), faster; SV counts (esp. translocations/dups) are preset-dependent
NCHR=8
TEMP_DIR="$TMPDIR"/syri_${SLURM_ARRAY_TASK_ID}   # $TMPDIR is the per-job scratch dir (see CLAUDE.md)

#### PIPELINE #####
mkdir -pv "$OUTROOT" "$TEMP_DIR"
trap 'rm -rf "$TEMP_DIR"' EXIT

# pairs.tsv: set  ref  qry
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
SET=$(echo "$LINE" | cut -f1)
REF=$(echo "$LINE" | cut -f2)
QRY=$(echo "$LINE" | cut -f3)
PAIR="${REF}_vs_${QRY}"
REF_FA=$(realpath "$NORM"/"$REF".fa)
QRY_FA=$(realpath "$NORM"/"$QRY".fa)

echo "[$PAIR] set=$SET ref=$REF_FA qry=$QRY_FA preset=$PRESET"
echo "[$PAIR] minimap2: $(which minimap2)  syri: $(which syri)  plotsr: $(which plotsr)"
for f in "$REF_FA" "$QRY_FA"; do
    [ -s "$f" ] || { echo "[$PAIR] ERROR: missing $f (run 001_chrom_normalize.sh)" >&2; exit 1; }
done

mkdir -pv "$OUTROOT"/"$PAIR"
OUT=$(realpath "$OUTROOT"/"$PAIR")   # absolute: syri's --dir chdirs, so relative input paths would break

# 1) align (BAM lives in $TMPDIR only — regenerable, and bigdata is nearly full)
BAM="$TEMP_DIR"/"$PAIR".bam
minimap2 -ax "$PRESET" --eqx -t "$SLURM_NTASKS" "$REF_FA" "$QRY_FA" \
    | samtools sort -@ 4 -o "$BAM" -
samtools index "$BAM"

# 2) SyRI
syri -c "$BAM" -r "$REF_FA" -q "$QRY_FA" -F B \
    --nc "$NCHR" --dir "$OUT" --prefix "${PAIR}."

if [ ! -s "$OUT"/"${PAIR}.syri.out" ] || [ ! -s "$OUT"/"${PAIR}.syri.summary" ]; then
    echo "[$PAIR] ERROR: syri output missing or empty in $OUT" >&2
    exit 1
fi

# 3) per-pair plotsr figure (cosmetic — warn, don't fail the pair)
printf '#file\tname\ttags\n%s\t%s\tlw:1.5\n%s\t%s\tlw:1.5\n' \
    "$REF_FA" "$REF" "$QRY_FA" "$QRY" > "$OUT"/genomes.txt
plotsr --sr "$OUT"/"${PAIR}.syri.out" --genomes "$OUT"/genomes.txt \
    -o "$OUT"/"${PAIR}.plotsr.png" -H 8 -W 6 --lf "$OUT"/plotsr.log \
    || echo "[$PAIR] WARN: plotsr failed (syri output is fine)" >&2

echo "[$PAIR] done: $OUT"
