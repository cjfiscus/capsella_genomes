#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=std/summarize_%j.stdout
#SBATCH --error=std/summarize_%j.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=0-04:00:00
#SBATCH --job-name="aln_summ"
#SBATCH -p koeniglab

# Summarize the 30 SyRI runs (tables) and draw one chained multi-genome plotsr
# figure per homeolog set (Co / Cr).

# software dependencies
## syri 1.8.2, plotsr 1.2.0, python3 + pandas (conda env syri_1.8.2)

source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate syri_1.8.2
set -euo pipefail
export PYTHONUNBUFFERED=1

# SET VARIABLES
PAIRS=../data/pairs.tsv
SYRI=../results/02_syri
NORM=../results/01_norm
OUTDIR=../results/03_summary

#### PIPELINE #####
mkdir -pv "$OUTDIR"
echo "python: $(which python3)  plotsr: $(which plotsr)"

# 1) tables: sv_summary.tsv, pair_metrics.tsv, sv_size_filtered.tsv
python3 003a_summarize.py "$PAIRS" "$SYRI" "$NORM" "$OUTDIR"

# 2) one chained plotsr figure per set. Genome order = order of first appearance
#    in pairs.tsv (ref is always the earlier genome), and each adjacent pair
#    (genome i vs i+1) supplies one syri.out, so the chain is
#    G1 - G2 - G3 - ... exactly as in data/pairs.tsv.
for SET in $(cut -f1 "$PAIRS" | sort -u); do
    GENOMES=$(awk -F'\t' -v s="$SET" '$1==s { for (i=2;i<=3;i++) if (!($i in seen)) { seen[$i]=1; print $i } }' "$PAIRS")
    GFILE="$OUTDIR"/genomes_"$SET".txt
    printf '#file\tname\ttags\n' > "$GFILE"
    SR_ARGS=()
    PREV=""
    for G in $GENOMES; do
        printf '%s\t%s\tlw:1.5\n' "$(realpath "$NORM"/"$G".fa)" "$G" >> "$GFILE"
        if [ -n "$PREV" ]; then
            SRF="$SYRI"/"${PREV}_vs_${G}"/"${PREV}_vs_${G}.syri.out"
            [ -s "$SRF" ] || { echo "ERROR: missing $SRF (is ${PREV},${G} adjacent in pairs.tsv?)" >&2; exit 1; }
            SR_ARGS+=(--sr "$(realpath "$SRF")")
        fi
        PREV="$G"
    done
    echo "[$SET] genomes: $(echo $GENOMES | tr '\n' ' ')"
    plotsr "${SR_ARGS[@]}" --genomes "$GFILE" \
        -o "$OUTDIR"/"${SET}_chain.png" -H 16 -W 9 --lf "$OUTDIR"/plotsr_"${SET}".log
    echo "[$SET] wrote $OUTDIR/${SET}_chain.png"
done

echo "done: $OUTDIR"
ls -la "$OUTDIR"
