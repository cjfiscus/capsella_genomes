#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=std/test_env_%j.stdout
#SBATCH --error=std/test_env_%j.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=0-00:10:00
#SBATCH --job-name="aln_envtest"
#SBATCH -p koeniglab

# Smoke test for conda env syri_1.8.2: every tool must resolve inside the env.

source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate syri_1.8.2
set -euo pipefail
export PYTHONUNBUFFERED=1

echo "CONDA_PREFIX=$CONDA_PREFIX"
for t in minimap2 samtools seqkit syri plotsr python3; do
    p=$(which "$t")
    echo "$t -> $p"
    case "$p" in
        "$CONDA_PREFIX"/*) ;;
        *) echo "ERROR: $t resolves outside the env" >&2; exit 1 ;;
    esac
done

minimap2 --version
samtools --version | head -n1
seqkit version
syri --version
plotsr --version
python3 -c "import syri, pysam, numpy, pandas; print('python imports ok')"
seqkit sort --help | grep -- '--natural-order'
# every syri flag 002_align_syri.sh uses must exist (read --help, don't assume)
for f in -c -r -q -F --nc --dir --prefix; do syri --help | grep -q -- "$f" || { echo "ERROR: syri has no $f" >&2; exit 1; }; done; echo 'syri flags ok'
