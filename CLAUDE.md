# capsella_genomes

## Project

Novel genome assemblies for the allopolyploid *Capsella bursa-pastoris* (Cbp)
and its two progenitor species: *C. rubella* (Cr, a selfer descended from the
actual progenitor *C. grandiflora*, Cg) and *C. orientalis* (Co). This repo
holds the scripts and per-analysis results for comparative/population
genomic work on these genomes.

## Repo layout

Every analysis lives under `analyses/NNN_name/` (zero-padded, sequential),
each with its own:

```
analyses/NNN_name/
  README        # what the analysis does, how to run it
  data/         # small inputs specific to this analysis (lists, maps)
  scripts/      # sbatch scripts, numbered in run order (000_..., 001_...)
  results/      # outputs (gitignored — large/regeneratable)
```

- New work goes in the next number, never mixed into an existing analysis
  folder unless it's a direct continuation.
- When moving/renaming files, use `git mv` (or `git rm` + `git add`) so git
  records renames rather than delete+add.
- `results/` directories and `scripts/**/std/` (sbatch log dirs) are
  gitignored — regenerate, don't commit.
- Existing analyses: `000_snp_calling` (GATK SNP calling/filtering across
  Cbp/Co/Cr), `001_edta` (EDTA TE annotation + joint/pan TE library),
  `002_alignment` (pairwise minimap2 + SyRI across the 9 genomes, each Cbp
  split into its Co and Cr subgenomes; branch `002_alignment`).

## Working rules

- **Source genome files are read-only.** Never modify (edit, rename, trim,
  reverse-complement, re-index, make writable) anything under
  `/rhome/cfisc004/shared/CAPSELLA_GENOMES/` or
  `/rhome/cfisc004/shared/GENOMES/CAPSELLA/`. Process *copies* into an
  analysis's `results/` (renaming or subsetting there is fine), and write
  `.fai`/index files next to the copies, never next to the sources. Our own
  assemblies are `chmod a-w` and verifiable with `md5sum -c
  /rhome/cfisc004/shared/CAPSELLA_GENOMES/genomes/fa.md5`. User-level
  permission rules and a PreToolUse hook also block writes to these paths.
- **Subgenome membership that is already known is given, not inferred.**
  For Cbp2-2 it is the existing `Cbp2-2Co_…` / `Cbp2-2Cr_…` FASTAs
  (`SCF_1`-`SCF_8` Co, `SCF_9`-`SCF_16` Cr). Never re-derive or change it from
  alignments; alignments may only add informational columns (homolog
  number, strand).
- **Ask before changing how chromosomes/subgenomes are named, assigned or
  oriented** in derived data, even in copies — say what you plan to do first.
  Record any such mapping in a table next to the output (e.g. a `chrmap.tsv`).

## Genomes — this project's assemblies

All at `/rhome/cfisc004/shared/CAPSELLA_GENOMES/genomes/`, described further
in `/rhome/cfisc004/shared/CAPSELLA_GENOMES/README`.

| Name | File | Notes |
|---|---|---|
| Cbp full genome | `Cbp2-2_plus_manual_correction.fasta` | 77 seqs, 278 Mb. Contains both subgenomes: `SCF_1`-`SCF_8` = Co subgenome, `SCF_9`-`SCF_16` = Cr subgenome, remaining scaffolds unplaced. |
| Cbp Co subgenome | `Cbp2-2Co_plus_manual_correction.fasta` | `SCF_1`-`SCF_8` only, extracted from Cbp. |
| Cbp Cr subgenome | `Cbp2-2Cr_plus_manual_correction.fasta` | `SCF_9`-`SCF_16` only, extracted from Cbp. |
| Co genome | `Co39_plus_manual_correction.fasta` | 32 seqs, 136 Mb. |
| Cr genome | `Cr145_plus_manual_correction.fasta` | 129 seqs, 159 Mb. |
| Cr (Slotte et al. 2013) | `Crubella_474_v1.fa` | **External** — older published Cr assembly, kept alongside this project's own genomes for convenience but not one of our assemblies. |
| Neslia paniculata | `NP_1.0.fa` | **External** — outgroup assembly from Slotte et al. 2013, not one of our assemblies. |

Annotations (AUGUSTUS gene models) at
`/rhome/cfisc004/shared/CAPSELLA_GENOMES/annotation/`: `Cbp2-2.v1.gff`,
`Co39.gff`, `Cr145.gff`, plus matching `*_pep.fa` protein fastas.

Other shared data in `/rhome/cfisc004/shared/CAPSELLA_GENOMES/`:
`alignments/` (whole-genome alignments, hal/maf), `vcf/` (filtered SNP
calls per species/subgenome, from `000_snp_calling`).

## Genomes — outside/comparison assemblies

At `/rhome/cfisc004/shared/GENOMES/CAPSELLA/OUTSIDE/`, one subdirectory per
source. Owned by `dkoenig`; group-readable.

| Directory | Genome file | Annotation file | Notes |
|---|---|---|---|
| `CAPSELLA_BP_CHINESE` | `GCA_036452645.1_ASM3645264v1_genomic.fna` | `Capsella_burse-pastoris.genes.gff3` | GFF uses seqids `CbpA1`-`CbpA8`/`CbpB1`-`CbpB8`; FASTA uses GenBank accessions `CP144183.1`-`CP144198.1`. Mapping (verified by matching max GFF coordinate to `_assembly_report.txt` chromosome lengths): `CbpA1..A8` → `CP144183.1..CP144190.1`, `CbpB1..B8` → `CP144191.1..CP144198.1`. |
| `CAPSELLA_BP_RUSSIAN` | `GCA_001974645.2_ASM197464v2/GCA_001974645.2_ASM197464v2_genomic.fna` | `..._genomic.gff` | 1253 seqs. No annotation on NCBI for this assembly; the GFF here is from Kasianov et al. 2024 (BMC Biology, BRAKER2+Liftoff), remapped to GenBank accessions — see `ANNOTATION_SOURCE.txt` in that dir. |
| `CAPSELLA_GRANDIFLORA_MOSHER` | `Cgrandiflora_hapA_RAM_v1.1.fa.gz` | `Cgrandiflora_hapA_RAM_exon_v3.2.gff3.gz` | Seqids `lcl_Cg*`; FASTA headers carry extra `[location=...]` text after the ID. |
| `CAPSELLA_ORIENTALIS_MOSHER` | `Corientalis_RAM_v1.1.fa.gz` | `Corientalis_RAM_exon_v3.2.gff3.gz` | Seqids `lcl_Co_*`. |
| `CAPSELLA_ORIENTALIS_RUSSIAN` | `Corientalis_genome.fasta.gz` | `Corientalis_genome_annotation.gff.gz` | Seqids `O1`, `O2`, ... |
| `CAPSELLA_RUBELLA_MOSHER` | `Crubella_RAM_v1.1.fa.gz` | `Crubella_RAM_exon_v3.2.gff3.gz` | Seqids `lcl_Cr_*`. |

## Cluster / SLURM

Shared CentOS/Rocky HPC system. All Capsella jobs run on the
**`koeniglab`** partition: 8 nodes (`i[09-12]`, `r[17-20]`), 64+ CPUs and
~450 GB RAM per node, 30-day max walltime (7-day default). `/scratch` is
node-local (~1.5 TB); use it for job-local temp/intermediate files, not
final outputs. `/bigdata` is the real working storage (currently ~94%
full, 11 TB free — be mindful of large intermediate files).

**Scratch dir gotcha:** `/scratch/cfisc004` itself is root-owned — `mkdir
/scratch/cfisc004/<anything>` fails with `Permission denied`, even inside a
running job. SLURM auto-provisions one writable per-job dir at exactly
`/scratch/cfisc004/$SLURM_JOB_ID` (each array task gets its own unique ID,
so no collisions) and exposes it as `$TMPDIR`/`$SCRATCH`. Always use
`$TMPDIR` (or a subdir under it) for scratch space in scripts — never
hand-build a `/scratch/cfisc004/<name>_<id>` path.

sbatch header style used in this repo (see `analyses/000_snp_calling/`):

```bash
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=<n>
#SBATCH --mem-per-cpu=<mem>      # or --mem=<total>
#SBATCH --output=std/<name>%j.stdout
#SBATCH --error=std/<name>%j.stderr
#SBATCH --mail-user=fiscuscj@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=<d>-00:00:00
#SBATCH --job-name="<name>"
#SBATCH -p koeniglab

# software dependencies
## <tool> <version>

# SET VARIABLES
...
```

Array jobs read a list file by line number:
`NAME=$(head -n "$SLURM_ARRAY_TASK_ID" "$LIST" | tail -n1 | cut -f...)`.

**Conda activation ordering gotcha:** in any sbatch script that uses conda,
`conda activate <env>` must run *before* `set -euo pipefail`, never after.
Conda's activate shell functions aren't strict-mode-safe — under
`set -u`/`set -e` they can silently fail to finish updating `PATH` while
`$CONDA_PREFIX` still looks correct, so bare commands resolve to the base
conda env instead of the activated one. Template:

```bash
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate <env-name>
set -euo pipefail
export PYTHONUNBUFFERED=1   # avoid swallowed output from Python subprocesses
```

**`pipefail` + early-exiting reader gotcha:** under `set -o pipefail`, a
pipeline like `seqkit sort ... | seqkit head -n N` (or `... | head -n N`)
can fail the whole job, because `head` exits early and SIGPIPEs the upstream
command. Write the intermediate list to a file (in `$TMPDIR`) and take the
top-N from the file, or use `sed -n 1p` / `awk 'NR<=N'` on small outputs.
Likewise, read a tool's `--help` before scripting flags — don't assume
options like a `--tmpdir` exist.

Useful `module avail` tools (as alternative to conda envs):
`gffread/0.12.7`, `seqkit/2.4.0`, `samtools`, `agat` (0.8, 1.4.0),
`cufflinks/2.2.1`.

## Software / conda envs

All under `/rhome/cfisc004/bigdata/.conda/envs/`. Check `conda env list`
before creating a new one — this account accumulates stale envs.

- **`EDTA_2.3.0`** — current EDTA (TE annotation), used by `001_edta`.
  Bundles `panEDTA.sh` in its `bin/`. Several bioconda packaging bugs were
  found and fixed for this env (RepeatMasker/RMBlast config, FamDB,
  blast/makeblastdb binary clobber, a TIR-Learner pandas-3.0 patch) — full
  details in this account's memory,
  `~/.claude/projects/-rhome-cfisc004/memory/edta_2.3.0_bioconda_fixes.md`.
  Working reference test script: `~/bigdata/software/EDTA_2.3.0_test/`.
  EDTA 2.3.0 auto-encodes/decodes overly-long sequence IDs itself
  (`EDTA.pl`, seqid codec step) — genome headers just need trimming to the
  first whitespace-delimited token, not manual shortening.
  Older `EDTA` (2.0.0) and `edta` envs are stale — don't use.
- **`syri_1.8.2`** — SyRI 1.8.2 + plotsr 1.2.0 + minimap2 2.31 + samtools +
  seqkit, used by `002_alignment` (smoke test:
  `analyses/002_alignment/scripts/test_env.sh`). Notes: syri has no
  `--tmpdir`, and `--dir` chdirs, so pass absolute input paths.
- Other envs of note: `samtools`, `seqkit`, `repeatmasker`, `BUSCO`,
  `cactus`, `anchorwave`, `orthofinder`, `gff3toolkit`, `vep`, `structure`.

## Contact / attribution

Git user: `cjfiscus`. SLURM job notifications: `fiscuscj@gmail.com`.
