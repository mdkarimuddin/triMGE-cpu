#!/bin/bash -l
#SBATCH --job-name=trimge_cpu_wsbs
#SBATCH --account=project_xxxxxxx
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/project_xxxxxxx/triMGE-cpu/logs/wsbs_cpu_%j.out
#SBATCH --error=/scratch/project_xxxxxxx/triMGE-cpu/logs/wsbs_cpu_%j.err
#
# Classify the WSBS parent contigs with triMGE-cpu.

set -euo pipefail

EXP=/scratch/project_xxxxxxx/nucleotide_15studies/trimge_wsbs_exploratory
CPU_REPO=/scratch/project_xxxxxxx/triMGE-cpu
INPUT_FASTA="${INPUT_FASTA:-${EXP}/input/wsbs_parent_contigs.fna}"
WORKDIR="${WORKDIR:-${EXP}/cpu_pipeline_out}"
MODELS="${MODELS:-${CPU_REPO}/models}"

mkdir -p "${CPU_REPO}/logs" "${WORKDIR}"

module purge 2>/dev/null || true
module load python-data

export PYTHONPATH="${CPU_REPO}/src:${PYTHONPATH:-}"
cd "${CPU_REPO}"

echo "=== triMGE-cpu ==="
echo "INPUT_FASTA=${INPUT_FASTA}"
echo "WORKDIR=${WORKDIR}"
echo "MODELS=${MODELS}"
echo "Host=$(hostname)  Date=$(date)"
grep -c '^>' "${INPUT_FASTA}" || true

python -m trimgecpu.easy_cli \
  --fasta "${INPUT_FASTA}" \
  --workdir "${WORKDIR}" \
  --models-dir "${MODELS}" \
  --verbose

echo "DONE $(date)"
echo "Predictions: ${WORKDIR}/02_predictions/predictions.tsv"
echo "Report:      ${WORKDIR}/03_report/"
