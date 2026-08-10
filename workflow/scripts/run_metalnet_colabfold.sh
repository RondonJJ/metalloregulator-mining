#!/bin/bash
# =====================================================================
# MetalNet2 en modo ColabFold, parametrizado.
# Basado en el script que ya te funciono. Se quitaron las directivas SGE
# (#$ ...) porque Snakemake envia el trabajo a la cola via el perfil.
#
# Uso:
#   run_metalnet_colabfold.sh PYTHON SCRIPT INPUT_FASTA OUTPUT_PAIRS MSA_SOURCE THREADS USE_GPU
# =====================================================================
set -euo pipefail

PYTHON="$1"          # ej: /users-d2/j.rondon/miniconda3/envs/metalnet/bin/python
SCRIPT="$2"          # ej: .../MetalNet2/model/scripts/run_prediction_workflow.py
INPUT_FASTA="$3"
OUTPUT_PAIRS="$4"
MSA_SOURCE="$5"      # colabfold
THREADS="$6"
USE_GPU="$7"         # 0 o 1

export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"

mkdir -p "$(dirname "$OUTPUT_PAIRS")"

CUDA_ARG=()
if [ "$USE_GPU" = "1" ]; then
    CUDA_ARG=(--cuda 1)   # acelera ESM2 si el nodo tiene GPU
fi

"$PYTHON" "$SCRIPT" \
    --input_fasta "$INPUT_FASTA" \
    --msa_source "$MSA_SOURCE" \
    --output_pairs "$OUTPUT_PAIRS" \
    --to_image \
    --keep_inter_files \
    "${CUDA_ARG[@]}"
