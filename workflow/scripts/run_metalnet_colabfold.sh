#!/bin/bash
# =====================================================================
# MetalNet2 en modo ColabFold, parametrizado.
# Basado en el script que ya te funciono. Se quitaron las directivas SGE
# (las maneja el wrapper de qsub / Snakemake).
#
# IMPORTANTE: activa el env 'metalnet' ANTES de correr, para que los
# subprocesos que MetalNet lanza con os.system() (p.ej. search_msa.py)
# tambien hereden ese entorno. Sin esto, el subproceso usa otro Python
# y falla con "No module named 'absl'".
#
# Uso:
#   run_metalnet_colabfold.sh PYTHON SCRIPT INPUT_FASTA OUTPUT_PAIRS MSA_SOURCE THREADS USE_GPU
# =====================================================================
set -euo pipefail

# ---- Activar el entorno de MetalNet ----
# Ajusta la ruta si tu miniconda esta en otro lugar.
source ~/miniconda3/etc/profile.d/conda.sh
conda activate metalnet

PYTHON="$1"          # (informativo) ruta al python de metalnet
SCRIPT="$2"          # .../MetalNet2/model/scripts/run_prediction_workflow.py
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

# Con el env ya activo, usamos 'python' (el de metalnet) en vez de la ruta absoluta,
# asi el proceso principal y sus subprocesos comparten el mismo interprete.
python "$SCRIPT" \
    --input_fasta "$INPUT_FASTA" \
    --msa_source "$MSA_SOURCE" \
    --output_pairs "$OUTPUT_PAIRS" \
    --to_image \
    --keep_inter_files \
    "${CUDA_ARG[@]}"
