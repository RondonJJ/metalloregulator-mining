#!/bin/bash
#$ -cwd
#$ -V
#$ -q h14.q
#$ -pe ompi511h14 4
#$ -N metalnet_pilot

# Prueba con pocas secuencias para confirmar que la API de ColabFold responde
# desde el nodo ANTES de lanzar las 161. Genera input/pilot.fasta si no existe.

PYTHON=/users-d2/j.rondon/miniconda3/envs/metalnet/bin/python
SCRIPTS_FILE=/users-d2/j.rondon/programas/MetalNet2/model/scripts/run_prediction_workflow.py
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4

mkdir -p ./output_pilot
# toma las 3 primeras secuencias del FASTA consolidado
awk '/^>/{n++} n<=3{print} n>3{exit}' ./input/candidates.fasta > ./input/pilot.fasta

$PYTHON $SCRIPTS_FILE \
    --input_fasta ./input/pilot.fasta \
    --msa_source colabfold \
    --output_pairs ./output_pilot/pred_pairs.tsv \
    --to_image \
    --keep_inter_files
