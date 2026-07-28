#!/bin/bash
#$ -cwd
#$ -V
#$ -q h14.q
#$ -pe ompi511h14 4
#$ -N metalnet_bita

PYTHON=/users-d2/j.rondon/miniconda3/envs/metalnet/bin/python
SCRIPTS_FILE=/users-d2/j.rondon/programas/MetalNet2/model/scripts/run_prediction_workflow.py

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

mkdir -p ./output

$PYTHON $SCRIPTS_FILE \
    --input_fasta ./input/candidates.fasta \
    --input_msa ./output/msa_files.tsv \
    --msa_source user_defined \
    --output_pairs ./output/pred_pairs.tsv \
    --to_image \
    --keep_inter_files
