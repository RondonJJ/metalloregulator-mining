#!/bin/bash
#$ -cwd
#$ -V
#$ -q h14.q
#$ -pe ompi511h14 4
#$ -N msa_bita
#$ -t 1-161
# ^^^ AJUSTA 161 al numero que imprime 00_split.sh (wc -l seq_list.txt)

# ===== SETEAR: base de datos de HHblits (UniRef30/BFD) =====
DB=/RUTA/A/TU/UniRef30_2023_02/UniRef30_2023_02
# ==========================================================
HHBLITS=hhblits            # o ruta absoluta al binario
CPU=4

f=$(sed -n "${SGE_TASK_ID}p" seq_list.txt)
id=$(basename "$f" .fasta)

# si ya existe un a3m no vacio, saltar (permite reintentos del array)
if [ -s "msa/${id}.a3m" ]; then
    echo "[$id] ya existe, salto"; exit 0
fi

$HHBLITS -i "$f" -oa3m "msa/${id}.a3m" -d "$DB" \
    -n 3 -e 1e-3 -cpu $CPU -maxfilt 100000 -diff inf -id 99 -cov 50
