#!/bin/bash
# Divide el FASTA consolidado en un archivo por secuencia (seqs/<seq_id>.fasta)
# y arma la lista para el array job. Corre en el nodo de login (es rapido).
set -e
mkdir -p seqs msa output
awk '/^>/{id=substr($1,2); gsub(/[^A-Za-z0-9_]/,"_",id); f="seqs/"id".fasta"} {print > f}' input/candidates.fasta
ls seqs/*.fasta > seq_list.txt
echo "Secuencias a procesar: $(wc -l < seq_list.txt)"
echo "Usa este numero como rango del array: -t 1-$(wc -l < seq_list.txt)"
