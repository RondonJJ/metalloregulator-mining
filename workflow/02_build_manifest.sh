#!/bin/bash
# Construye el manifiesto que MetalNet2 espera en modo user_defined.
# Solo incluye secuencias que efectivamente tienen un a3m no vacio.
set -e
mkdir -p output
printf "seq_id\tmsa_file\n" > output/msa_files.tsv
n=0
for a in msa/*.a3m; do
    [ -s "$a" ] || continue
    id=$(basename "$a" .a3m)
    printf "%s\t%s\n" "$id" "$(readlink -f "$a")" >> output/msa_files.tsv
    n=$((n+1))
done
echo "MSAs incluidos en el manifiesto: $n"
