"""
Recolecta los FASTA *_idseqsclustered de todas las familias y muestras en un
unico candidates.fasta, reescribiendo cada header a: muestra__familia__id_original
Genera ademas provenance.tsv para reconstruir el origen despues de MetalNet2.

NOTA: el separador es '__' (doble guion bajo), NO '|'. MetalNet2 arma comandos de
shell con estos nombres, y el '|' se interpreta como tuberia y rompe el paso de MSA.
"""
from pathlib import Path
import csv

SEP = "__"   # separador seguro para shell y sistema de archivos

suffix    = snakemake.params.suffix                # noqa: F821
out_fasta = snakemake.output.fasta                 # noqa: F821
out_prov  = snakemake.output.prov                  # noqa: F821
in_dirs   = snakemake.input.dirs                    # noqa: F821

Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]   # primer token del header
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


rows = []
n_seqs = 0
with open(out_fasta, "w") as out:
    for d in in_dirs:
        sample = Path(d).name
        # cada familia es un subdirectorio con su *_idseqsclustered.fasta
        for fa in sorted(Path(d).rglob(f"*{suffix}.fasta")):
            family = fa.parent.name
            for orig_id, seq in read_fasta(fa):
                new_id = f"{sample}{SEP}{family}{SEP}{orig_id}"
                out.write(f">{new_id}\n{seq}\n")
                rows.append((new_id, sample, family, orig_id))
                n_seqs += 1

with open(out_prov, "w", newline="") as ph:
    w = csv.writer(ph, delimiter="\t")
    w.writerow(["seq_id", "sample", "family", "orig_id"])
    w.writerows(rows)

print(f"[collect_candidates] {n_seqs} secuencias de {len(in_dirs)} muestra(s) -> {out_fasta}")
