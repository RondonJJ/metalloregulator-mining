#!/usr/bin/env python3
"""
prepare_metalnet_input.py

Consolida los FASTA de candidatas de BITACORA (organizados en
<root>/<organismo>/<familia>/*.fasta) en una unica entrada lista para MetalNet2:

  - candidates.fasta : una copia por SECUENCIA UNICA, con seq_id limpio y unico
                       (apto como nombre de archivo, ya que MetalNet nombra sus
                        salidas segun el seq_id).
  - membership.tsv   : mapea cada seq_id -> accesion original, longitud, %X y
                       TODAS las combinaciones (organismo, familia) de las que
                       proviene. Sirve para re-atribuir los resultados despues.

Sin dependencias externas (solo Python 3 estandar).

Uso:
  python prepare_metalnet_input.py -i Results_Bita_-5 -o metalnet_input
"""
import argparse, glob, os, re
from collections import defaultdict


def parse_fasta(fp):
    h, s = None, []
    for line in open(fp):
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if h is not None:
                yield h, "".join(s)
            h, s = line[1:].strip(), []
        else:
            s.append(line.strip())
    if h is not None:
        yield h, "".join(s)


def safe_id(acc):
    """WP_003877217.1 -> WP_003877217_1 (seguro y unico como nombre de archivo)."""
    s = re.sub(r"[^A-Za-z0-9_]", "_", acc)
    return re.sub(r"_+", "_", s).strip("_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_root", required=True,
                    help="Directorio raiz de BITACORA (con <org>/<fam>/*.fasta)")
    ap.add_argument("-o", "--outdir", default="metalnet_input")
    ap.add_argument("--glob", default="*/*/*.fasta",
                    help="Patron relativo a input_root para hallar los FASTA")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.input_root, args.glob)))
    # seq_content -> dict(accs=set, members=set[(org,fam)])
    by_seq = {}
    order = []
    n_records = 0
    for f in files:
        parts = f.split(os.sep)
        org, fam = parts[-3], parts[-2]
        for h, s in parse_fasta(f):
            n_records += 1
            acc = h.split()[0]
            seq = s.rstrip("*").upper()
            if not seq:
                continue
            if seq not in by_seq:
                by_seq[seq] = {"accs": [], "members": set()}
                order.append(seq)
            if acc not in by_seq[seq]["accs"]:
                by_seq[seq]["accs"].append(acc)
            by_seq[seq]["members"].add((org, fam))

    # asignar seq_id unico
    used = set()
    rows = []
    fasta_path = os.path.join(args.outdir, "candidates.fasta")
    with open(fasta_path, "w") as fh:
        for seq in order:
            info = by_seq[seq]
            base = safe_id(info["accs"][0])
            sid, i = base, 2
            while sid in used:                 # colision: misma accesion, otra secuencia
                sid = "%s_v%d" % (base, i); i += 1
            used.add(sid)
            fh.write(">%s\n" % sid)
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j + 60] + "\n")
            xf = seq.count("X") / len(seq)
            orgs = sorted({o for o, _ in info["members"]})
            fams = sorted({fa for _, fa in info["members"]})
            rows.append([sid, ";".join(info["accs"]), len(seq),
                         round(xf * 100, 2), len(orgs), len(fams),
                         ";".join(orgs), ";".join(fams)])

    mem_path = os.path.join(args.outdir, "membership.tsv")
    with open(mem_path, "w") as fh:
        fh.write("seq_id\taccessions\tlength\tpct_X\tn_org\tn_fam\t"
                 "organisms\tfamilies\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print("=== Consolidacion BITACORA -> MetalNet2 ===")
    print("FASTA de entrada procesados: %d" % len(files))
    print("Registros totales:           %d" % n_records)
    print("Secuencias unicas:           %d" % len(order))
    print("-" * 44)
    print("FASTA consolidado: %s" % fasta_path)
    print("Tabla de mapeo:    %s" % mem_path)


if __name__ == "__main__":
    main()
