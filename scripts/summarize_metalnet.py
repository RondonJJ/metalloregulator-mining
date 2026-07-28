#!/usr/bin/env python3
"""
summarize_metalnet.py

Resume la salida cruda de MetalNet2 (pred_pairs.tsv) en un reporte por proteina:
cuales tienen un sitio de union a metal de alta confianza (filter_by_graph==1),
que residuos lo forman (numeracion 1-based, la misma del .gv y de la proteina),
la composicion CHED y las probabilidades. Une con membership.tsv para atribuir
cada proteina a su organismo y familia de BITACORA.

Sin dependencias externas.

Uso:
  python summarize_metalnet.py -p output/pred_pairs.tsv -m membership.tsv \
      -o output/site_summary.tsv --pairs-out output/confident_pairs.tsv
"""
import argparse
from collections import defaultdict


def load_membership(path):
    mem = {}
    if not path:
        return mem
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            mem[r[idx["seq_id"]]] = (r[idx.get("organisms", -1)] if "organisms" in idx else "",
                                     r[idx.get("families", -1)] if "families" in idx else "")
    return mem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--pred_pairs", required=True)
    ap.add_argument("-m", "--membership", default=None)
    ap.add_argument("-o", "--out", default="site_summary.tsv")
    ap.add_argument("--pairs-out", default=None,
                    help="(opcional) TSV con los pares de alta confianza en detalle")
    args = ap.parse_args()

    mem = load_membership(args.membership)

    # leer pred_pairs
    all_ids = []
    conf_edges = defaultdict(list)   # seq_id -> [(r1,pos1,r2,pos2,prob)]
    with open(args.pred_pairs) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if not r or len(r) < len(header):
                continue
            sid = r[ix["seq_id"]]
            if sid not in all_ids:
                all_ids.append(sid)
            if r[ix["filter_by_graph"]] == "1":
                p1 = int(r[ix["resi_seq_posi_1"]]) + 1     # -> 1-based
                p2 = int(r[ix["resi_seq_posi_2"]]) + 1
                conf_edges[sid].append(
                    (r[ix["resi_1"]], p1, r[ix["resi_2"]], p2,
                     float(r[ix["prob"]])))

    # armar filas de resumen
    rows = []
    for sid in all_ids:
        edges = conf_edges.get(sid, [])
        residues = {}
        probs = []
        for a, pa, b, pb, pr in edges:
            residues[pa] = a
            residues[pb] = b
            probs.append(pr)
        res_sorted = sorted(residues.items())
        res_str = ",".join("%s%d" % (aa, pos) for pos, aa in res_sorted)
        comp = defaultdict(int)
        for _, aa in res_sorted:
            comp[aa] += 1
        comp_str = ",".join("%d%s" % (comp[a], a) for a in "CHED" if comp[a])
        org, fam = mem.get(sid, ("", ""))
        rows.append({
            "seq_id": sid,
            "has_site": "yes" if edges else "no",
            "n_pairs": len(edges),
            "n_residues": len(residues),
            "site_residues": res_str if res_str else "-",
            "composition": comp_str if comp_str else "-",
            "max_prob": round(max(probs), 3) if probs else 0.0,
            "mean_prob": round(sum(probs) / len(probs), 3) if probs else 0.0,
            "organisms": org, "families": fam,
        })

    # ordenar: primero las que tienen sitio, mas residuos, mayor prob
    rows.sort(key=lambda x: (x["has_site"] == "yes", x["n_residues"],
                             x["max_prob"]), reverse=True)

    cols = ["seq_id", "has_site", "n_pairs", "n_residues", "site_residues",
            "composition", "max_prob", "mean_prob", "organisms", "families"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    if args.pairs_out:
        with open(args.pairs_out, "w") as fh:
            fh.write("seq_id\tresi_1\tpos_1\tresi_2\tpos_2\tprob\n")
            for sid in all_ids:
                for a, pa, b, pb, pr in conf_edges.get(sid, []):
                    fh.write("%s\t%s\t%d\t%s\t%d\t%.3f\n"
                             % (sid, a, pa, b, pb, pr))

    n_site = sum(1 for r in rows if r["has_site"] == "yes")
    print("Proteinas analizadas:        %d" % len(rows))
    print("Con sitio de alta confianza: %d" % n_site)
    print("Sin sitio predicho:          %d" % (len(rows) - n_site))
    print("Resumen: %s" % args.out)
    if args.pairs_out:
        print("Pares detallados: %s" % args.pairs_out)


if __name__ == "__main__":
    main()
