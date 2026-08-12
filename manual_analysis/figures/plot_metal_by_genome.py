#!/usr/bin/env python3
"""
plot_metal_by_genome.py

Genera la figura "cobertura potencial de metales por genoma" (SVG editable):
filas = metales (agrupados), columnas = genomas, y dentro de cada celda un
"chip" por cada familia presente en ese genoma capaz de sensar ese metal
(relleno oscuro = sitio CHED confirmado; claro = candidato sin sitio).

Los metales por familia son los inductores POTENCIALES descritos en la
literatura; la especificidad no está resuelta por homología. El estado por
(familia, genoma) se puede derivar de results/full_run/counts_by_genome_family.tsv
(site>0 -> 'site'; total>0 y site==0 -> 'cand'; total==0 -> 'abs'); aquí se deja
explícito para que la figura sea autocontenida y reproducible.

Uso:
  python plot_metal_by_genome.py -o metal_by_genome.svg
"""
import argparse

FAMMETALS = {
    "ArsR": ["Co","Ni","Cu","Zn","Cd","Hg","Pb","Bi","As","Sb"],
    "MerR": ["Fe","Co","Ni","Cu","Zn","Ag","Cd","Au","Hg","Pb"],
    "Fur":  ["Mn","Fe","Zn","Ni","Co"], "DtxR": ["Mn","Fe"], "Rrf2": ["Fe"],
    "CsoR": ["Co","Ni","Cu","Zn"], "CopY": ["Cu","Zn"], "MarR": ["Cu","Zn"],
    "GntR": ["Cu","Zn"], "LysR": ["Mo","W"], "TetR": ["Fe","Cu","Zn","As"],
}
ABBR = {"ArsR":"Ars","MerR":"Mer","Fur":"Fur","DtxR":"Dtx","Rrf2":"Rrf","CsoR":"Cso",
        "CopY":"Cop","MarR":"Mar","GntR":"Gnt","LysR":"Lys","TetR":"Tet"}
FAMORDER = ["ArsR","MerR","Fur","DtxR","Rrf2","CsoR","CopY","MarR","GntR","LysR","TetR"]
# estado por (familia, genoma): site | cand | abs
STATUS = {
 "ArsR":{"Ma":"site","Mt":"site","Vc":"cand","Vv":"cand"},
 "MerR":{"Ma":"site","Mt":"cand","Vc":"site","Vv":"site"},
 "Fur":{"Ma":"site","Mt":"site","Vc":"site","Vv":"site"},
 "DtxR":{"Ma":"site","Mt":"site","Vc":"abs","Vv":"abs"},
 "Rrf2":{"Ma":"cand","Mt":"cand","Vc":"site","Vv":"site"},
 "CsoR":{"Ma":"site","Mt":"site","Vc":"abs","Vv":"abs"},
 "CopY":{"Ma":"cand","Mt":"cand","Vc":"abs","Vv":"abs"},
 "MarR":{"Ma":"cand","Mt":"cand","Vc":"cand","Vv":"site"},
 "GntR":{"Ma":"site","Mt":"site","Vc":"site","Vv":"site"},
 "LysR":{"Ma":"cand","Mt":"abs","Vc":"abs","Vv":"abs"},
 "TetR":{"Ma":"cand","Mt":"cand","Vc":"abs","Vv":"cand"},
}
GROUPS = [("Transición esenciales",["Fe","Mn","Zn","Ni","Co","Cu"]),
          ("Metales pesados tóxicos",["Cd","Pb","Hg","Ag","Au"]),
          ("Metaloides",["As","Sb","Bi"]), ("Oxianiones",["Mo","W"])]
GENOMES = [("Ma","M. avium"),("Mt","M. tuberculosis"),("Vc","V. cholerae"),("Vv","V. vulnificus")]


def build_svg():
    metal2fam = {}
    for f in FAMORDER:
        for m in FAMMETALS[f]:
            metal2fam.setdefault(m, []).append(f)
    Lx=150; Gx=Lx; cw=104; rh=34; y0=58
    nrows=sum(len(ms) for _,ms in GROUPS)
    W=Gx+cw*4+12; H=y0+nrows*rh+92
    s=[f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg">',
       '<style>text{font-family:Arial}</style>',
       f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>',
       f'<text x="22" y="26" font-size="13" font-weight="500" fill="#1a1a1a">Cobertura potencial de metales por genoma</text>']
    for gi,(gk,gl) in enumerate(GENOMES):
        cx=Gx+gi*cw+cw/2
        s.append(f'<text x="{cx:.0f}" y="{y0-8}" text-anchor="middle" font-size="11.5" font-style="italic" fill="#1a1a1a">{gl}</text>')
        s.append(f'<line x1="{Gx+gi*cw}" y1="{y0-4}" x2="{Gx+gi*cw}" y2="{y0+nrows*rh}" stroke="#e1e0d9"/>')
    s.append(f'<line x1="{Gx+4*cw}" y1="{y0-4}" x2="{Gx+4*cw}" y2="{y0+nrows*rh}" stroke="#e1e0d9"/>')

    def chip(x,y,fam,status):
        fill,tc=("#185FA5","#ffffff") if status=="site" else ("#B5D4F4","#0C447C")
        return (f'<rect x="{x}" y="{y}" width="22" height="13" rx="3" fill="{fill}"/>'
                f'<text x="{x+11}" y="{y+9.6:.0f}" text-anchor="middle" font-size="8" fill="{tc}">{ABBR[fam]}</text>')
    ri=0
    for gname,mets in GROUPS:
        gy0=y0+ri*rh; gyc=gy0+len(mets)*rh/2
        s.append(f'<text x="14" y="{gyc:.0f}" text-anchor="middle" font-size="9.5" fill="#7a7768" transform="rotate(-90 14 {gyc:.0f})">{gname}</text>')
        s.append(f'<line x1="24" y1="{gy0}" x2="{Gx+4*cw}" y2="{gy0}" stroke="#c3c2b7"/>')
        for m in mets:
            ry=y0+ri*rh
            s.append(f'<text x="{Lx-12}" y="{ry+rh/2+4:.0f}" text-anchor="end" font-size="11.5" font-weight="500" fill="#1a1a1a">{m}</text>')
            for gi,(gk,gl) in enumerate(GENOMES):
                present=[f for f in metal2fam[m] if STATUS[f][gk]!="abs"]
                cellx=Gx+gi*cw+6; n=len(present); nr=(n+3)//4 if n else 0
                starty=ry+(rh-(nr*16-3))/2 if nr else ry
                for k,f in enumerate(present):
                    s.append(chip(cellx+(k%4)*24, starty+(k//4)*16, f, STATUS[f][gk]))
            ri+=1
    s.append(f'<line x1="24" y1="{y0+nrows*rh}" x2="{Gx+4*cw}" y2="{y0+nrows*rh}" stroke="#c3c2b7"/>')
    ly=y0+nrows*rh+22
    s.append(chip(Lx-6,ly-10,"MerR","site")); s.append(f'<text x="{Lx+22}" y="{ly}" font-size="10.5" fill="#555">familia con sitio CHED (confirmada)</text>')
    s.append(chip(Lx+254,ly-10,"MerR","cand")); s.append(f'<text x="{Lx+282}" y="{ly}" font-size="10.5" fill="#555">familia candidata (sin sitio)</text>')
    s.append(f'<text x="{Lx-6}" y="{ly+20}" font-size="10" fill="#777">Cada chip = una familia presente en ese genoma capaz de sensar ese metal. Celda vacía = metal no cubierto.</text>')
    key="Ars=ArsR  Mer=MerR  Fur=Fur  Dtx=DtxR  Rrf=Rrf2  Cso=CsoR  Cop=CopY  Mar=MarR  Gnt=GntR  Lys=LysR  Tet=TetR"
    s.append(f'<text x="{Lx-6}" y="{ly+38}" font-size="9.5" fill="#777">{key}</text>')
    s.append('</svg>')
    return "".join(s)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="metal_by_genome.svg")
    args = ap.parse_args()
    open(args.out, "w").write(build_svg())
    print("escrito", args.out)
