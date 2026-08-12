import re
# ---- datos ----
fammetals={"ArsR":["Co","Ni","Cu","Zn","Cd","Hg","Pb","Bi","As","Sb"],
"MerR":["Fe","Co","Ni","Cu","Zn","Ag","Cd","Au","Hg","Pb"],"Fur":["Mn","Fe","Zn","Ni","Co"],
"DtxR":["Mn","Fe"],"Rrf2":["Fe"],"CsoR":["Co","Ni","Cu","Zn"],"CopY":["Cu","Zn"],
"MarR":["Cu","Zn"],"GntR":["Cu","Zn"],"LysR":["Mo","W"],"TetR":["Fe","Cu","Zn","As"]}
abbr={"ArsR":"Ars","MerR":"Mer","Fur":"Fur","DtxR":"Dtx","Rrf2":"Rrf","CsoR":"Cso",
"CopY":"Cop","MarR":"Mar","GntR":"Gnt","LysR":"Lys","TetR":"Tet"}
famorder=["ArsR","MerR","Fur","DtxR","Rrf2","CsoR","CopY","MarR","GntR","LysR","TetR"]
st={"ArsR":{"Ma":"site","Mt":"site","Vc":"cand","Vv":"cand"},
"MerR":{"Ma":"site","Mt":"cand","Vc":"site","Vv":"site"},
"Fur":{"Ma":"site","Mt":"site","Vc":"site","Vv":"site"},
"DtxR":{"Ma":"site","Mt":"site","Vc":"abs","Vv":"abs"},
"Rrf2":{"Ma":"cand","Mt":"cand","Vc":"site","Vv":"site"},
"CsoR":{"Ma":"site","Mt":"site","Vc":"abs","Vv":"abs"},
"CopY":{"Ma":"cand","Mt":"cand","Vc":"abs","Vv":"abs"},
"MarR":{"Ma":"cand","Mt":"cand","Vc":"cand","Vv":"site"},
"GntR":{"Ma":"site","Mt":"site","Vc":"site","Vv":"site"},
"LysR":{"Ma":"cand","Mt":"abs","Vc":"abs","Vv":"abs"},
"TetR":{"Ma":"cand","Mt":"cand","Vc":"abs","Vv":"cand"}}
metal2fam={}
for f in famorder:
    for m in fammetals[f]: metal2fam.setdefault(m,[]).append(f)
genomes=[("Ma","M. av."),("Mt","M. tb."),("Vc","V. ch."),("Vv","V. vu.")]

# ---- layout panel E (mm), YOFF baja todo bajo la Fig5 ----
YOFF=184.0
blocks=[("Fe","Mn","Zn","Ni","Co","Cu","Cd","Pb"),("Hg","Ag","Au","As","Sb","Bi","Mo","W")]
group_of={m:g for g,ms in [("esenciales",["Fe","Mn","Zn","Ni","Co","Cu"]),
 ("tóxicos",["Cd","Pb","Hg","Ag","Au"]),("metaloides",["As","Sb","Bi"]),("oxianiones",["Mo","W"])] for m in ms}
bx=[7.0,110.0]; labw=10.0; cw=21.0; rh=7.0; y_hdr=YOFF+11; y0=YOFF+15
CH_W,CH_H=4.9,2.9
S=[]
def T(x,y,txt,fs,anchor="start",fill="#1a1a1a",bold=False,italic=False):
    st_=f' font-weight="bold"' if bold else ''
    it=' font-style="italic"' if italic else ''
    S.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{fs}" text-anchor="{anchor}" fill="{fill}"{st_}{it} font-family="Arial">{txt}</text>')
def chip(x,y,fam,status):
    fill,tc=("#185FA5","#ffffff") if status=="site" else ("#B5D4F4","#0C447C")
    S.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{CH_W}" height="{CH_H}" rx="0.7" fill="{fill}"/>')
    S.append(f'<text x="{x+CH_W/2:.2f}" y="{y+CH_H-0.85:.2f}" font-size="1.9" text-anchor="middle" fill="{tc}" font-family="Arial">{abbr[fam]}</text>')
# panel letter + titulo
T(2.5, YOFF+6.0, "E", 4.5, bold=True)
T(9, YOFF+5.8, "Cobertura potencial de metales por genoma", 3.1, bold=True)
# por bloque
for bi,mets in enumerate(blocks):
    X=bx[bi]
    # encabezados de genoma
    for gi,(gk,gl) in enumerate(genomes):
        cx=X+labw+gi*cw+cw/2
        T(cx, y_hdr, gl, 2.5, anchor="middle", italic=True)
    # separador superior
    S.append(f'<line x1="{X+0.5:.1f}" y1="{y0-2:.1f}" x2="{X+labw+4*cw:.1f}" y2="{y0-2:.1f}" stroke="#9a988c" stroke-width="0.3"/>')
    prevg=None
    for ri,m in enumerate(mets):
        ry=y0+ri*rh
        # etiqueta de grupo cuando cambia
        g=group_of[m]
        if g!=prevg:
            T(X+0.3, ry+2.2, g, 2.0, fill="#7a7768", italic=True)
            if ri>0: S.append(f'<line x1="{X+0.5:.1f}" y1="{ry-0.5:.1f}" x2="{X+labw+4*cw:.1f}" y2="{ry-0.5:.1f}" stroke="#cfcdc2" stroke-width="0.25"/>')
            prevg=g
        T(X+labw-1, ry+rh/2+1.0, m, 3.0, anchor="end", bold=True)
        for gi,(gk,gl) in enumerate(genomes):
            present=[f for f in metal2fam[m] if st[f][gk]!="abs"]
            n=len(present); nr=(n+3)//4 if n else 0
            cellx=X+labw+gi*cw+0.7
            starty=ry+(rh-(nr*(CH_H+0.4)-0.4))/2 if nr else ry
            for k,f in enumerate(present):
                cx=cellx+(k%4)*(CH_W+0.15); cy=starty+(k//4)*(CH_H+0.4)
                chip(cx,cy,f,st[f][gk])
    S.append(f'<line x1="{X+0.5:.1f}" y1="{y0+8*rh:.1f}" x2="{X+labw+4*cw:.1f}" y2="{y0+8*rh:.1f}" stroke="#9a988c" stroke-width="0.3"/>')
# leyenda
ly=y0+8*rh+6
chip(7,ly-2.4,"MerR","site"); T(13,ly,"familia con sitio CHED (confirmada)",2.3)
chip(70,ly-2.4,"MerR","cand"); T(76,ly,"familia candidata (sin sitio)",2.3)
T(7,ly+4,"Genomas: Ma=M. avium · Mt=M. tuberculosis · Vc=V. cholerae · Vv=V. vulnificus",2.1,fill="#555")
T(7,ly+7.6,"Familias: Ars=ArsR Mer=MerR Fur=Fur Dtx=DtxR Rrf=Rrf2 Cso=CsoR Cop=CopY Mar=MarR Gnt=GntR Lys=LysR Tet=TetR",2.1,fill="#555")
T(7,ly+11.2,"Cada chip = familia presente capaz de sensar ese metal; los metales son los inductores potenciales de la familia (especificidad no resuelta).",2.0,fill="#555")
panelE="\n".join(S)
Htot=ly+14

inner=open('fig5_inner.svg').read()
out=(f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
     f'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
     f'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.0.dtd" '
     f'width="210mm" height="{Htot:.1f}mm" viewBox="0 0 210 {Htot:.1f}">\n'
     f'<rect x="0" y="0" width="210" height="{Htot:.1f}" fill="#ffffff"/>\n'
     f'<svg x="0" y="0" width="210" height="177" viewBox="0 0 210 177" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.0.dtd">\n{inner}\n</svg>\n'
     f'{panelE}\n</svg>\n')
open('Fig5_with_panelE.svg','w').write(out)
print("combinado escrito. Altura total:", round(Htot,1),"mm  (A4=297mm)")
print("panel E ocupa y:",YOFF,"a",round(Htot,1))
