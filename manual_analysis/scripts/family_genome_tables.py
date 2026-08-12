import glob
from collections import defaultdict
def parse(fp):
    h,s=None,[]
    for line in open(fp):
        line=line.rstrip()
        if not line: continue
        if line.startswith(">"):
            if h: yield h,"".join(s)
            h,s=line[1:].split()[0],[]
        else: s.append(line)
    if h: yield h,"".join(s)
seq2id={s:i for i,s in parse("metalnet_input/candidates.fasta")}
has={}
with open("full/site_summary.tsv") as fh:
    hdr=fh.readline().rstrip().split("\t"); ix={c:i for i,c in enumerate(hdr)}
    for line in fh:
        r=line.rstrip().split("\t"); has[r[ix["seq_id"]]]=(r[ix["has_site"]]=="yes")

# familia -> cognado metal (aproximado, para reporte)
metal={"ecFur":"Fe","ecZur":"Zn","ecZntR":"Zn","DtxR_bsMntR":"Mn/Fe","CsoR_bsCsoR":"Cu",
       "PbrR":"Pb","BsCzrA":"Zn/Co","mtNmtR":"Ni/Co","Rrf2_ecIscR":"Fe-S","MarR_AdcR":"Zn",
       "CopY_saMecI":"Cu","LysR_ModE":"Mo","TetR_SczA":"Zn","GntR_LldR":"(vario)"}

gf=defaultdict(set); gfh=defaultdict(set)
fam_all=defaultdict(set); fam_h=defaultdict(set)
org_all=defaultdict(set); org_h=defaultdict(set)
for f in sorted(glob.glob("bita/Results_Bita_-5/*/*/*.fasta")):
    org=f.split("/")[-3].replace("_e-5",""); fam=f.split("/")[-2]
    for h,s in parse(f):
        sid=seq2id.get(s.rstrip("*").upper())
        if not sid: continue
        hit=has.get(sid,False)
        gf[(org,fam)].add(sid); fam_all[fam].add((org,sid)); org_all[org].add(sid)
        if hit: gfh[(org,fam)].add(sid); fam_h[fam].add((org,sid)); org_h[org].add(sid)
orgs=sorted(org_all); fams=sorted(fam_all)

with open("full/report/family_success.tsv","w") as o:
    o.write("family\tcognate_metal\tcandidates\twith_site\tpct\n")
    for fam in sorted(fams,key=lambda x:-(len(fam_h[x])/len(fam_all[x]))):
        n=len(fam_all[fam]); hh=len(fam_h[fam])
        o.write(f"{fam}\t{metal.get(fam,'?')}\t{n}\t{hh}\t{round(100*hh/n)}\n")

with open("full/report/genome_success.tsv","w") as o:
    o.write("genome\tcandidates\twith_site\tpct\n")
    for org in orgs:
        n=len(org_all[org]); hh=len(org_h[org])
        o.write(f"{org}\t{n}\t{hh}\t{round(100*hh/n)}\n")

with open("full/report/genome_family_matrix.tsv","w") as o:
    o.write("family\t"+"\t".join(orgs)+"\n")
    for fam in fams:
        cells=[]
        for org in orgs:
            n=len(gf[(org,fam)]); hh=len(gfh[(org,fam)])
            cells.append(f"{hh}/{n}" if n else "NA")
        o.write(fam+"\t"+"\t".join(cells)+"\n")
print("tablas escritas en full/report/")
