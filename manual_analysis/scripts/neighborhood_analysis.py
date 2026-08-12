import re
from collections import defaultdict

GFF={"Mavi_homi":"nb/neihgborhood_data/gffs/Mavi_hom.gff",
     "MtubH37Rv":"nb/neihgborhood_data/gffs/H37Rv.gff",
     "VchoRFB16":"nb/neihgborhood_data/gffs/Vcho.gff",
     "VvulNBRC":"nb/neihgborhood_data/gffs/Vvul.gff"}
def attr(a,k):
    m=re.search(rf'(?:^|;){k}=([^;]+)',a); return m.group(1) if m else ""

CATS={
 'metal_efflux_transport':{
  'sym':['copa','copb','copc','copd','copz','znta','zntb','zupt','znua','znub','znuc','nika','nikb',
         'nikc','nikd','nike','feob','feoa','arsb','acr3','czca','czcb','czcc','czcd','cnra','cnrb','cnrc',
         'cora','mgta','mgtb','mgte','mnth','sita','sitb','sitc','sitd','troa','moda','modb','modc',
         'adca','adcb','adcc','dmef','yiip','zosa','ctpv','ctpc','ctpg','ctpj','pbra','pbrt','pbrb','nrsd'],
  'phr':['p-type atpase','cation diffusion facilitator','cation efflux','heavy metal translocating',
         'zinc transporter','zinc-transporting','zinc abc','nickel transporter','nickel abc','nickel/cobalt',
         'cobalt transporter','cobalt-zinc-cadmium','manganese transporter','manganese abc','copper-transporting',
         'copper-translocating','copper abc','copper resistance','copper efflux','copper-exporting','molybdate',
         'tungstate','cadmium','mercuric transport','heavy metal','metal abc transporter','metal transporter',
         'divalent metal','divalent cation','cation-transporting','metal cation','arsenical','arsenite efflux',
         'magnesium/cobalt','lead-translocating']},
 'metal_resistance_enzyme':{
  'sym':['arsc','arsi','arsm','arsh','mera','merb','merc','merp','cadi','cueo','czcm'],
  'phr':['arsenate reductase','arsenite methyltransferase','mercuric reductase','alkylmercury',
         'tellurite','copper oxidase','multicopper oxidase','cupric reductase']},
 'iron_siderophore':{
  'sym':['feca','fecb','fecc','fecd','fece','fhua','fhub','fhuc','fhud','feua','feub','iuca','iucb','iucc',
         'iucd','irga','tonb','exbb','exbd','vcta','vctc','vctd','vctg','vctp','viua','viub','hmut','hmuu',
         'hmuv','entc','ente','entf','fepa','fepb','fepc','fepd','fepg','irtb'],
  'phr':['siderophore','ferric','ferrous','iron uptake','iron chelate','iron abc','iron transport',
         'iron-regulated','iron-siderophore','enterobactin','vibrioferrin','vibriobactin','mycobactin',
         'hemin','heme transport','heme abc','heme uptake','tonb-dependent','ferrichrome','ferric uptake',
         'iron complex','ferrioxamine']},
 'metal_storage_chelation':{
  'sym':['dps','ftna','ftnb','bfr','bfd','mymt'],
  'phr':['ferritin','bacterioferritin','metallothionein','metallochaperone','copper chaperone','iron storage']},
}
METAL_CATS=set(CATS)

def categorize(product,gene):
    txt=((gene or '')+' '+(product or '')).lower()
    for cat,kw in CATS.items():
        for s in kw['sym']:
            if re.search(rf'\b{s}\b',txt): return cat
        for p in kw['phr']:
            if p in txt: return cat
    if any(k in txt for k in ['transcriptional regulator','repressor','activator','dna-binding','helix-turn-helix']):
        return 'other_regulator'
    if any(k in txt for k in ['hypothetical','duf','unknown','uncharacterized']):
        return 'hypothetical_uncharacterized'
    if any(k in txt for k in ['transporter','permease','efflux','abc transporter']):
        return 'generic_transport_nonmetal'
    return 'other'

def load_gff(fp):
    bc=defaultdict(list)
    for line in open(fp):
        if line.startswith("#"): continue
        c=line.rstrip("\n").split("\t")
        if len(c)<9 or c[2]!="CDS": continue
        pid=attr(c[8],"protein_id") or attr(c[8],"Name")
        if not pid: continue
        bc[c[0]].append([int(c[3]),int(c[4]),c[6],pid,attr(c[8],"product"),attr(c[8],"gene")])
    for k in bc: bc[k].sort()
    return bc
GEN={g:load_gff(fp) for g,fp in GFF.items()}
def gap(a1,a2,b1,b2):
    if a2>=b1 and b2>=a1: return 0
    return (b1-a2-1) if b1>a2 else (a1-b2-1)
def neigh(genome,acc):
    bc=GEN[genome]
    for contig,genes in bc.items():
        for i,g in enumerate(genes):
            if g[3]==acc:
                R=g; out=[]
                for j in range(max(0,i-6),min(len(genes),i+7)):
                    if j==i: continue
                    G=genes[j]; d=gap(R[0],R[1],G[0],G[1]); cat=categorize(G[4],G[5])
                    out.append({'gap':d,'gene':G[5],'prod':G[4],'cat':cat,
                                'side':'up' if G[1]<R[0] else 'down','metal':cat in METAL_CATS})
                out.sort(key=lambda x:x['gap'])
                return {'contig':contig,'strand':R[2],'prod':R[4],'neigh':out}
    return None

# re-validar el falso positivo y los efectivos
print("Chequeo polymerase:", categorize("RNA polymerase sigma factor SigB","sigB"), "(debe ser other_regulator/other)")
for g,a,lbl in [("MtubH37Rv","NP_217227.1","IdeR/DtxR"),("VvulNBRC","WP_017419581.1","MarR/dmeF"),
                ("VvulNBRC","WP_017422578.1","ZntR/vctC")]:
    nb=neigh(g,a); mets=[f"{(n['gene'] or n['prod'][:30])}({n['gap']}bp,{n['cat']})" for n in nb['neigh'] if n['metal'] and n['gap']<=500]
    print(f"  {lbl} {a}: metálicos≤500 = {mets if mets else 'ninguno'}")

def parse_master():
    rows=[]
    with open("full/report/master_candidates.tsv") as fh:
        hdr=fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r=line.rstrip("\n").split("\t"); rows.append(dict(zip(hdr,r)))
    return rows

rows=parse_master(); allrec=[]
for r in rows:
    nb=neigh(r["genome"],r["accession"]); rec=dict(r)
    if nb is None:
        rec.update(dict(contig="NA",up1="",down1="",n500=0,m200="NA",m500="NA",mname="",nm="",nmg=""))
    else:
        ne=nb["neigh"]; w500=[n for n in ne if n["gap"]<=500]
        m200=[n for n in w500 if n["metal"] and n["gap"]<=200]; m500=[n for n in w500 if n["metal"]]
        mets=[n for n in ne if n["metal"]]
        up1=next((n for n in ne if n["side"]=="up"),None); dn1=next((n for n in ne if n["side"]=="down"),None)
        f=lambda n: f"{(n['gene'] or (n['prod'] or '')[:30])}|{n['gap']}bp|{n['cat']}" if n else ""
        rec.update(dict(contig=nb["contig"],up1=f(up1),down1=f(dn1),n500=len(w500),
            m200="yes" if m200 else "no",m500="yes" if m500 else "no",
            mname=";".join(sorted(set((n['gene'] or (n['prod'] or '')[:30]) for n in m500))),
            nm=(mets[0]['gene'] or (mets[0]['prod'] or '')[:30]) if mets else "",
            nmg=mets[0]['gap'] if mets else ""))
    allrec.append(rec)

def eff(site,metal): return "NA" if site=="NA" else ("yes" if site=="yes" and metal=="yes" else "no")

# neighborhood_table
nc=["genome","accession","families","has_metalnet_site","contig","up1","down1","n500","m200","m500","mname","nm","nmg"]
with open("full/report/neighborhood_table.tsv","w") as o:
    o.write("genome\taccession\tfamilies\thas_metalnet_site\tcontig\tupstream_1\tdownstream_1\tn_neighbors_500bp\tmetal_within_200bp\tmetal_within_500bp\tmetal_neighbor\tnearest_metal_gene\tnearest_metal_gap_bp\n")
    for r in allrec: o.write("\t".join(str(r.get(c,"")) for c in nc)+"\n")

# master completo
mc=["genome","accession","families","n_families","has_metalnet_site","site_residues","composition","max_prob",
    "n_neighbor_genes_500bp","metal_gene_within_200bp","metal_gene_within_500bp","metal_neighbor","effective_200bp","effective_500bp"]
with open("full/report/master_candidates_completed.tsv","w") as o:
    o.write("\t".join(mc)+"\n")
    for r in allrec:
        r["n_neighbor_genes_500bp"]=r["n500"]; r["metal_gene_within_200bp"]=r["m200"]
        r["metal_gene_within_500bp"]=r["m500"]; r["metal_neighbor"]=r["mname"]
        r["effective_200bp"]=eff(r["has_metalnet_site"],r["m200"]); r["effective_500bp"]=eff(r["has_metalnet_site"],r["m500"])
        o.write("\t".join(str(r.get(c,"")) for c in mc)+"\n")

# RESUMEN + 2x2
site=lambda r: r["has_metalnet_site"]=="yes"
mn=lambda r: r["m500"]=="yes"
n=len(allrec)
s_only=sum(1 for r in allrec if site(r) and not mn(r))
m_only=sum(1 for r in allrec if mn(r) and not site(r))
both=sum(1 for r in allrec if site(r) and mn(r))
neither=sum(1 for r in allrec if not site(r) and not mn(r))
print("\n=== 2x2  (sitio MetalNet  x  gen metálico ≤500bp) ===")
print(f"                        metal+    metal-")
print(f"  sitio+ (MetalNet)       {both:>3}      {s_only:>3}")
print(f"  sitio-                  {m_only:>3}      {neither:>3}")
print(f"\n  EFECTIVOS (ambos, ≤500bp): {both}   |   (≤200bp): {sum(1 for r in allrec if site(r) and r['m200']=='yes')}")
print("\n=== Lista de EFECTIVOS (sitio + vecino metálico ≤500bp) ===")
for r in allrec:
    if site(r) and mn(r):
        print(f"  {r['genome']:<11} {r['accession']:<16} [{r['families']}]")
        print(f"        sitio={r['site_residues']}  ({r['composition']}, p={r['max_prob']})")
        print(f"        vecino metálico={r['mname']}  | flanco: {r['up1']}  //  {r['down1']}")
print("\nArchivos: neighborhood_table.tsv, master_candidates_completed.tsv")
