"""
Anota la salida de MetalNet2 y produce tablas a nivel par, proteina (regulador
unico), familia y genoma.

Claves del conteo (segun el analisis original):
- 'perfil' = FPDB/HMM (BsCzrA, mtNmtR, ecFur...). 'familia' = grupo biologico
  (ArsR, Fur, MerR...). El mapeo esta en config/profiles.tsv.
- Un REGULADOR es una proteina unica (orig_id) dentro de un genoma. Si la misma
  proteina la detectan 2 perfiles de la misma familia, cuenta UNA vez.
- 'sitio CHDE' = la proteina tiene >=1 par predicho por MetalNet. Por defecto SIN
  filtro (como en el analisis original); se puede exigir prob/filtro por config.
"""
import pandas as pd
from pathlib import Path

SEP = "__"

pairs_path    = snakemake.input.pairs      # noqa: F821
prov_path     = snakemake.input.prov       # noqa: F821
profiles_path = snakemake.input.profiles   # noqa: F821

out_pairs   = snakemake.output.pairs       # noqa: F821
out_perprot = snakemake.output.perprot     # noqa: F821
out_perfam  = snakemake.output.perfam      # noqa: F821
out_pergen  = snakemake.output.pergen      # noqa: F821

min_prob      = float(snakemake.params.min_prob)       # noqa: F821
require_graph = bool(snakemake.params.require_graph)   # noqa: F821

for p in [out_pairs, out_perprot, out_perfam, out_pergen]:
    Path(p).parent.mkdir(parents=True, exist_ok=True)

# ---- cargar ----
pairs = pd.read_csv(pairs_path, sep="\t")
prov  = pd.read_csv(prov_path, sep="\t")          # seq_id, sample, family(=perfil), orig_id
prof  = pd.read_csv(profiles_path, sep="\t")      # profile, family

# en provenance, la columna 'family' es en realidad el PERFIL -> renombrar
prov = prov.rename(columns={"family": "profile"})

def add_family(df):
    df = df.merge(prof, on="profile", how="left")
    miss = df["family"].isna()
    if miss.any():
        faltan = sorted(df.loc[miss, "profile"].unique())
        print(f"[parse_metalnet] AVISO: perfiles sin familia en profiles.tsv: {faltan}")
        df.loc[miss, "family"] = df.loc[miss, "profile"]
    return df

prov = add_family(prov)

# ---- tabla de PARES anotada ----
id_col = "seq_id" if "seq_id" in pairs.columns else pairs.columns[0]
pairs = pairs.rename(columns={id_col: "seq_id"})
parts = pairs["seq_id"].astype(str).str.split(SEP, n=2, expand=True)
if parts.shape[1] == 3:
    pairs["sample"], pairs["profile"], pairs["orig_id"] = parts[0], parts[1], parts[2]
pairs = add_family(pairs)
pairs.to_csv(out_pairs, sep="\t", index=False)

# ---- definir "sitio" segun filtros (por defecto: ninguno) ----
site_mask = pd.Series(True, index=pairs.index)
if min_prob > 0 and "prob" in pairs.columns:
    site_mask &= pairs["prob"].astype(float) >= min_prob
if require_graph and "filter_by_graph" in pairs.columns:
    site_mask &= pairs["filter_by_graph"].astype(int) == 1
pairs_site = pairs[site_mask].copy()

# proteinas con sitio: (sample, orig_id) con >=1 par que pasa el filtro
if len(pairs_site):
    site_prot = set(map(tuple, pairs_site[["sample", "orig_id"]].drop_duplicates().values.tolist()))
    npairs_prot = (pairs_site.groupby(["sample", "orig_id"]).size()
                   .rename("n_pairs").reset_index())
else:
    site_prot = set()
    npairs_prot = pd.DataFrame(columns=["sample", "orig_id", "n_pairs"])

# ---- resumen por PROTEINA (regulador unico) ----
fam_by_prot = (prov.groupby(["sample", "orig_id"])["family"]
               .agg(lambda s: ",".join(sorted(set(s)))).rename("families").reset_index())
prof_by_prot = (prov.groupby(["sample", "orig_id"])["profile"]
                .agg(lambda s: ",".join(sorted(set(s)))).rename("profiles").reset_index())
perprot = fam_by_prot.merge(prof_by_prot, on=["sample", "orig_id"], how="left")
perprot = perprot.merge(npairs_prot, on=["sample", "orig_id"], how="left")
perprot["n_pairs"] = perprot["n_pairs"].fillna(0).astype(int)
perprot["has_chde_site"] = perprot["n_pairs"] > 0
perprot = perprot.sort_values(["sample", "families", "orig_id"])
perprot.to_csv(out_perprot, sep="\t", index=False)

# ---- resumen por FAMILIA (hits unicos) ----
reg_fam = prov.drop_duplicates(["sample", "family", "orig_id"]).copy()
reg_fam["has_site"] = [ (s, o) in site_prot for s, o in zip(reg_fam["sample"], reg_fam["orig_id"]) ]
perfam = (reg_fam.groupby(["sample", "family"])
          .agg(n_regulators=("orig_id", "nunique"),
               n_with_chde_site=("has_site", "sum"))
          .reset_index())
perfam["n_with_chde_site"] = perfam["n_with_chde_site"].astype(int)
perfam["n_without_site"] = perfam["n_regulators"] - perfam["n_with_chde_site"]
perfam = perfam.sort_values(["sample", "family"])
perfam.to_csv(out_perfam, sep="\t", index=False)

# ---- resumen por GENOMA (totales unicos) ----
pergen = (reg_fam.groupby("sample")
          .agg(n_regulators=("orig_id", "nunique")).reset_index())
with_site = (reg_fam[reg_fam["has_site"]].groupby("sample")["orig_id"].nunique()
             .rename("n_with_chde_site").reset_index())
pergen = pergen.merge(with_site, on="sample", how="left")
pergen["n_with_chde_site"] = pergen["n_with_chde_site"].fillna(0).astype(int)
pergen.to_csv(out_pergen, sep="\t", index=False)

tot = int(pergen["n_regulators"].sum())
tot_site = int(pergen["n_with_chde_site"].sum())
print(f"[parse_metalnet] reguladores unicos: {tot} | con sitio CHDE: {tot_site}")
print(f"[parse_metalnet] filtro sitio -> min_prob={min_prob}, require_graph={require_graph}")
