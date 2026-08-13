"""
Anota la salida de MetalNet2 y produce tablas a nivel par, proteina, perfil,
familia, familia x metal, y genoma.

- 'perfil' = FPDB/HMM (BsCzrA...). 'familia' = grupo biologico (ArsR...).
  Mapeo perfil->familia y metal en config/profiles.tsv.
- REGULADOR = proteina unica (orig_id) por genoma. Perfiles de la misma familia
  se colapsan (cuenta 1 vez a nivel familia/genoma).
- 'sitio CHDE' = proteina con >=1 par que pasa el criterio. Por defecto exige
  filter_by_graph==1 (el veredicto de MetalNet, que reproduce el analisis original).
"""
import pandas as pd
from pathlib import Path

SEP = "__"

pairs_path    = snakemake.input.pairs      # noqa: F821
prov_path     = snakemake.input.prov       # noqa: F821
profiles_path = snakemake.input.profiles   # noqa: F821

out_pairs       = snakemake.output.pairs        # noqa: F821
out_perprot     = snakemake.output.perprot      # noqa: F821
out_perprof     = snakemake.output.perprof      # noqa: F821
out_perfam      = snakemake.output.perfam       # noqa: F821
out_perfammetal = snakemake.output.perfam_metal # noqa: F821
out_pergen      = snakemake.output.pergen       # noqa: F821
out_famids      = snakemake.output.famids       # noqa: F821

min_prob      = float(snakemake.params.min_prob)       # noqa: F821
require_graph = bool(snakemake.params.require_graph)   # noqa: F821

for p in [out_pairs, out_perprot, out_perprof, out_perfam, out_perfammetal, out_pergen, out_famids]:
    Path(p).parent.mkdir(parents=True, exist_ok=True)

pairs = pd.read_csv(pairs_path, sep="\t")
prov  = pd.read_csv(prov_path, sep="\t")
prof  = pd.read_csv(profiles_path, sep="\t")
if "metal" not in prof.columns:
    prof["metal"] = ""

prov = prov.rename(columns={"family": "profile"})

def add_family(df):
    df = df.merge(prof, on="profile", how="left")
    miss = df["family"].isna()
    if miss.any():
        print(f"[parse_metalnet] AVISO: perfiles sin familia: {sorted(df.loc[miss,'profile'].unique())}")
        df.loc[miss, "family"] = df.loc[miss, "profile"]
    df["metal"] = df["metal"].fillna("")
    return df

prov = add_family(prov)

# ---- PARES anotados ----
id_col = "seq_id" if "seq_id" in pairs.columns else pairs.columns[0]
pairs = pairs.rename(columns={id_col: "seq_id"})
parts = pairs["seq_id"].astype(str).str.split(SEP, n=2, expand=True)
if parts.shape[1] == 3:
    pairs["sample"], pairs["profile"], pairs["orig_id"] = parts[0], parts[1], parts[2]
pairs = add_family(pairs)
pairs.to_csv(out_pairs, sep="\t", index=False)

# ---- criterio de sitio ----
site_mask = pd.Series(True, index=pairs.index)
if min_prob > 0 and "prob" in pairs.columns:
    site_mask &= pairs["prob"].astype(float) >= min_prob
if require_graph and "filter_by_graph" in pairs.columns:
    site_mask &= pairs["filter_by_graph"].astype(int) == 1
pairs_site = pairs[site_mask].copy()

if len(pairs_site):
    site_prot = set(map(tuple, pairs_site[["sample", "orig_id"]].drop_duplicates().values.tolist()))
    npairs_prot = pairs_site.groupby(["sample", "orig_id"]).size().rename("n_pairs").reset_index()
else:
    site_prot = set()
    npairs_prot = pd.DataFrame(columns=["sample", "orig_id", "n_pairs"])

def has(s, o):
    return (s, o) in site_prot

# ---- por PROTEINA ----
fam_by_prot = (prov.groupby(["sample", "orig_id"])["family"]
               .agg(lambda s: ",".join(sorted(set(s)))).rename("families").reset_index())
prof_by_prot = (prov.groupby(["sample", "orig_id"])["profile"]
                .agg(lambda s: ",".join(sorted(set(s)))).rename("profiles").reset_index())
perprot = fam_by_prot.merge(prof_by_prot, on=["sample", "orig_id"]).merge(
    npairs_prot, on=["sample", "orig_id"], how="left")
perprot["n_pairs"] = perprot["n_pairs"].fillna(0).astype(int)
perprot["has_chde_site"] = perprot["n_pairs"] > 0
perprot.sort_values(["sample", "families", "orig_id"]).to_csv(out_perprot, sep="\t", index=False)

# ---- por PERFIL ----
reg_prof = prov.drop_duplicates(["sample", "profile", "family", "metal", "orig_id"]).copy()
reg_prof["site"] = [has(s, o) for s, o in zip(reg_prof["sample"], reg_prof["orig_id"])]
perprof = (reg_prof.groupby(["sample", "profile", "family", "metal"])
           .agg(n_regulators=("orig_id", "nunique"), n_with_chde_site=("site", "sum")).reset_index())
perprof["n_with_chde_site"] = perprof["n_with_chde_site"].astype(int)
perprof["frac_with_site"] = (perprof["n_with_chde_site"] / perprof["n_regulators"]).round(3)
perprof.sort_values(["sample", "family", "profile"]).to_csv(out_perprof, sep="\t", index=False)

# ---- por FAMILIA (hits unicos) ----
reg_fam = prov.drop_duplicates(["sample", "family", "orig_id"]).copy()
reg_fam["site"] = [has(s, o) for s, o in zip(reg_fam["sample"], reg_fam["orig_id"])]
perfam = (reg_fam.groupby(["sample", "family"])
          .agg(n_regulators=("orig_id", "nunique"), n_with_chde_site=("site", "sum")).reset_index())
perfam["n_with_chde_site"] = perfam["n_with_chde_site"].astype(int)
perfam["n_without_site"] = perfam["n_regulators"] - perfam["n_with_chde_site"]
perfam.sort_values(["sample", "family"]).to_csv(out_perfam, sep="\t", index=False)

# ---- FAMILIA x METAL (colapsando perfiles; metales = union de los perfiles) ----
fam_metals = {}
for _, r in prof.iterrows():
    met = str(r.get("metal", "") or "")
    if not met or met.lower() == "nan":
        continue
    for m in met.split("/"):
        m = m.strip()
        if m:
            fam_metals.setdefault(r["family"], set()).add(m)

rows = []
for _, r in perfam.iterrows():
    mets = sorted(fam_metals.get(r["family"], [])) or ["(sin metal)"]
    for m in mets:
        rows.append((r["sample"], m, r["family"], r["n_regulators"], r["n_with_chde_site"]))
perfam_metal = pd.DataFrame(rows, columns=["sample", "metal", "family", "n_regulators", "n_with_chde_site"])
perfam_metal["frac_with_site"] = (perfam_metal["n_with_chde_site"] / perfam_metal["n_regulators"]).round(3)
perfam_metal.to_csv(out_perfammetal, sep="\t", index=False)

# ---- LISTA DE IDs POR FAMILIA (numeracion por proteina, reinicia por genoma) ----
# Da el "ArsR1, ArsR2, ..." que permite rastrear cada proteina; la numeracion
# reinicia en cada genoma porque las proteinas son distintas.
famids = reg_fam.sort_values(["sample", "family", "orig_id"]).copy()
famids["family_index"] = famids.groupby(["sample", "family"]).cumcount() + 1
famids["label"] = famids["family"] + famids["family_index"].astype(str)
famids = famids.rename(columns={"site": "has_chde_site"})
# perfiles que detectaron cada proteina (por si una proteina la ven 2 perfiles)
prof_by = (prov.groupby(["sample", "family", "orig_id"])["profile"]
           .agg(lambda s: ",".join(sorted(set(s)))).rename("profiles").reset_index())
famids = famids.merge(prof_by, on=["sample", "family", "orig_id"], how="left")
famids = famids[["sample", "family", "family_index", "label", "orig_id", "profiles", "has_chde_site"]]
famids.to_csv(out_famids, sep="\t", index=False)

# ---- por GENOMA ----
pergen = reg_fam.groupby("sample").agg(n_regulators=("orig_id", "nunique")).reset_index()
ws = reg_fam[reg_fam["site"]].groupby("sample")["orig_id"].nunique().rename("n_with_chde_site").reset_index()
pergen = pergen.merge(ws, on="sample", how="left")
pergen["n_with_chde_site"] = pergen["n_with_chde_site"].fillna(0).astype(int)
pergen.to_csv(out_pergen, sep="\t", index=False)

print(f"[parse_metalnet] reguladores unicos: {int(pergen['n_regulators'].sum())} | "
      f"con sitio CHDE: {int(pergen['n_with_chde_site'].sum())} | "
      f"(min_prob={min_prob}, require_graph={require_graph})")
