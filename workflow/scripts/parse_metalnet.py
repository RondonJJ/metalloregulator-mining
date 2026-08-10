"""
Anota la salida de MetalNet2 (pred_pairs.tsv): recupera muestra/familia desde el
header (muestra|familia|id_original) y la cruza con provenance.tsv.
"""
import pandas as pd
from pathlib import Path

pairs_path = snakemake.input.pairs   # noqa: F821
prov_path  = snakemake.input.prov    # noqa: F821
out_path   = snakemake.output.table  # noqa: F821

pairs = pd.read_csv(pairs_path, sep="\t")
prov  = pd.read_csv(prov_path, sep="\t")

# >>> EDITAR (2/3): nombre real de la columna de ID de proteina en pred_pairs.tsv <<<
# Se intenta detectar automaticamente; si falla, fija id_col a mano.
id_col = None
for c in ["seq_id", "protein", "id", "name", "query", "sequence"]:
    if c in pairs.columns:
        id_col = c
        break
if id_col is None:
    id_col = pairs.columns[0]
    print(f"[parse_metalnet] AVISO: uso la primera columna '{id_col}' como ID.")

pairs = pairs.rename(columns={id_col: "seq_id"})

# recuperar muestra|familia|id_original desde el header
parts = pairs["seq_id"].astype(str).str.split("|", n=2, expand=True)
if parts.shape[1] == 3:
    pairs["sample"]  = parts[0]
    pairs["family"]  = parts[1]
    pairs["orig_id"] = parts[2]

merged = pairs.merge(prov, on="seq_id", how="left", suffixes=("", "_prov"))

Path(out_path).parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(out_path, sep="\t", index=False)
print(f"[parse_metalnet] {len(merged)} filas anotadas -> {out_path}")
