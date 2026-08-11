"""
Figuras a partir de las tablas por familia y por genoma.
- regulators_per_family.png : reguladores por familia, apilando con/sin sitio CHDE
  (un subplot por genoma). Equivale al panel A + la distincion de color del panel D.
- regulators_per_genome.png : total de reguladores por genoma y cuantos con sitio
  (equivale al panel B).
Punto de partida: adapta colores/estilo a tu figura final.
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

perfam = pd.read_csv(snakemake.input.perfam, sep="\t")   # noqa: F821
pergen = pd.read_csv(snakemake.input.pergen, sep="\t")   # noqa: F821

out_family = snakemake.output.perfamily_png              # noqa: F821
out_genome = snakemake.output.pergenome_png              # noqa: F821
for p in [out_family, out_genome]:
    Path(p).parent.mkdir(parents=True, exist_ok=True)

C_SITE = "#1f4e79"     # con sitio CHDE (oscuro)
C_NONE = "#a6c8e8"     # sin sitio (claro)

# ---------- Figura 1: reguladores por familia (por genoma) ----------
samples = sorted(perfam["sample"].unique())
n = len(samples)
fig, axes = plt.subplots(n, 1, figsize=(9, 3.2 * n), squeeze=False)
for ax, s in zip(axes[:, 0], samples):
    sub = perfam[perfam["sample"] == s].sort_values("family")
    fams = sub["family"].tolist()
    ax.bar(fams, sub["n_with_chde_site"], color=C_SITE, label="con sitio CHDE")
    ax.bar(fams, sub["n_without_site"], bottom=sub["n_with_chde_site"],
           color=C_NONE, label="sin sitio")
    ax.set_title(s)
    ax.set_ylabel("Reguladores (hits unicos)")
    ax.tick_params(axis="x", rotation=45)
    # etiqueta del total encima de cada barra
    for x, tot in zip(fams, sub["n_regulators"]):
        ax.text(x, tot, str(int(tot)), ha="center", va="bottom", fontsize=8)
axes[0, 0].legend(frameon=False)
fig.suptitle("Reguladores por familia")
fig.tight_layout()
fig.savefig(out_family, dpi=200)
plt.close(fig)

# ---------- Figura 2: total por genoma ----------
fig, ax = plt.subplots(figsize=(1.6 * max(len(pergen), 2) + 2, 4.5))
x = range(len(pergen))
w = 0.4
ax.bar([i - w/2 for i in x], pergen["n_regulators"], width=w,
       color=C_NONE, label="total reguladores")
ax.bar([i + w/2 for i in x], pergen["n_with_chde_site"], width=w,
       color=C_SITE, label="con sitio CHDE")
ax.set_xticks(list(x))
ax.set_xticklabels(pergen["sample"], rotation=0)
ax.set_ylabel("Reguladores (hits unicos)")
ax.set_title("Total por genoma")
for i, (a, b) in enumerate(zip(pergen["n_regulators"], pergen["n_with_chde_site"])):
    ax.text(i - w/2, a, str(int(a)), ha="center", va="bottom", fontsize=8)
    ax.text(i + w/2, b, str(int(b)), ha="center", va="bottom", fontsize=8)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(out_genome, dpi=200)
plt.close(fig)

print(f"[plots] figuras -> {out_family} , {out_genome}")
