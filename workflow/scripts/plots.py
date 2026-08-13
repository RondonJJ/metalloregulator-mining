"""
Figuras del pipeline. SIEMPRE en SVG y PDF (nunca PNG). Textos en INGLES.

- metal_family_matrix : figura estilo Fig5 panel D. Una sola grilla metal x genoma;
  cada celda lista las FAMILIAS (nombre completo) que sensan ese metal, en azul
  oscuro si tienen >=1 sitio CHDE en ese genoma, claro si no. Superindice = nº de
  reguladores de esa familia en ese genoma (permite rastrear cuanto se repite una
  familia entre metales). El vinculo familia->metales sale de config/family_metals.tsv.
- regulators_per_family  : reguladores por familia (con/sin sitio CHDE).
- regulators_per_genome  : total por genoma.

Editables: METAL_CATEGORIES y FAMILY_ORDER (abajo). El mapeo familia->metales se
edita en config/family_metals.tsv.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from math import ceil
from pathlib import Path

# ---- Font settings ----
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.weight"] = "heavy"

METAL_CATEGORIES = [
    ("Essential metals", ["Fe", "Mn", "Zn", "Ni", "Co", "Cu"]),
    ("Xenobiotics",      ["Cd", "Pb", "Hg", "Ag", "Au"]),
    ("Metalloids",       ["As", "Sb", "Bi"]),
    ("Oxyanions",        ["Mo", "W"]),
]
# orden de las familias dentro de cada celda (como en la Fig5)
FAMILY_ORDER = ["ArsR", "MerR", "Fur", "CsoR", "CopY", "MarR", "GntR",
                "TetR", "DtxR", "Rrf2", "LysR", "NikR"]
PER_LINE = 4                 # familias por linea dentro de una celda
C_DARK = "#1f4e79"; C_LIGHT = "#c6dbef"
T_DARK = "white";   T_LIGHT = "#12335a"

perfam       = pd.read_csv(snakemake.input.perfam, sep="\t")         # noqa: F821
pergen       = pd.read_csv(snakemake.input.pergen, sep="\t")         # noqa: F821
family_metals_df = pd.read_csv(snakemake.input.family_metals, sep="\t")  # noqa: F821

out = snakemake.output  # noqa: F821
for p in out:
    Path(p).parent.mkdir(parents=True, exist_ok=True)


def save(fig, svg_path, pdf_path):
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


# ========== Figura 1: matriz metal x genoma (estilo Fig5 panel D) ==========
family_metals = {r["family"]: set(str(r["metals"]).split(","))
                 for _, r in family_metals_df.iterrows()}
# stats por (sample, family) -> (n_regulators, n_with_site)
stats = {(r["sample"], r["family"]): (int(r["n_regulators"]), int(r["n_with_chde_site"]))
         for _, r in perfam.iterrows()}
samples = sorted(perfam["sample"].unique())

sensed = set().union(*family_metals.values()) if family_metals else set()
rows = [(cat, m) for cat, mets in METAL_CATEGORIES for m in mets if m in sensed]
n_rows, n_cols = len(rows), len(samples)

fig, ax = plt.subplots(figsize=(n_cols * 2.5 + 2.4, n_rows * 0.62 + 1.6))
ax.set_xlim(-0.05, n_cols); ax.set_ylim(0, n_rows + 0.9); ax.axis("off")

# encabezados de columna (genomas)
for j, s in enumerate(samples):
    ax.text(j + 0.5, n_rows + 0.4, s, ha="center", va="center",
            fontsize=10, fontstyle="italic", fontweight="bold")
# lineas verticales
for j in range(n_cols + 1):
    ax.plot([j, j], [0, n_rows], color="#dddddd", lw=0.6, zorder=0)

# filas
for i, (cat, metal) in enumerate(rows):
    y_top = n_rows - i
    ax.text(-0.12, y_top - 0.5, metal, ha="right", va="center", fontsize=9, fontweight="bold")
    for j, s in enumerate(samples):
        fams = [f for f in FAMILY_ORDER
                if metal in family_metals.get(f, set()) and stats.get((s, f), (0, 0))[0] > 0]
        if not fams:
            continue
        lines = ceil(len(fams) / PER_LINE)
        for k, f in enumerate(fams):
            line, col = k // PER_LINE, k % PER_LINE
            x = j + (col + 0.5) * (1.0 / PER_LINE)
            y = y_top - (line + 0.5) * (1.0 / max(lines, 1))
            n_reg, n_site = stats[(s, f)]
            dark = n_site > 0
            ax.text(x, y, rf"$\mathregular{{{f}}}^{{{n_reg}}}$", ha="center", va="center",
                    fontsize=6.2, color=T_DARK if dark else T_LIGHT,
                    bbox=dict(boxstyle="round,pad=0.25",
                              facecolor=C_DARK if dark else C_LIGHT, edgecolor="none"))

# separadores horizontales por categoria + etiquetas
ax.plot([0, n_cols], [n_rows, n_rows], color="#333333", lw=1.0)
prev = None
cat_rows = {}
for i, (cat, metal) in enumerate(rows):
    cat_rows.setdefault(cat, []).append(i)
    if cat != prev:
        ax.plot([0, n_cols], [n_rows - i, n_rows - i], color="#333333", lw=1.0)
        prev = cat
ax.plot([0, n_cols], [0, 0], color="#333333", lw=1.0)
for cat, idxs in cat_rows.items():
    i0, i1 = min(idxs), max(idxs)
    y = (n_rows - i0) - (i1 - i0 + 1) / 2
    ax.text(-0.66, y, cat, ha="center", va="center", rotation=90,
            fontsize=9, fontstyle="italic", fontweight="bold")

legend = [Patch(facecolor=C_DARK, label="Family with CHDE site"),
          Patch(facecolor=C_LIGHT, label="Family without CHDE site")]
fig.legend(handles=legend, loc="lower center", ncol=2, frameon=False,
           fontsize=9, bbox_to_anchor=(0.5, -0.01))
ax.set_title("Metalloregulator families per metal and genome\n(superscript = number of regulators)",
             fontsize=12, pad=16)
save(fig, out.metal_matrix_svg, out.metal_matrix_pdf)


# ========== Figura 2: reguladores por familia ==========
C_SITE = "#1f4e79"; C_NONE = "#a6c8e8"
samples_f = sorted(perfam["sample"].unique())
n = len(samples_f)
fig, axes = plt.subplots(n, 1, figsize=(9, 3.0 * n), squeeze=False)
for ax, s in zip(axes[:, 0], samples_f):
    d = perfam[perfam["sample"] == s].sort_values("family")
    fams = d["family"].tolist()
    ax.bar(fams, d["n_with_chde_site"], color=C_SITE, label="with CHDE site")
    ax.bar(fams, d["n_without_site"], bottom=d["n_with_chde_site"], color=C_NONE, label="without site")
    ax.set_title(s); ax.set_ylabel("Regulators (unique hits)")
    ax.tick_params(axis="x", rotation=45)
    for x, t in zip(fams, d["n_regulators"]):
        ax.text(x, t, str(int(t)), ha="center", va="bottom", fontsize=8)
axes[0, 0].legend(frameon=False)
fig.suptitle("Metalloregulators per family")
fig.tight_layout()
save(fig, out.perfamily_svg, out.perfamily_pdf)


# ========== Figura 3: total por genoma ==========
fig, ax = plt.subplots(figsize=(1.6 * max(len(pergen), 2) + 2, 4.5))
x = range(len(pergen)); w = 0.4
ax.bar([i - w/2 for i in x], pergen["n_regulators"], width=w, color=C_NONE, label="total")
ax.bar([i + w/2 for i in x], pergen["n_with_chde_site"], width=w, color=C_SITE, label="with CHDE site")
ax.set_xticks(list(x)); ax.set_xticklabels(pergen["sample"])
ax.set_ylabel("Regulators (unique hits)"); ax.set_title("Total per genome")
for i, (a, b) in enumerate(zip(pergen["n_regulators"], pergen["n_with_chde_site"])):
    ax.text(i - w/2, a, str(int(a)), ha="center", va="bottom", fontsize=8)
    ax.text(i + w/2, b, str(int(b)), ha="center", va="bottom", fontsize=8)
ax.legend(frameon=False)
fig.tight_layout()
save(fig, out.pergenome_svg, out.pergenome_pdf)

print("[plots] metal_family_matrix + per_family + per_genome (SVG + PDF, ingles)")
