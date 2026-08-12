"""
Figuras del pipeline. SIEMPRE en SVG y PDF (nunca PNG). Textos en INGLES.

- regulators_per_family_metal : heatmap familia x metal, un panel por genoma,
  metales agrupados por categoria, escala por tramos (estilo Fig5). Sugiere con
  que familias cuenta cada genoma para sensar cada metal. Perfiles colapsados a
  familia; metal = union de los perfiles.
- regulators_per_family : reguladores por familia (con/sin sitio CHDE).
- regulators_per_genome : total por genoma.

Para reordenar/recategorizar metales, edita METAL_CATEGORIES abajo.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import matplotlib.transforms as mtransforms
from pathlib import Path

# --- categorias de metales (editable). Solo se muestran los presentes en los datos ---
METAL_CATEGORIES = [
    ("Essential metals", ["Fe", "Mn", "Zn", "Ni", "Co", "Cu"]),
    ("Xenobiotics",      ["Cd", "Pb", "Hg", "Ag", "Au"]),
    ("Metalloids",       ["As", "Sb", "Bi"]),
    ("Oxyanions",        ["Mo", "W"]),
    ("Other",            ["Fe-S", "vario"]),
]
# etiquetas de metal en ingles para el eje
METAL_LABELS = {"vario": "various", "(sin metal)": "unassigned"}

# escala por tramos (estilo Fig5): 0 / 1-33 / 34-66 / 67-99 / 100 %
BIN_COLORS = ["#eef4fb", "#c6dbef", "#6baed6", "#2171b5", "#08306b"]
BIN_BOUNDS = [0.0, 0.0001, 0.335, 0.665, 0.9999, 1.0001]
BIN_LABELS = ["0%", "1–33%", "34–66%", "67–99%", "100%"]
NA_COLOR = "#d9d9d9"

C_SITE = "#1f4e79"
C_NONE = "#a6c8e8"

perfam       = pd.read_csv(snakemake.input.perfam, sep="\t")        # noqa: F821
pergen       = pd.read_csv(snakemake.input.pergen, sep="\t")        # noqa: F821
perfam_metal = pd.read_csv(snakemake.input.perfam_metal, sep="\t")  # noqa: F821

out = snakemake.output  # noqa: F821
for p in out:
    Path(p).parent.mkdir(parents=True, exist_ok=True)


def save(fig, svg_path, pdf_path):
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


# ========== Figura 1: familia x metal (un panel por genoma) ==========
present_metals = set(perfam_metal["metal"])
metal_order, cat_of = [], {}
for cat, mets in METAL_CATEGORIES:
    for m in mets:
        if m in present_metals and m not in metal_order:
            metal_order.append(m); cat_of[m] = cat
for m in sorted(present_metals):          # cualquiera no listado -> Other
    if m not in metal_order:
        metal_order.append(m); cat_of[m] = "Other"

families = sorted(perfam_metal["family"].unique())
samples = sorted(perfam_metal["sample"].unique())

cmap = ListedColormap(BIN_COLORS); cmap.set_bad(NA_COLOR)
norm = BoundaryNorm(BIN_BOUNDS, cmap.N)

ncol = 2
nrow = int(np.ceil(len(samples) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(1.0 * len(families) * ncol + 3,
                                              0.45 * len(metal_order) * nrow + 3),
                         squeeze=False)

for idx, s in enumerate(samples):
    ax = axes[idx // ncol][idx % ncol]
    sub = perfam_metal[perfam_metal["sample"] == s]
    lut = {(r["metal"], r["family"]): (r["frac_with_site"], int(r["n_with_chde_site"]), int(r["n_regulators"]))
           for _, r in sub.iterrows()}
    frac = np.full((len(metal_order), len(families)), np.nan)
    annot = np.empty((len(metal_order), len(families)), dtype=object); annot[:] = ""
    for i, m in enumerate(metal_order):
        for j, fam in enumerate(families):
            if (m, fam) in lut:
                f, nw, nt = lut[(m, fam)]
                frac[i, j] = f
                annot[i, j] = f"{nw}/{nt}"
    masked = np.ma.masked_invalid(frac)
    ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(families))); ax.set_xticklabels(families, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(metal_order)))
    ax.set_yticklabels([METAL_LABELS.get(m, m) for m in metal_order], fontsize=8)
    for i in range(len(metal_order)):
        for j in range(len(families)):
            if annot[i, j]:
                dark = (not np.isnan(frac[i, j])) and frac[i, j] >= 0.66
                ax.text(j, i, annot[i, j], ha="center", va="center", fontsize=6,
                        color="white" if dark else "#222222")
    # separadores entre categorias
    for i, m in enumerate(metal_order):
        if i != 0 and cat_of[m] != cat_of[metal_order[i - 1]]:
            ax.axhline(i - 0.5, color="black", lw=0.8)
    # etiquetas de categoria (solo columna izquierda)
    if idx % ncol == 0:
        trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
        spans, start = [], 0
        for k in range(1, len(metal_order) + 1):
            if k == len(metal_order) or cat_of[metal_order[k]] != cat_of[metal_order[start]]:
                spans.append((cat_of[metal_order[start]], start, k - 1)); start = k
        for cat, a, b in spans:
            ax.text(-0.32, (a + b) / 2, cat, transform=trans, rotation=90,
                    va="center", ha="center", fontsize=8, fontweight="bold", clip_on=False)
    ax.set_title(s, fontsize=10)

# apagar ejes sobrantes
for k in range(len(samples), nrow * ncol):
    axes[k // ncol][k % ncol].axis("off")

legend_patches = [Patch(facecolor=c, edgecolor="#999999", label=l)
                  for c, l in zip(BIN_COLORS, BIN_LABELS)]
legend_patches.append(Patch(facecolor=NA_COLOR, edgecolor="#999999", label="n/a"))
fig.legend(handles=legend_patches, title="% with CHDE site", loc="lower center",
           ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Metalloregulators per family and metal", fontsize=13)
fig.tight_layout(rect=[0, 0.03, 1, 0.97])
save(fig, out.fam_metal_svg, out.fam_metal_pdf)


# ========== Figura 2: reguladores por familia ==========
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

print("[plots] SVG + PDF generados (ingles, metales por categoria, escala por tramos)")
