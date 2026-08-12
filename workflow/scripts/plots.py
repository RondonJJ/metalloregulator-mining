"""
Figuras del pipeline. SIEMPRE en SVG y PDF (nunca PNG).
- regulators_per_family_metal : heatmap familia x metal, coloreado por fraccion
  con sitio CHDE. Sugiere con que familias cuenta el genoma para sensar cada metal.
  Los perfiles se colapsan a familia; el metal es la union de los perfiles.
- regulators_per_family : reguladores por familia (con/sin sitio CHDE).
- regulators_per_genome : total por genoma.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

perfam       = pd.read_csv(snakemake.input.perfam, sep="\t")        # noqa: F821
pergen       = pd.read_csv(snakemake.input.pergen, sep="\t")        # noqa: F821
perfam_metal = pd.read_csv(snakemake.input.perfam_metal, sep="\t")  # noqa: F821

out = snakemake.output  # noqa: F821
for p in out:
    Path(p).parent.mkdir(parents=True, exist_ok=True)

C_SITE = "#1f4e79"
C_NONE = "#a6c8e8"
CMAP = LinearSegmentedColormap.from_list("blues", ["#eef4fb", "#1f4e79"])
METAL_ORDER = ["Fe", "Mn", "Zn", "Ni", "Co", "Cu", "Pb", "Fe-S", "vario", "(sin metal)"]


def save(fig, svg_path, pdf_path):
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


# ---------- Figura 1: familia x metal ----------
samples = sorted(perfam_metal["sample"].unique())
for s in samples:
    sub = perfam_metal[perfam_metal["sample"] == s]
    metals = [m for m in METAL_ORDER if m in set(sub["metal"])]
    metals += sorted(set(sub["metal"]) - set(metals))
    families = sorted(sub["family"].unique())

    frac = np.full((len(metals), len(families)), np.nan)
    annot = np.empty((len(metals), len(families)), dtype=object)
    annot[:] = ""
    for _, r in sub.iterrows():
        i, j = metals.index(r["metal"]), families.index(r["family"])
        frac[i, j] = r["frac_with_site"]
        annot[i, j] = f"{int(r['n_with_chde_site'])}/{int(r['n_regulators'])}"

    masked = np.ma.masked_invalid(frac)
    cmap = CMAP.copy()
    cmap.set_bad("#f0f0f0")

    fig, ax = plt.subplots(figsize=(0.85 * len(families) + 2.5, 0.6 * len(metals) + 2))
    im = ax.imshow(masked, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(families)))
    ax.set_xticklabels(families, rotation=45, ha="right")
    ax.set_yticks(range(len(metals)))
    ax.set_yticklabels(metals)
    for i in range(len(metals)):
        for j in range(len(families)):
            if annot[i, j]:
                val = frac[i, j]
                color = "white" if (not np.isnan(val) and val >= 0.55) else "#222222"
                ax.text(j, i, annot[i, j], ha="center", va="center", fontsize=8, color=color)
    ax.set_title(f"Reguladores por familia y metal — {s}\n(celda: con sitio / total)")
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("fraccion con sitio CHDE")
    save(fig, out.fam_metal_svg, out.fam_metal_pdf)

# ---------- Figura 2: reguladores por familia ----------
samples_f = sorted(perfam["sample"].unique())
n = len(samples_f)
fig, axes = plt.subplots(n, 1, figsize=(9, 3.2 * n), squeeze=False)
for ax, s in zip(axes[:, 0], samples_f):
    d = perfam[perfam["sample"] == s].sort_values("family")
    fams = d["family"].tolist()
    ax.bar(fams, d["n_with_chde_site"], color=C_SITE, label="con sitio CHDE")
    ax.bar(fams, d["n_without_site"], bottom=d["n_with_chde_site"], color=C_NONE, label="sin sitio")
    ax.set_title(s); ax.set_ylabel("Reguladores (hits unicos)")
    ax.tick_params(axis="x", rotation=45)
    for x, t in zip(fams, d["n_regulators"]):
        ax.text(x, t, str(int(t)), ha="center", va="bottom", fontsize=8)
axes[0, 0].legend(frameon=False)
fig.suptitle("Reguladores por familia")
fig.tight_layout()
save(fig, out.perfamily_svg, out.perfamily_pdf)

# ---------- Figura 3: total por genoma ----------
fig, ax = plt.subplots(figsize=(1.6 * max(len(pergen), 2) + 2, 4.5))
x = range(len(pergen)); w = 0.4
ax.bar([i - w/2 for i in x], pergen["n_regulators"], width=w, color=C_NONE, label="total")
ax.bar([i + w/2 for i in x], pergen["n_with_chde_site"], width=w, color=C_SITE, label="con sitio CHDE")
ax.set_xticks(list(x)); ax.set_xticklabels(pergen["sample"])
ax.set_ylabel("Reguladores (hits unicos)"); ax.set_title("Total por genoma")
for i, (a, b) in enumerate(zip(pergen["n_regulators"], pergen["n_with_chde_site"])):
    ax.text(i - w/2, a, str(int(a)), ha="center", va="bottom", fontsize=8)
    ax.text(i + w/2, b, str(int(b)), ha="center", va="bottom", fontsize=8)
ax.legend(frameon=False)
fig.tight_layout()
save(fig, out.pergenome_svg, out.pergenome_pdf)

print("[plots] SVG + PDF generados")
