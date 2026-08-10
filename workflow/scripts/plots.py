"""
Figuras basicas a partir de la tabla anotada. Ampliar segun tus graficos previos.
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

table_path = snakemake.input.table   # noqa: F821
out_png    = snakemake.output[0]     # noqa: F821

df = pd.read_csv(table_path, sep="\t")
Path(out_png).parent.mkdir(parents=True, exist_ok=True)

# >>> EDITAR (3/3): adapta este grafico (y agrega los tuyos) a tus columnas <<<
if "family" in df.columns and len(df):
    counts = df.groupby("family").size().sort_values(ascending=False)
    ax = counts.plot(kind="bar", figsize=(8, 5))
    ax.set_ylabel("Pares / sitios predichos")
    ax.set_xlabel("Familia")
    ax.set_title("Sitios de union a metal predichos por familia")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
else:
    # placeholder para no romper el DAG si aun no hay datos
    plt.figure()
    plt.text(0.5, 0.5, "Sin datos de familia", ha="center", va="center")
    plt.axis("off")
    plt.savefig(out_png, dpi=200)

print(f"[plots] figura -> {out_png}")
