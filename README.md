# metalloregulator-mining

Downstream pipeline to **mine, confirm and characterise bacterial metalloregulators**
from genome annotations, starting from candidate transcription factors recovered
by SSN-derived HMM profiles + [BITACORA](https://github.com/molevol-ub/bitacora),
adding **coevolution-based metal-site prediction** ([MetalNet](https://github.com/wangchulab/MetalNet))
and **genomic-context** evidence.

This repository holds the code, result tables and figures for the genome-mining
part of the study on *Mycobacterium tuberculosis* H37Rv, *M. avium* subsp.
*hominissuis*, *Vibrio cholerae* RFB16 and *V. vulnificus* NBRC 15645.

> It does **not** re-implement BITACORA, HMMER, HHblits/ColabFold or MetalNet —
> those are external tools. This repo covers everything **after** BITACORA:
> consolidating candidates, running MetalNet, predicting metal-coordination
> sites, and cross-referencing with genomic neighbourhood.

---

## Pipeline overview

```
SSN clusters ──▶ HMM profiles ──▶ BITACORA search        (external; not in this repo)
                                        │
                                        ▼
                        [ prepare_metalnet_input.py ]     consolidate + dedup candidates
                                        │  candidates.fasta + membership.tsv
                                        ▼
                        [ workflow/*.sh (SGE) ]           MSAs + MetalNet2 (ColabFold or HHblits)
                                        │  pred_pairs.tsv
                                        ▼
                        [ summarize_metalnet.py ]         per-protein CHED metal-site calls
                                        │  site_summary.tsv / confident_pairs.tsv
                                        ▼
                        [ neighborhood_analysis.py ]      metal-related genes in the vicinity
                                        │  neighborhood_table.tsv + master table
                                        ▼
                        [ family_genome_tables.py ]       counts by family × genome
                        [ figures/plot_metal_by_genome.py ]  metal-coverage figure
```

---

## Repository layout

```
scripts/     Python analysis (no external deps beyond the std lib)
  prepare_metalnet_input.py   BITACORA FASTAs -> deduplicated candidates.fasta + membership.tsv
  summarize_metalnet.py       MetalNet pred_pairs.tsv -> per-protein site summary
  neighborhood_analysis.py    genome GFF + candidates -> genomic-neighbourhood / "effective" table
  family_genome_tables.py     site/total counts per family and genome

workflow/    Cluster (SGE/qsub) scripts to run MSAs + MetalNet2
  run_metalnet_colabfold_PILOT.sh   3-sequence test (ColabFold MSA route)
  run_metalnet_colabfold.sh         full run (ColabFold MSA route, recommended)
  00_split.sh 01_make_msas.sh 02_build_manifest.sh 03_run_metalnet.sh   HHblits route (plan B)

results/     Result tables produced for the paper
  candidates.fasta, membership.tsv
  full_run/  site_summary, confident_pairs, master_candidates(_completed),
             neighborhood_table, counts_by_genome_family, family/genome_success, genome_family_matrix
  pilot_run/ site_summary, confident_pairs

figures/     Editable SVG figures + their generators
  metal_by_genome.svg          metal × genome coverage panel
  plot_metal_by_genome.py      generator for the above
  Fig5_with_panelE.svg         full Figure 5 (A–D) + panel E merged
  build_fig5_panelE.py         merges panel E under an existing Fig5 base SVG

docs/
  INSTRUCTIVO.md               step-by-step how-to (Spanish)
```

---

## Requirements

- **Python ≥ 3.8** — analysis scripts use only the standard library (no pip installs needed).
- **External tools** (for the `workflow/` stage, on your cluster):
  MetalNet2, and either ColabFold (remote MMseqs2 API, needs internet on the node)
  or HHblits + a UniRef30/BFD database. SGE/`qsub` for job submission.
- Figures are plain SVG; edit in Inkscape/Illustrator or render with any SVG tool.

See `requirements.txt` for optional Python helpers.

---

## Quick start (analysis only, from BITACORA output)

```bash
# 1) consolidate BITACORA candidates (<root>/<genome>/<family>/*.fasta)
python scripts/prepare_metalnet_input.py -i Results_Bita_-5 -o metalnet_input

# 2) run MetalNet2 on the cluster (see docs/INSTRUCTIVO.md and workflow/)
#    -> produces output/pred_pairs.tsv

# 3) summarise metal-coordination sites
python scripts/summarize_metalnet.py -p output/pred_pairs.tsv \
    -m metalnet_input/membership.tsv \
    -o results/site_summary.tsv --pairs-out results/confident_pairs.tsv

# 4) genomic neighbourhood (needs whole-genome GFFs) — edit paths inside the script
python scripts/neighborhood_analysis.py

# 5) family × genome counts and the metal-coverage figure
python scripts/family_genome_tables.py
python figures/plot_metal_by_genome.py -o figures/metal_by_genome.svg
```

Full details, cluster commands and interpretation notes are in
[`docs/INSTRUCTIVO.md`](docs/INSTRUCTIVO.md).

---

## Interpreting the results

- `confident_pairs.tsv` — the high-confidence coevolving CHED residue pairs
  (`filter_by_graph == 1` in MetalNet output); positions are **1-based**.
- `site_summary.tsv` — one row per protein: whether a metal site was predicted,
  which residues, CHED composition.
- `master_candidates_completed.tsv` — one row per (genome, accession) with the
  MetalNet call **and** the genomic-neighbourhood flags, plus an `effective_*`
  column (predicted site **and** an adjacent metal-homeostasis gene).

**Caveats.** A predicted site supports metalloregulator assignment but the
absence of one is not proof of non-binding (divergent/carboxylate sites,
global regulators such as Fur/Zur). Homology and neighbourhood **approximate**
the sensed metal; they do not resolve specificity. See the paper for details.

---

## Citing

If you use this code, please cite the associated paper (see `CITATION.cff`)
and the external tools it builds on: BITACORA, MetalNet, HMMER, and
ColabFold/HHblits.

## License

Code released under the MIT License (see `LICENSE`).
