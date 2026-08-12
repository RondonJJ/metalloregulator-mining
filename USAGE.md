# Guía de uso — metalloregulator-mining

Pipeline reproducible para minar metalorreguladores por familia a partir de
perfiles HMM y una base de secuencias de clusters (SSN):

**BITACORA → recolección → MetalNet2 → análisis → figuras**

Esta guía asume que empezás de cero. Está pensada para el uso propio en un
servidor con sistema de colas SGE/Grid Engine; al final hay notas para adaptarlo
a terceros.

---

## 1. Idea general y modelo de trabajo

El pipeline está orquestado con **Snakemake**. Cada etapa declara sus entradas y
salidas, y Snakemake solo re-ejecuta lo que cambió.

Hay dos máquinas con roles distintos, y conviene tenerlo claro desde el principio:

- **Tu laptop**: editás archivos (config, scripts) y usás `git`. Acá **no** corre
  nada pesado; no hace falta instalar BITACORA ni MetalNet.
- **El servidor**: acá están instalados los programas pesados y acá corrés el
  pipeline. El repo se clona/actualiza vía GitHub.

**El repositorio es solo texto** (Snakefile, scripts, configs). Los programas
(BITACORA, MetalNet2), los entornos conda y los genomas **viven en el servidor** y
se quedan ahí; el repo solo los apunta por ruta. GitHub es el puente entre laptop
y servidor: editás en la laptop → `push` → `pull` en el servidor.

---

## 2. Requisitos

### Lo que el pipeline instala solo (vía `--use-conda`)

- **BITACORA**: BLAST+ y HMMER (definidos en `workflow/envs/bitacora.yaml`).
- **Análisis**: python, pandas, matplotlib (`workflow/envs/analysis.yaml`).

### Lo que tenés que tener instalado vos (en el servidor)

- **Snakemake 7** en su propio entorno conda (ver §3). **No** sirve la 3.x.
- **BITACORA** (carpeta `Scripts/`): https://github.com/molevol-ub/bitacora
  Incluye la carpeta `Tools/` para reformatear anotaciones/proteínas.
- **MetalNet2** en su propio entorno conda, ya funcional. Esta es la dependencia
  más difícil de instalar; una vez que anda, el pipeline la usa por ruta absoluta.
- **Acceso a internet en el nodo de cómputo** (MetalNet usa la API de ColabFold
  para generar los MSA).
- GPU: **opcional**. Por defecto ESM2 corre en CPU (más lento, pero funciona).
- SGE/Grid Engine: en nuestro caso sí; el pipeline también corre local con
  `--cores` sin cola.

---

## 3. Instalar Snakemake (una sola vez, en el servidor)

Snakemake va en su **propio** entorno, separado de `base` y de `metalnet`:

```bash
conda create -n smk -c conda-forge -c bioconda 'snakemake>=7,<8' -y
conda activate smk
snakemake --version      # debe decir 7.x
```

Notas:
- Forzá los canales `conda-forge` y `bioconda`, o conda puede instalar una versión
  antiquísima (3.x) que **no** es compatible con este pipeline.
- Si conda tarda mucho resolviendo, instalá `mamba` (`conda install -n base -c
  conda-forge mamba -y`) y usá `mamba create ...`.
- Cada vez que entres al servidor a correr: `conda activate smk`.

---

## 4. Estructura del repositorio

```
metalloregulator-mining/
├── config/
│   ├── config.yaml       # rutas y parámetros (editá aquí, no en los scripts)
│   ├── samples.tsv        # genomas a minar (una fila por genoma)
│   └── profiles.tsv       # mapeo perfil HMM -> familia -> metal
├── workflow/
│   ├── Snakefile          # las 5 etapas
│   ├── scripts/
│   │   ├── runBITACORA_templated.sh    # BITACORA v1.4 parametrizado por entorno
│   │   ├── run_metalnet_colabfold.sh   # wrapper de MetalNet2
│   │   ├── collect_candidates.py
│   │   ├── parse_metalnet.py
│   │   └── plots.py
│   └── envs/
│       ├── bitacora.yaml
│       └── analysis.yaml
├── resources/             # inputs (NO se versionan): query_db, genomas
└── results/               # salidas (NO se versionan)
```

`resources/` y `results/` están en `.gitignore`: los genomas y las salidas no se
suben a GitHub.

---

## 5. Preparación de datos (paso manual previo)

Antes de correr el pipeline hay que preparar y **verificar** los archivos de cada
genoma. Esto es manual a propósito, porque depende del genoma.

### 5.1 Query Database (perfiles de la SSN)

La carpeta `query_db` (BITACORA la llama FPDB) tiene, por cada perfil, dos archivos:
`<PERFIL>_db.fasta` y `<PERFIL>_db.hmm`. Ejemplo:

```
databases/
├── BsCzrA_db.fasta      BsCzrA_db.hmm
├── mtNmtR_db.fasta      mtNmtR_db.hmm
├── ecFur_db.fasta       ecFur_db.hmm
└── ...
```

### 5.2 Archivos por genoma

BITACORA en modo genoma necesita **tres** archivos por genoma:

- genoma en FASTA (`.fna` / `.fasta`)
- GFF3/GTF **reformateado** (`.gff3`)
- proteínas predichas del GFF (FASTA)

Si el GFF o las proteínas no están en el formato que BITACORA espera, usá los
scripts de la carpeta **`Tools/`** de BITACORA para reformatearlos/inferirlos.
Chequealo antes de lanzar: un GFF mal formado hace fallar la etapa de BITACORA en
`check_data`.

### 5.3 `config/samples.tsv`

Una fila por genoma, con las rutas **absolutas** en el servidor (pueden estar
fuera del repo):

```
sample      genome_fasta               genome_gff                proteins_fasta            evalue
MtubH37Rv   /ruta/.../H37Rv.fna        /ruta/.../H37Rv.gff3      /ruta/.../H37Rv.faa       1e-5
Mavi_homi   /ruta/.../Mavi_hom.fna     /ruta/.../Mavi_hom.gff    /ruta/.../Mavi_hom.faa    1e-5
```

(Separado por tabuladores. Cada genoma puede tener su propio E-value.)

### 5.4 `config/profiles.tsv`

Mapea cada perfil HMM a su familia biológica y su metal. Es lo que colapsa los
perfiles a familia (p. ej. BsCzrA + mtNmtR → ArsR) y da las etiquetas de metal
de las figuras. Debe cubrir **todos** los perfiles de tu `query_db`:

```
profile         family   metal
BsCzrA          ArsR     Zn/Co
mtNmtR          ArsR     Ni/Co
ecFur           Fur      Fe
ecZur           Fur      Zn
...
```

Si un perfil aparece en los datos pero no está acá, el análisis lo procesa igual
pero cae en un fallback (usa el nombre del perfil como familia, sin metal) y lo
avisa por consola.

---

## 6. Configuración (`config/config.yaml`)

Todo lo ajustable vive acá; no hay nada hardcodeado en los scripts.

**BITACORA**
- `bitacora_scriptdir`: carpeta `Scripts/` de tu instalación de BITACORA.
- `query_db`: carpeta con las FPDB (§5.1).
- `bitacora_gemomap`: jar de GeMoMa (solo si `bitacora_gemoma: "T"`).
- `bitacora_threads`, `bitacora_gemoma`, `bitacora_maxintron`,
  `bitacora_genomicblastp`, `bitacora_addfilter`, `bitacora_filterlength`,
  `bitacora_retainnonfilter`, `bitacora_clean`: parámetros de BITACORA.
- `bitacora_keep_suffix`: sufijo de los archivos de BITACORA que se usan
  (`_genomic_and_annotated_proteins_trimmed_idseqsclustered`). No cambiar salvo
  que BITACORA cambie sus nombres de salida.

**MetalNet2**
- `metalnet_python`: ruta absoluta al `python` de tu entorno `metalnet`.
- `metalnet_script`: ruta a `run_prediction_workflow.py` de MetalNet2.
- `metalnet_msa_source`: `colabfold` (necesita internet en el nodo).
- `metalnet_use_gpu`: `false` (CPU) o `true` (agrega `--cuda 1`).
- `metalnet_threads`.

**Análisis / figuras**
- `profiles_map`: `config/profiles.tsv`.
- `metalnet_min_prob`: umbral de probabilidad para contar un par (default `0.0`).
- `metalnet_require_graph_filter`: si `true`, un "sitio CHDE" exige
  `filter_by_graph == 1` (el veredicto de MetalNet). **Dejar en `true`** — es el
  criterio correcto; con `false` cuenta cualquier par candidato y todo da 100%.

---

## 7. Cómo correr

Siempre desde la raíz del repo, con el entorno `smk` activo.

### 7.1 Validar el plan (dry-run, no ejecuta nada)

```bash
conda activate smk
cd /ruta/al/repo/metalloregulator-mining
snakemake -n
```

Debería mostrar la cadena `bitacora → collect_candidates → metalnet →
parse_metalnet → plots`. El dry-run chequea que existan los archivos de entrada,
así que primero tienen que estar cargados `samples.tsv`, `query_db` y los genomas.

### 7.2 Ejecutar (local, sin cola)

```bash
snakemake --use-conda --conda-frontend conda --cores 10
```

- `--use-conda`: crea los entornos de BITACORA y análisis automáticamente.
- `--conda-frontend conda`: necesario si no tenés `mamba` instalado.
- Sin un target al final, corre hasta el objetivo final (`rule all`).

Para correr **una sola etapa**, pedí su salida:

```bash
snakemake --use-conda --conda-frontend conda --cores 10 results/bitacora/MtubH37Rv
```

### 7.3 Ejecutar en SGE (por cola)

Para tareas pesadas, envolvé Snakemake en un script de `qsub`. Ejemplo (ajustá
cola y `-pe` a tu servidor):

```bash
cat > run_smk.sh << 'EOF'
#!/bin/bash -x
#$ -cwd
#$ -V
#$ -q h14.q
#$ -pe ompi511h14 10
#$ -N smk_pipeline
source ~/miniconda3/etc/profile.d/conda.sh
conda activate smk
cd /ruta/al/repo/metalloregulator-mining
snakemake --use-conda --conda-frontend conda --cores 10
EOF

qsub run_smk.sh
```

El `#$ -V` es importante: exporta tu entorno al trabajo, para que el
`conda activate` funcione en el nodo.

**Para la etapa de MetalNet** (necesita internet), usá un `qsub` con la cola/nodo
que tenga salida a internet, pidiendo su salida:
`snakemake ... results/metalnet/output/pred_pairs.tsv`.

### 7.4 Estrategia recomendada por tandas

1. **BITACORA** de todos los genomas (cómputo puro, sin internet):
   `snakemake ... results/bitacora/<genoma>` para cada uno.
2. **Recolección + MetalNet + análisis + figuras** (MetalNet es la etapa larga y
   necesita internet): un `qsub` pidiendo el objetivo final.

Seguí el progreso con `tail -f logs/<etapa>.log` y `qstat`.

---

## 8. Salidas

### Tablas (`results/tables/`)
- `metalnet_pairs_annotated.tsv`: un renglón por par de residuos predicho, con
  genoma/perfil/familia.
- `regulators_per_protein.tsv`: un renglón por proteína (regulador único).
- `regulators_per_profile.tsv`: por perfil (comparable a figuras por perfil).
- `regulators_per_family.tsv`: por familia (perfiles colapsados).
- `regulators_per_family_metal.tsv`: por familia y metal.
- `regulators_per_genome.tsv`: totales por genoma.

### Figuras (`results/figures/`, SVG **y** PDF)
- `regulators_per_family_metal.*`: heatmap familia × metal, un panel por genoma.
- `regulators_per_family.*`: reguladores por familia (con/sin sitio CHDE).
- `regulators_per_genome.*`: total por genoma.

---

## 9. Conceptos clave del conteo

- **Perfil vs familia**: un "perfil" es una FPDB/HMM (BsCzrA, mtNmtR…); una
  "familia" es el grupo biológico (ArsR…). Algunas familias tienen 2 perfiles.
  El mapeo está en `profiles.tsv`.
- **Hits únicos**: un regulador es una proteína única por genoma. Si la misma
  proteína la detectan 2 perfiles de la misma familia, cuenta **una** vez.
- **Sitio CHDE**: una proteína "tiene sitio" si MetalNet predice ≥1 par que pasa
  el filtro de grafo (`filter_by_graph == 1`). Es el veredicto del propio modelo.

---

## 10. Agregar más genomas

1. Prepará los 3 archivos del genoma nuevo (§5.2) y agregá su fila a
   `config/samples.tsv`.
2. Si el genoma trae perfiles nuevos, agregalos a `config/profiles.tsv`.
3. Corré BITACORA del genoma nuevo, luego el resto.

**Importante**: al cambiar `samples.tsv`, la recolección debería rehacerse, pero
por cómo Snakemake maneja los timestamps a veces no lo detecta solo. Si ves que el
`candidates.fasta` no incluye el genoma nuevo, forzá el rehacer borrando desde la
recolección para abajo (sin tocar BITACORA, que es lo caro):

```bash
rm -rf results/metalnet results/tables results/figures
# NO borres results/bitacora/
snakemake -n     # verificá que NO reaparezca 'bitacora'
```

Comprobá los genomas del candidates tras la recolección:

```bash
grep ">" results/metalnet/input/candidates.fasta | sed 's/>//; s/__.*//' | sort -u
```

---

## 11. Solución de problemas (errores reales que vimos)

**`snakemake: command not found`**
No activaste el entorno o no está instalado. `conda activate smk` (§3).

**`snakemake --version` dice 3.x**
Versión demasiado vieja. Reinstalá con canales correctos:
`conda create -n smk -c conda-forge -c bioconda 'snakemake>=7,<8'` (§3).

**`The 'mamba' command is not available`**
Agregá `--conda-frontend conda` al comando, o instalá mamba.

**`--use-conda` descarga BLAST/HMMER aunque ya los tengo**
Es a propósito: crea entornos aislados para reproducibilidad (una sola vez, se
cachea). Si preferís usar los tuyos del sistema, corré **sin** `--use-conda`
(deben estar en el PATH: `makeblastdb`, `tblastn`, `blastp`, `hmmsearch`,
`hmmbuild`).

**MetalNet: `ModuleNotFoundError: No module named 'absl'`**
Un subproceso de MetalNet no heredaba el entorno `metalnet`. El wrapper
`run_metalnet_colabfold.sh` ya lo activa (`conda activate metalnet`). Verificá que
la ruta a tu conda en ese script sea correcta (`~/miniconda3/etc/profile.d/conda.sh`).

**MetalNet: `sh: 1: <familia>: not found` y falta `coevo.csv` / archivos `.uni.a3m`**
El `|` en los headers rompe el shell. El pipeline usa `__` como separador
(`sample__perfil__proteina`); si ves esto, revisá que `collect_candidates.py` use
`__`, no `|`, y regenerá `candidates.fasta`.

**Todas las familias dan 100% / N/N (todo tiene sitio)**
`metalnet_require_graph_filter` está en `false`. Ponelo en `true` (§6) y rehacé
el análisis.

**El heatmap muestra `(sin metal)` en todas las filas**
`config/profiles.tsv` es una versión vieja sin la columna `metal`, o está
desincronizado con el servidor. Actualizalo (§5.4).

**BITACORA muere en `check_data`**
Formato de GFF o proteínas. Reformateá con la carpeta `Tools/` de BITACORA (§5.2).

**Después de agregar genomas, la tabla sigue mostrando solo algunos**
La recolección no se rehízo. Ver §10 (borrar desde `results/metalnet` para abajo).

---

## 12. Reproducibilidad y límites

- **ColabFold no es determinista en el tiempo**: los MSA vienen por API y su base
  cambia, así que dos corridas separadas por meses pueden dar números finos
  distintos. El patrón se mantiene; los valores exactos pueden moverse. Es
  inherente al modo `colabfold`.
- Para dejar registro del entorno de MetalNet:
  `conda env export --no-builds -n metalnet > workflow/envs/metalnet.lock.yaml`.
- **El vínculo familia→metal** (`profiles.tsv`) es curado por vos, no una
  predicción. Para el metal predicho por sitio, MetalNet trae
  `predict_metal_type.py` (se podría enganchar como capa extra).
- El eje de metales y la vecindad genómica de figuras tipo panel C/D son
  curación/otra capa; el pipeline automatiza los conteos y el sitio CHDE.

---

## 13. Notas para publicar / terceros (pendiente)

- **MetalNet en contenedor** (Apptainer/Singularity): es lo que elimina el dolor
  de instalación y hace el pipeline portable de verdad.
- **Perfil de SGE** para que Snakemake mande cada regla como un trabajo aparte.
- Hacer `config/samples.tsv` dependencia explícita de `collect_candidates` para
  que agregar genomas invalide el `candidates.fasta` automáticamente.
- Traducir esta guía al inglés y agregar una licencia.

---

## 14. Ciclo de trabajo con git (laptop ↔ servidor)

```bash
# --- laptop: editás y subís ---
git add -A
git commit -m "mensaje"
git push

# --- servidor: bajás y corrés ---
git pull
conda activate smk
snakemake -n
```

Para evitar desincronizaciones (nos pasó con `profiles.tsv`), cuando toques varios
archivos conviene sincronizar todo el repo por git en vez de copiar archivos
sueltos.
