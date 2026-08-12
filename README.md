# metalloregulator-mining

Pipeline reproducible para minar metalorreguladores por familia
(**BITACORA → MetalNet2 → análisis → figuras**), a partir de perfiles HMM
y una base de secuencias de clusters obtenidos por SSN.

## Idea de diseño

Cada etapa declara sus entradas y salidas; Snakemake solo re-corre lo que cambió.
El mismo `Snakefile` funciona **local** (sin cola) o en tu **server SGE** (con un
perfil). El paso pesado (MetalNet2) está aislado como una *costura*: se puede
correr aparte y el pipeline retoma solo.

```
config/config.yaml     -> rutas y parametros (edita aqui, no en los scripts)
config/samples.tsv     -> genomas a minar
workflow/Snakefile     -> las 5 etapas
workflow/scripts/      -> tus scripts (incluye run_metalnet_colabfold.sh)
workflow/envs/         -> entornos conda por etapa (versiones fijadas)
resources/             -> inputs (query_db, genomas)  [no se versionan]
results/               -> salidas                       [no se versionan]
```

## Paso previo MANUAL (antes de correr el pipeline)

Preparar y **verificar** los archivos de cada genoma con `Tools/` de BITACORA
(formato de GFF, inferir/reformatear proteínas, coincidencia de IDs entre GFF y
proteínas). Cuando los 3 archivos por muestra estén listos, cárgalos en
`config/samples.tsv`.

## Antes de correr: 3 cosas a editar

1. **`workflow/Snakefile`, regla `bitacora`** — la invocación real de
   `runBITACORA_genome_mode.sh` (editando las variables del script maestro).
   La salida debe caer dentro de `results/bitacora/{sample}`.
2. **`workflow/scripts/parse_metalnet.py`** — el nombre real de la columna de ID
   de proteína en `pred_pairs.tsv` (hay autodetección, confírmala).
3. **`workflow/scripts/plots.py`** — adapta/añade tus gráficos previos.

En `config/config.yaml` revisa además las rutas de tu entorno de MetalNet2.

## Cómo correr

Local (una laptop, sin cola):

```bash
snakemake --use-conda --cores 4
```

En tu server SGE (la cola es un perfil opcional, no un requisito del pipeline):

```bash
snakemake --use-conda --profile workflow/profiles/sge
```

Solo hasta preparar el input de MetalNet2 (para quien no tenga GPU/internet):

```bash
snakemake --use-conda --cores 4 results/metalnet/input/candidates.fasta
```

## MetalNet2: lo que hay que saber

- Usa **tu entorno conda ya funcional** por ruta absoluta (ver `config.yaml`),
  no se recrea desde un yaml. El arreglo definitivo de portabilidad es empaquetarlo
  en un contenedor **Apptainer/Singularity** (construir la imagen una vez y compartirla).
- `--msa_source colabfold` **necesita internet en el nodo** y es no-determinista
  (la base remota cambia con el tiempo). `--keep_inter_files` cachea los MSA.
- **GPU es opcional**: por defecto corre ESM2 en CPU (más lento). Para acelerar,
  poné `metalnet_use_gpu: true` en `config.yaml`.

## Para publicar / dejar registro

```bash
conda env export --no-builds -n metalnet > workflow/envs/metalnet.lock.yaml  # congela tu env
snakemake --report report.html                                              # DAG + versiones
```
