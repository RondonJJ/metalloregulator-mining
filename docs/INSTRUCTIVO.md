# Instructivo — Análisis de metaloproteínas con MetalNet2 en el clúster

**Pipeline:** candidatas de BITACORA → MetalNet2 → sitios de unión a metal
**Clúster:** SGE (cola `h14.q`, entorno paralelo `ompi511h14`)
**Última actualización de referencia:** corrida piloto validada (ruta ColabFold)

---

## 0. Qué hace este pipeline (en una frase)

Toma las proteínas candidatas que salieron de BITACORA (organizadas por
organismo y familia de regulador), las consolida en una única entrada sin
redundancia, y corre **MetalNet2** para predecir **qué candidatas son
metaloproteínas y qué residuos forman el sitio de unión a metal**. Al final
genera un resumen legible por proteína.

> **Ojo con lo que NO responde:** MetalNet dice *si* hay sitio y *dónde* está,
> pero **no** dice *qué* metal une (ver sección 8). La familia de BITACORA es
> una *hipótesis* de metal, no una prueba.

---

## 1. Requisitos y rutas fijas

| Cosa | Valor en tu instalación |
|---|---|
| Entorno conda | `metalnet` |
| Python del entorno | `/users-d2/j.rondon/miniconda3/envs/metalnet/bin/python` |
| Motor de MetalNet2 | `/users-d2/j.rondon/programas/MetalNet2/model/scripts/run_prediction_workflow.py` |
| Carpeta de la app | `~/programas/MetalNet2/app/` |
| Cola / PE | `h14.q` / `ompi511h14` |

**Dos rutas para generar los MSAs:**

- **ColabFold (la que usamos, recomendada):** MetalNet2 pide los MSAs a la API
  de ColabFold. Requiere **internet en el nodo de cómputo** (ya confirmado en
  h14). No necesita bases de datos locales.
- **HHblits (plan B):** genera los MSAs localmente. Requiere una **base de datos
  tipo UniRef30/BFD** en disco. Úsala solo si la API está caída o te limita.

---

## 2. Estructura de directorios

Trabajá en una carpeta nueva dentro de `app/`, por ejemplo `bita_run/`:

```
~/programas/MetalNet2/app/bita_run/
├── input/
│   └── candidates.fasta        # entrada consolidada (161 secuencias únicas)
├── membership.tsv              # mapea seq_id -> organismo/familia
├── scripts/                    # todos los scripts de abajo
└── output/                     # lo genera MetalNet (MSAs, ESM2, predicciones)
```

---

## 3. Los archivos que armamos

| Archivo | Para qué sirve |
|---|---|
| `prepare_metalnet_input.py` | Consolida los FASTA de BITACORA en `candidates.fasta` + `membership.tsv` (dedup + IDs limpios). |
| `scripts/run_metalnet_colabfold_PILOT.sh` | Corrida de prueba (3 secuencias) para validar la API antes de lanzar todo. |
| `scripts/run_metalnet_colabfold.sh` | **Corrida completa por ColabFold** (ruta principal). |
| `summarize_metalnet.py` | Convierte `pred_pairs.tsv` en un resumen por proteína (`site_summary.tsv`). |
| `scripts/00_split.sh` | (HHblits) Divide el FASTA en uno por secuencia + lista para el array. |
| `scripts/01_make_msas.sh` | (HHblits) Array job que genera un MSA por secuencia. |
| `scripts/02_build_manifest.sh` | (HHblits) Arma el manifiesto `msa_files.tsv`. |
| `scripts/03_run_metalnet.sh` | (HHblits) Corrida de MetalNet en modo `user_defined`. |

---

## 4. Flujo de trabajo — ruta ColabFold (recomendada)

### Paso 0 — Consolidar la salida de BITACORA (una sola vez)

Desde donde tengas los resultados de BITACORA (carpeta con estructura
`<organismo>/<familia>/*.fasta`):

```bash
python prepare_metalnet_input.py -i Results_Bita_-5 -o metalnet_input
```

Genera `metalnet_input/candidates.fasta` (secuencias únicas, IDs seguros) y
`metalnet_input/membership.tsv`. Copiá el FASTA a `bita_run/input/` y el
`membership.tsv` a `bita_run/`.

> Por qué este paso: la misma proteína suele ser recuperada por varias
> familias/organismos. Correr duplicados es tirar cómputo, y los IDs repetidos
> harían colisionar los archivos de salida (MetalNet nombra todo por `seq_id`).

### Paso 1 — Piloto (validar la API)

```bash
cd ~/programas/MetalNet2/app/bita_run/
qsub scripts/run_metalnet_colabfold_PILOT.sh
```

Corre 3 secuencias. Si al terminar aparece `output_pilot/pred_pairs.tsv` y se
pobló `output_pilot/msa/`, la conexión a ColabFold funciona. **No sigas si el
piloto falla.**

### Paso 2 — Corrida completa

```bash
qsub scripts/run_metalnet_colabfold.sh
```

Procesa las 161 en una sola pasada. Tiempo estimado: **~5–8 h** (dominado por
la descarga de MSAs). Dejá que corra de una noche. Resultados en
`output/pred_pairs.tsv` + un `.gv.pdf` por proteína con sitio.

Monitoreo:

```bash
qstat                       # estado del job
tail -f metalnet_bita_cf.o<JOBID>   # log de salida
```

### Paso 3 — Resumen legible

```bash
python summarize_metalnet.py \
    -p output/pred_pairs.tsv \
    -m membership.tsv \
    -o output/site_summary.tsv \
    --pairs-out output/confident_pairs.tsv
```

- `site_summary.tsv`: una fila por proteína — si tiene sitio, residuos, composición CHED, probabilidad, organismo y familia.
- `confident_pairs.tsv`: los pares de alta confianza en detalle, con posiciones 1-based (listas para mapear sobre estructura).

---

## 5. Ruta alternativa — HHblits (plan B)

Solo si no podés usar ColabFold. Reemplaza los pasos 1–2 de arriba:

```bash
# 0) dividir e indexar (nodo de login, rápido)
bash scripts/00_split.sh          # imprime el rango del array, ej. -t 1-161

# 1) editar scripts/01_make_msas.sh:
#    - poner la ruta de la DB en la variable DB=...
#    - ajustar el "#$ -t 1-161" al número que imprimió 00_split.sh
qsub scripts/01_make_msas.sh      # array job: un MSA por secuencia

# 2) armar el manifiesto cuando terminen todos
bash scripts/02_build_manifest.sh # genera output/msa_files.tsv

# 3) correr MetalNet en modo user_defined
qsub scripts/03_run_metalnet.sh
```

Después seguís con el **Paso 3** (resumen) igual que en la ruta ColabFold.

---

## 6. Cómo leer los resultados

### `pred_pairs.tsv` (salida cruda de MetalNet)

Cada fila es un par de residuos CHED evaluado. Columnas clave:

- `resi_1`, `resi_seq_posi_1`, `resi_2`, `resi_seq_posi_2`: los dos residuos y sus posiciones.
- `prob`: probabilidad del modelo de que el par una metal.
- **`filter_by_graph`**: `1` = el par pasó el filtro de red/grafo → **sitio de alta confianza**. `0` = descartado. **Esta es la columna que importa.**

> ⚠️ **Numeración:** en `pred_pairs.tsv` las posiciones son **0-based**. El
> residuo real (el que ves en el `.gv` y en la proteína) es **posición + 1**.
> El script `summarize_metalnet.py` ya hace esta corrección; `confident_pairs.tsv`
> sale directamente en 1-based.

### `site_summary.tsv` (lo que vas a mirar la mayoría del tiempo)

Una fila por proteína, ordenada por confianza: `has_site` (yes/no), residuos del
sitio (ej. `H108,H111,D115`), composición (ej. `3H,1E,1D`), `max_prob`, y el
organismo/familia de BITACORA.

### `<seq_id>.gv` / `<seq_id>.gv.pdf`

Visualización del sitio: nodos = residuos (numeración real 1-based), aristas =
pares coevolucionantes. Un grafo vacío = no se predijo sitio.

---

## 7. Ajustes útiles

### Subir de 4 a 12 threads

En `run_metalnet_colabfold.sh`, cambiá:

```bash
#$ -pe ompi511h14 12          # antes 4

export OMP_NUM_THREADS=12     # antes 4
export MKL_NUM_THREADS=12
export OPENBLAS_NUM_THREADS=12
```

**Salvedades:**
- En la ruta ColabFold el grueso del tiempo es la **descarga de MSAs (red)**, que los threads NO aceleran. Los 12 hilos solo ayudan a ESM2 (CPU) y al cálculo de coevolución.
- Asegurate de que el PE `ompi511h14` te dé los 12 núcleos en **un solo nodo** (los threads son de memoria compartida; repartidos entre nodos no sirven).
- En la ruta HHblits sí conviene subir `-cpu` y `CPU` a 12 (HHblits es intensivo en CPU de verdad).

### GPU

Si el nodo tiene GPU, descomentá `--cuda 1` en el script para acelerar ESM2.
En CPU corre igual, solo más lento.

### Reanudar una corrida cortada

Los scripts dejan `--keep_inter_files` activo: los MSAs quedan en `output/msa/`.
Si el job se cae, no perdés lo ya calculado. Podés incluso armar el manifiesto
desde esos `.a3m` (`02_build_manifest.sh`) y terminar en modo `user_defined`.

### Buena convivencia con la API

Lanzá la corrida completa como **un job único**, no muchas copias en paralelo,
para no saturar la API pública de ColabFold.

---

## 8. Límites de interpretación (leer antes de sacar conclusiones)

- MetalNet responde **"¿es metaloproteína y dónde está el sitio?"**, no **"¿qué metal?"**. La identidad del metal solo la infiere en casos donde la topología coincide con sitios conocidos, y en el piloto no la emitió.
- El mismo conjunto de ligandos His/Asp/Glu/Cys puede coordinar Zn, Ni, Co, Mn, Fe, Cd… La identidad la definen geometría fina y residuos de segunda esfera, que MetalNet no modela.
- La **familia de BITACORA es una hipótesis** del metal cognado del query, no una asignación. (Ejemplo real: `WP_003877217_1` fue recuperada por BsCzrA *y* mtNmtR a la vez.)
- Para asignar metal harían falta: estructura (AlphaFold) + geometría del sitio, ubicación filogenética con sensores caracterizados, mapeo sobre motivos conocidos de cada familia, y en última instancia validación experimental.

---

## 9. Chuleta rápida

```bash
# --- una sola vez: consolidar BITACORA ---
python prepare_metalnet_input.py -i Results_Bita_-5 -o metalnet_input
# copiar metalnet_input/candidates.fasta -> bita_run/input/
# copiar metalnet_input/membership.tsv   -> bita_run/

cd ~/programas/MetalNet2/app/bita_run/

# --- validar y correr ---
qsub scripts/run_metalnet_colabfold_PILOT.sh    # 1) piloto (3 seqs)
qsub scripts/run_metalnet_colabfold.sh          # 2) completo (161 seqs, ~5-8 h)

# --- resumir ---
python summarize_metalnet.py -p output/pred_pairs.tsv -m membership.tsv \
    -o output/site_summary.tsv --pairs-out output/confident_pairs.tsv

# resultado principal: output/site_summary.tsv
# clave de lectura: filter_by_graph==1 = sitio de alta confianza
```
