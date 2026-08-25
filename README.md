# SQANTI-MOD

Pipeline for the integration of long-reads modification calls in uBAM files into the SQANTI classification file of a reconstructed transcriptome.

It consists of three basic steps:

1. `Genomic coverage filtering`: We use modkit, ONT's tool for working with modified bases, to generate a BED file of the modified sites and we filter out calls with low genomic coverage based on a specified threshold(with default parameters we filter out sites with < 5 reads) 
2. `Read to transcript projection`: proyecta las modificaciones desde las coordenadas genómicas a las coordenadas del transcrito usando el BAM, el GTF y la asociación lectura-transcrito de IsoQuant.
3. `collapse_to_tx_and_merge.py`: agrupa las posiciones modificadas por transcrito, aplica filtros de cobertura y ocupancia, y añade la información a `classification.txt` de SQANTI3.

## Requisitos

- Linux o un entorno HPC con SLURM para ejecutar los scripts `.sh`.
- Python 3.
- `pysam`, `pandas` y `numpy`.
- [`modkit`](https://nanoporetech.github.io/modkit/) disponible en el `PATH` o indicado mediante una ruta absoluta.
- BAM alineado e indexado (`sample.bam` y `sample.bam.bai`) con tags de modificaciones.
- Genoma de referencia en FASTA.
- Anotación GTF con atributos `transcript_id`.
- Archivo de asociación lectura-transcrito generado por IsoQuant, en TSV o TSV comprimido con gzip.
- Archivo `classification.txt` generado por SQANTI3.

Instalación de las dependencias Python en el entorno activo:

```bash
python -m pip install pysam pandas numpy
```

## Configuración

Edita [`config.sh`](config.sh) con las rutas de tu proyecto. La plantilla [`config.ini.example`](config.ini.example) contiene además los nombres de los parámetros usados por el pipeline.

Los parámetros principales son:

| Parámetro | Valor predeterminado | Descripción |
|---|---:|---|
| `MIN_GENOMIC_COV` | `5` | Cobertura mínima para conservar un sitio en el BED. |
| `PROB_LIM` | `0.95` | Probabilidad mínima de modificación. |
| `MIN_TX_COV` | `20` | Cobertura mínima de una posición en el transcrito. |
| `MIN_TX_OCC` | `0.2` | Ocupancia mínima de una modificación. |
| `MODS` | `m5C m6A 2OmeC pseU 2OmeU inosine 2OmeA` | Modificaciones solicitadas a `modkit`. |

Antes de ejecutar el segundo y el tercer paso, define las variables específicas de cada muestra (`BED`, `BAM`, `GTF`, `ASSOC`, `OUT_TSV` y `SQANTI_CLASS`) en el entorno o en el script que envuelva cada trabajo.

## Ejecución por etapas

Desde el directorio del proyecto:

```bash
source config.sh
mkdir -p "$BED_OUTDIR" "$TSV_OUTDIR" "$CLASS_OUTDIR" "$LOG_DIR"
```

### 1. Generar y filtrar BED

El script recorre los BAM que coinciden con `BAM_DIR` y `BAM_SUFFIX`:

```bash
sbatch 01_bedmethyl.sh
```

Produce:

```text
data/beds/<muestra>.bed
data/beds/<muestra>_filtered.bed
logs/pileup_<muestra>.txt
```

### 2. Convertir modificaciones a coordenadas de transcrito

Configura las variables de la muestra y ejecuta:

```bash
sbatch 02_mods_to_CIGAR_transcripts.sh
```

El resultado es un TSV por lectura con, entre otras, las columnas `transcript_id`, `transcript_length`, `tx_start`, `tx_end` y `positions_<modificación>`.

### 3. Colapsar por transcrito y combinar con SQANTI3

```bash
sbatch 03_collapse_to_tx_and_merge.sh
```

El archivo final es un TSV que conserva las columnas de SQANTI3 y añade, para cada modificación, las posiciones, coberturas y ocupancias observadas.

## Ejecución directa para una muestra

También existe un orquestador en Python que ejecuta la generación del BED y la extracción transcriptómica:

```bash
python sqanti-mod.py \
  --bam /ruta/sample.bam \
  --gtf /ruta/sample.gtf \
  --classification /ruta/classification.txt \
  --reference /ruta/genome.fa \
  --tsv /ruta/read_to_transcript.tsv.gz \
  --prefix sample \
  --output-dir output_sample \
  --modkit /ruta/modkit
```

Opciones útiles:

```text
--mods MOD [MOD ...]       Modificaciones a analizar
--min-genomic-cov INTEGER  Cobertura genómica mínima (5)
--prob-lim FLOAT           Umbral de probabilidad (0.95)
--min-tx-cov INTEGER       Cobertura mínima por transcrito (20)
--min-tx-occ FLOAT         Ocupancia mínima (0.2)
```

Comprueba siempre que el BAM tenga su índice `.bai` y que el nombre de cada `transcript_id` coincida entre el GTF, el archivo de IsoQuant y la clasificación de SQANTI3.

## Estructura del proyecto

```text
01_bedmethyl.sh                      Generación y filtrado de BED
02_mods_to_CIGAR_transcripts.sh      Proyección a coordenadas de transcrito
03_collapse_to_tx_and_merge.sh       Colapso y unión con SQANTI3
saca_tags_transcript_models.py       Procesamiento por lectura
collapse_to_tx_and_merge.py          Resumen por transcrito y merge
sqanti-mod.py                        Orquestador para una muestra
config.sh                            Configuración de ejecución
config.ini.example                   Plantilla de parámetros
```

## Notas

- Los scripts `.sh` están preparados para un entorno con módulos, Conda y SLURM; adapta las directivas `#SBATCH` a tu clúster.
- La extracción utiliza el número de CPUs definido en `SLURM_CPUS_PER_TASK`, con valor de respaldo `4` cuando se ejecuta fuera de SLURM.
- Las coordenadas transcriptómicas generadas por el pipeline son 0-based.
- Revisa las rutas y variables de cada muestra antes de lanzar trabajos en serie.