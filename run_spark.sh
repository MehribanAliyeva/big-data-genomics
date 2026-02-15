#!/usr/bin/env bash
set -euo pipefail

export ROOT_DIR="${ROOT_DIR:-/data/genomics}"
export RUN_TABLE_PATH="${RUN_TABLE_PATH:-$ROOT_DIR/SraRunTable.tsv}"

export REF_FA_GZ_URL="${REF_FA_GZ_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Aegilops_tauschii/latest_assembly_versions/GCF_002575655.3_Aet_v6.0/GCF_002575655.3_Aet_v6.0_genomic.fna.gz}"
export REF_FA_NAME="${REF_FA_NAME:-GCF_002575655.3_Aet_v6.0_genomic.fna}"

export SAMPLE_COL="${SAMPLE_COL:-29}"
export MAX_SAMPLES="${MAX_SAMPLES:-0}"

export SPARK_PARTITIONS="${SPARK_PARTITIONS:-4}"
export THREADS_FASTQ="${THREADS_FASTQ:-4}"
export THREADS_FASTQC="${THREADS_FASTQC:-4}"
export THREADS_BWA="${THREADS_BWA:-8}"
export THREADS_SORT="${THREADS_SORT:-8}"

export COMPILE_SUMMARY="${COMPILE_SUMMARY:-1}"
export KEEP_WORKDIR="${KEEP_WORKDIR:-0}"
export KEEP_FASTQ="${KEEP_FASTQ:-0}"
export RUN_FASTQC="${RUN_FASTQC:-1}"
export SPARK_DRIVER_MEMORY="${SPARK_DRIVER_MEMORY:-8g}"
export SPARK_EXECUTOR_MEMORY="${SPARK_EXECUTOR_MEMORY:-8g}"
export SPARK_EXECUTOR_CORES="${SPARK_EXECUTOR_CORES:-4}"
export SPARK_EXECUTOR_INSTANCES="${SPARK_EXECUTOR_INSTANCES:-1}"

export SRA_CONFIG_DIR="${SRA_CONFIG_DIR:-$ROOT_DIR/ncbi}"


SPARK_MASTER="${SPARK_MASTER:-local[*]}"
SPARK_APP_NAME="${SPARK_APP_NAME:-PlantGenomicsSparkPipeline}"

SPARK_EXTRA_CONF="${SPARK_EXTRA_CONF:-}"

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
MAIN_PY="${MAIN_PY:-$REPO_DIR/main.py}"

echo "Repo:            $REPO_DIR"
echo "Main:            $MAIN_PY"
echo "Spark master:    $SPARK_MASTER"
echo "ROOT_DIR:        $ROOT_DIR"
echo "RUN_TABLE_PATH:  $RUN_TABLE_PATH"
echo "MAX_SAMPLES:     $MAX_SAMPLES"
echo "SPARK_PARTS:     $SPARK_PARTITIONS"

spark-submit \
  --name "$SPARK_APP_NAME" \
  --master "$SPARK_MASTER" \
  --driver-memory "$SPARK_DRIVER_MEMORY" \
  --executor-memory "$SPARK_EXECUTOR_MEMORY" \
  --executor-cores "$SPARK_EXECUTOR_CORES" \
  --num-executors "$SPARK_EXECUTOR_INSTANCES" \
  $SPARK_EXTRA_CONF \
  "$MAIN_PY"

