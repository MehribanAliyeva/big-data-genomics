# Scalable Quality-Aware Mapping and Gene-Level Anomaly Detection in Plant Genome Sequencing Using Apache Spark

This repository implements a **Spark-based genomics pipeline** for scalable processing of plant whole-genome sequencing (WGS) data.  
The pipeline performs **quality-aware mapping**, extracts **sample-level metrics**, and prepares structured data suitable for **anomaly detection**.

The workflow mirrors a traditional bash-based bioinformatics pipeline but adds:
- Parallel execution using Apache Spark
- Structured metric extraction (JSON / Parquet / CSV)
- Resume-safe execution
- Scalability to larger machines

This README assumes that **`SraRunTable.tsv` already exists** in the genomics root directory.

---

## System Requirements

- Linux (Ubuntu 20.04+ recommended)
- Java 11
- Python 3.8+
- Apache Spark (standalone)
- Disk space: **hundreds of GB recommended for WGS**

---

## Install System Dependencies

```bash
# 1) System deps
sudo apt-get update
sudo apt-get install -y \
  curl wget git unzip tar gzip \
  build-essential \
  python3 python3-pip \
  openjdk-11-jre-headless \
  bwa samtools fastqc

mkdir -p ~/tools && cd ~/tools
wget -O sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.tar.gz
SRA_DIR=$(find ~/tools -maxdepth 1 -type d -name "sratoolkit.*-ubuntu64" | head -n 1)
echo "export PATH=\$PATH:$SRA_DIR/bin" >> ~/.bashrc

cd ~/tools
wget -O spark.tgz https://archive.apache.org/dist/spark/spark-3.4.2/spark-3.4.2-bin-hadoop3.tgz
tar -xzf spark.tgz
mv spark-3.4.2-bin-hadoop3 spark
echo 'export SPARK_HOME=$HOME/tools/spark' >> ~/.bashrc
echo 'export PATH=$PATH:$SPARK_HOME/bin:$SPARK_HOME/sbin' >> ~/.bashrc
echo 'export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64' >> ~/.bashrc

source ~/.bashrc
spark-submit --version
prefetch --version || true
curl -Ls https://astral.sh/uv/install.sh | bash
source ~/.bashrc
uv --version
cd /path/to/this/repo

uv sync
uv sync
bash run_spark.sh

```


