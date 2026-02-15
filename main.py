
import os
import sys
import time
import json
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any, Iterator, Tuple

from pyspark.sql import SparkSession

def env_int(name: str, default: int) -> int:
    v = os.environ.get(name, str(default)).strip()
    try:
        return int(v)
    except ValueError:
        raise ValueError(f"Env var {name} must be int, got {v!r}")


def env_str(name: str, default: str) -> str:
    v = os.environ.get(name, default)
    if v is None:
        return default
    return str(v).strip()


# Paths / inputs
ROOT_DIR        = Path(env_str("ROOT_DIR", "/data/genomics")).expanduser().resolve()
RUN_TABLE_PATH  = Path(env_str("RUN_TABLE_PATH", str(ROOT_DIR / "SraRunTable.tsv"))).expanduser().resolve()
SAMPLE_COL      = env_int("SAMPLE_COL", 29)  # 0-based; your shell uses $30 => python uses 29

REF_FA_GZ_URL   = env_str(
    "REF_FA_GZ_URL",
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Aegilops_tauschii/latest_assembly_versions/"
    "GCF_002575655.3_Aet_v6.0/GCF_002575655.3_Aet_v6.0_genomic.fna.gz",
)
REF_FA_NAME     = env_str("REF_FA_NAME", "GCF_002575655.3_Aet_v6.0_genomic.fna")

# Spark behavior
MAX_SAMPLES     = env_int("MAX_SAMPLES", 0)  # 0 => all
SPARK_PARTS     = env_int("SPARK_PARTITIONS", 1)

# Tool threads
THREADS_FASTQ   = env_int("THREADS_FASTQ", 4)
THREADS_FASTQC  = env_int("THREADS_FASTQC", 4)
THREADS_BWA     = env_int("THREADS_BWA", 4)
THREADS_SORT    = env_int("THREADS_SORT", THREADS_BWA)

COMPILE_SUMMARY = env_int("COMPILE_SUMMARY", 1)  # 1 => read per-sample JSON into Parquet+CSV at end
KEEP_WORKDIR    = env_int("KEEP_WORKDIR", 0)     # 1 => don't delete per-sample work dir
KEEP_FASTQ      = env_int("KEEP_FASTQ", 0)       # 1 => don't delete per-sample fastq after mapping

# SRA stability
SRA_CONFIG_DIR  = Path(env_str("SRA_CONFIG_DIR", str(ROOT_DIR / "ncbi"))).expanduser().resolve()

RUN_FASTQC      = env_int("RUN_FASTQC", 1)


def run(cmd, log: Path | None = None, env: dict | None = None, check: bool = True) -> int:
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("a", encoding="utf-8") as lf:
            lf.write(f"\n$ {' '.join(cmd)}\n")
            lf.flush()
            p = subprocess.Popen(cmd, stdout=lf, stderr=lf, env=env)
            rc = p.wait()
        if check and rc != 0:
            tail = ""
            try:
                lines = log.read_text(encoding="utf-8", errors="ignore").splitlines()
                tail = "\n".join(lines[-60:])
            except Exception:
                pass
            raise RuntimeError(f"Command failed ({rc}): {' '.join(cmd)}\n--- log tail ---\n{tail}")
        return rc

    p = subprocess.run(cmd, env=env)
    if check and p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}")
    return p.returncode


def ensure_sra_config():
    SRA_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    os.environ["VDB_CONFIG"] = str(SRA_CONFIG_DIR / "user-settings.mkfg")
    os.environ["VDB_CONFIG_PATH"] = str(SRA_CONFIG_DIR)

    cfg = Path(os.environ["VDB_CONFIG"])
    if not cfg.exists():
        try:
            subprocess.run(["vdb-config", "--restore-defaults"], check=False)
        except FileNotFoundError:
            # vdb-config may not be present; prefetch/fasterq-dump can still work.
            pass


def ensure_reference(root: Path) -> Path:
    ref_dir = root / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)

    ref_fa_gz = ref_dir / Path(REF_FA_GZ_URL).name
    ref_fa = ref_dir / REF_FA_NAME

    if not ref_fa.exists():
        if not ref_fa_gz.exists():
            run(["wget", "-nc", "-O", str(ref_fa_gz), REF_FA_GZ_URL], log=None)
        run(["gunzip", "-f", str(ref_fa_gz)], log=None)

    expected = [".bwt", ".sa", ".pac", ".ann", ".amb"]
    if not all(Path(str(ref_fa) + suf).exists() for suf in expected):
        run(["bwa", "index", str(ref_fa)], log=None)

    # samtools faidx
    if not Path(str(ref_fa) + ".fai").exists():
        run(["samtools", "faidx", str(ref_fa)], log=None)

    return ref_fa


def parse_flagstat_text(text: str) -> Dict[str, Any]:
    stats = {"total": 0, "mapped": 0, "properly_paired": 0, "singletons": 0, "duplicate_reads": 0}
    for line in text.splitlines():
        if " in total" in line:
            stats["total"] = int(line.split()[0])
        elif " mapped (" in line and "mate mapped" not in line:
            stats["mapped"] = int(line.split()[0])
        elif " properly paired (" in line:
            stats["properly_paired"] = int(line.split()[0])
        elif " singletons (" in line:
            stats["singletons"] = int(line.split()[0])
        elif " duplicates" in line:
            stats["duplicate_reads"] = int(line.split()[0])
    return stats


def parse_samtools_stats(stats_txt: str) -> Dict[str, Any]:
    out = {"average_length": None, "insert_size_mean": None, "insert_size_sd": None}
    for line in stats_txt.splitlines():
        if line.startswith("SN\taverage length:\t"):
            try:
                out["average_length"] = float(line.split("\t")[-1].strip())
            except ValueError:
                pass
        if line.startswith("SN\tinsert size average:\t"):
            try:
                out["insert_size_mean"] = float(line.split("\t")[-1].strip())
            except ValueError:
                pass
        if line.startswith("SN\tinsert size standard deviation:\t"):
            try:
                out["insert_size_sd"] = float(line.split("\t")[-1].strip())
            except ValueError:
                pass
    return out


def parse_idxstats(idxstats_txt: str) -> Dict[str, Any]:
    total_mapped = 0
    ref_with_mapped = 0
    for line in idxstats_txt.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            mapped = int(parts[2])
        except ValueError:
            continue
        total_mapped += mapped
        if mapped > 0:
            ref_with_mapped += 1
    return {"idxstats_total_mapped": total_mapped, "idxstats_ref_with_mapped": ref_with_mapped}


def parse_fastqc_summary(zip_path: Path) -> Dict[str, int]:
    import zipfile
    counts = {"PASS": 0, "WARN": 0, "FAIL": 0}
    with zipfile.ZipFile(zip_path) as z:
        with z.open("summary.txt") as f:
            for line in f:
                status = line.decode().split("\t")[0]
                if status in counts:
                    counts[status] += 1
    return counts


def write_json_atomic(out_path: Path, payload: Dict[str, Any]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + f".tmp.{os.getpid()}")
    tmp_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.replace(out_path)



def process_sample(row: Tuple[str, str], root: Path, ref_fa: Path) -> Optional[Dict[str, Any]]:
    srr, sample = row
    ensure_sra_config()

    logs = root / "logs"
    cram = root / "cram"
    fastqc_out = root / "fastqc"
    metrics_dir = root / "metrics" / "samples"

    # per-sample sandbox
    work = root / "work" / srr
    sra_dir = work / "sra"
    fastq_dir = work / "fastq"
    tmp_dir = work / "tmp"

    for d in [logs, cram, fastqc_out, metrics_dir, work, sra_dir, fastq_dir, tmp_dir]:
        d.mkdir(parents=True, exist_ok=True)

    log = logs / f"{srr}.log"
    out_json = metrics_dir / f"{srr}.json"

    # Resume: if JSON already exists, skip
    if out_json.exists():
        return None

    cram_file = cram / f"{srr}.cram"
    crai_file = cram / f"{srr}.cram.crai"
    have_alignment = cram_file.exists() and crai_file.exists()

    r1 = fastq_dir / f"{srr}_1.fastq"
    r2 = fastq_dir / f"{srr}_2.fastq"

    if not have_alignment:
        # 1) Download SRA
        if not (sra_dir / srr).exists():
            run(["prefetch", srr, "-O", str(sra_dir)], log)

        # 2) FASTQ extraction
        if not (r1.exists() and r2.exists()):
            run(
                ["fasterq-dump", str(sra_dir / srr),
                 "-O", str(fastq_dir),
                 "--split-files",
                 "--threads", str(THREADS_FASTQ)],
                log
            )

        # 3) FastQC (optional)
        fastqc_zip = fastqc_out / f"{srr}_1_fastqc.zip"
        if RUN_FASTQC and not fastqc_zip.exists():
            run(["fastqc", "-o", str(fastqc_out), "-t", str(THREADS_FASTQC), str(r1), str(r2)], log)

        # 4) Mapping → CRAM (pipefail so BWA errors are not hidden)
        rg = f"@RG\\tID:{srr}\\tSM:{sample}\\tPL:ILLUMINA"
        tmp_prefix = tmp_dir / f"{srr}_{int(time.time())}"
        cmd = (
            "set -euo pipefail; "
            f"bwa mem -t {THREADS_BWA} -R '{rg}' '{ref_fa}' '{r1}' '{r2}' | "
            f"samtools sort -@ {THREADS_SORT} -O cram -T '{tmp_prefix}' "
            f"-o '{cram_file}' --reference '{ref_fa}'"
        )
        run(["bash", "-lc", cmd], log)

        # 5) Index CRAM
        run(["samtools", "index", str(cram_file), "--reference", str(ref_fa)], log)

    # -------------------------
    # Metrics (write immediately)
    # -------------------------
    flagstat_txt = subprocess.check_output(
        ["samtools", "flagstat", str(cram_file), "--reference", str(ref_fa)],
        text=True
    )
    (logs / f"{srr}.flagstat.txt").write_text(flagstat_txt, encoding="utf-8")
    fs = parse_flagstat_text(flagstat_txt)

    stats_txt = subprocess.check_output(
        ["samtools", "stats", str(cram_file), "--reference", str(ref_fa)],
        text=True
    )
    (logs / f"{srr}.samtools_stats.txt").write_text(stats_txt, encoding="utf-8")
    st = parse_samtools_stats(stats_txt)

    idxstats_txt = subprocess.check_output(["samtools", "idxstats", str(cram_file)], text=True)
    (logs / f"{srr}.idxstats.txt").write_text(idxstats_txt, encoding="utf-8")
    ix = parse_idxstats(idxstats_txt)

    # FastQC summary if present
    qc = {"PASS": 0, "WARN": 0, "FAIL": 0}
    fastqc_zip = fastqc_out / f"{srr}_1_fastqc.zip"
    if RUN_FASTQC and fastqc_zip.exists():
        qc = parse_fastqc_summary(fastqc_zip)

    total = fs["total"]
    mapped = fs["mapped"]
    properly_paired = fs["properly_paired"]

    payload = {
        "sample_id": srr,
        "sample_name": sample,

        "total_reads": total,
        "mapped_reads": mapped,
        "mapping_rate": (mapped / total) if total else 0.0,

        "properly_paired_reads": properly_paired,
        "properly_paired_rate": (properly_paired / total) if total else 0.0,

        "duplicate_reads": fs["duplicate_reads"],
        "duplicate_rate": (fs["duplicate_reads"] / total) if total else 0.0,

        "singletons": fs["singletons"],
        "singleton_rate": (fs["singletons"] / total) if total else 0.0,

        "average_read_length": st["average_length"],
        "insert_size_mean": st["insert_size_mean"],
        "insert_size_sd": st["insert_size_sd"],

        "idxstats_total_mapped": ix["idxstats_total_mapped"],
        "idxstats_ref_with_mapped": ix["idxstats_ref_with_mapped"],

        "fastqc_pass": qc["PASS"],
        "fastqc_warn": qc["WARN"],
        "fastqc_fail": qc["FAIL"],

        "cram_path": str(cram_file),
        "timestamp": time.time(),
    }

    write_json_atomic(out_json, payload)

    # Cleanup
    if not KEEP_WORKDIR:
        shutil.rmtree(work, ignore_errors=True)
    if not KEEP_FASTQ:
        try:
            r1.unlink(missing_ok=True)
            r2.unlink(missing_ok=True)
        except Exception:
            pass

    return payload


def map_partitions(it: Iterator[Tuple[str, str]], root: Path, ref_fa: Path) -> Iterator[int]:
    """
    Sequential within partition; concurrency is controlled by number of partitions.
    """
    for row in it:
        srr = row[0] if row else "UNKNOWN"
        try:
            out = process_sample(row, root, ref_fa)
            if out is not None:
                yield 1
        except Exception as e:
            (root / "logs" / f"{srr}.FAILED").write_text(str(e), encoding="utf-8")
            continue


def main() -> int:
    ROOT_DIR.mkdir(parents=True, exist_ok=True)

    spark = SparkSession.builder.appName(env_str("SPARK_APP_NAME", "PlantGenomicsSparkPipeline")).getOrCreate()
    sc = spark.sparkContext

    ref_fa = ensure_reference(ROOT_DIR)

    print("ROOT_DIR:", ROOT_DIR)
    print("RUN_TABLE_PATH:", RUN_TABLE_PATH)
    print("REF_FA:", ref_fa)
    print("SAMPLE_COL (0-based):", SAMPLE_COL)
    print("MAX_SAMPLES:", MAX_SAMPLES, "(0 => all)")
    print("SPARK_PARTITIONS:", SPARK_PARTS)

    if not RUN_TABLE_PATH.exists():
        print(f"ERROR: Missing run table at: {RUN_TABLE_PATH}")
        return 1

    lines = sc.textFile(str(RUN_TABLE_PATH))

    def parse(line: str):
        parts = line.split("\t")
        if not parts or parts[0] == "Run":
            return None
        srr = parts[0].strip()
        if not srr:
            return None
        sample = parts[SAMPLE_COL].strip() if len(parts) > SAMPLE_COL and parts[SAMPLE_COL].strip() else "UNKNOWN"
        return (srr, sample)

    pairs = lines.map(parse).filter(lambda x: x is not None)

    if MAX_SAMPLES > 0:
        pairs = sc.parallelize(pairs.take(MAX_SAMPLES))

    pairs = pairs.repartition(SPARK_PARTS)

    # Trigger processing without collecting everything
    completed = pairs.mapPartitions(lambda it: map_partitions(it, ROOT_DIR, ref_fa)).sum()
    print("Samples completed (newly written):", completed)

    if COMPILE_SUMMARY:
        sample_json_dir = ROOT_DIR / "metrics" / "samples"
        if sample_json_dir.exists():
            df = spark.read.json(str(sample_json_dir / "*.json"))
            df.write.mode("overwrite").parquet(str(ROOT_DIR / "metrics" / "sample_metrics.parquet"))
            df.coalesce(1).write.mode("overwrite").option("header", "true").csv(
                str(ROOT_DIR / "metrics" / "sample_metrics_csv")
            )
            df.show(truncate=False)

    spark.stop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
