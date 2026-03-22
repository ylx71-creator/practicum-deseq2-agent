import re
from pathlib import Path

import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# =========================
# Paths
# =========================
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "output"
OUT_DIR.mkdir(exist_ok=True)

SERIES_PATH = DATA_DIR / "GSE164641_series_matrix.txt"
COUNTS_PATH = DATA_DIR / "GSE164641_raw_counts_GRCh38.p13_NCBI.tsv"


MASTER_OUT = OUT_DIR / "GSE164641_master_dataframe.csv"
NORMED_OUT = OUT_DIR / "deseq2_normalized_counts.csv"
DE_OUT = OUT_DIR / "deseq2_statistical_results.csv"


# =========================
# Part 1 helpers
# =========================
def _read_lines(fp: Path) -> list[str]:
    with open(fp, "r", encoding="utf-8", errors="ignore") as f:
        return f.readlines()


def _extract_gsm_list(lines: list[str]) -> list[str]:
    gsm_line = next((ln for ln in lines if ln.startswith("!Sample_geo_accession")), None)
    if gsm_line is None:
        raise ValueError("Cannot find '!Sample_geo_accession' in series matrix.")
    return re.findall(r"(GSM\d+)", gsm_line)


def _extract_tyrer_cuzick_risk_percent(lines: list[str]) -> list[float]:
    # Find the line containing: life time risk (tyrer-cuzick score): XX.XX%
    risk_line = None
    for ln in lines:
        low = ln.lower()
        if "life time risk" in low and "tyrer-cuzick" in low:
            risk_line = ln
            break
    if risk_line is None:
        raise ValueError("Cannot find 'life time risk (tyrer-cuzick score)' line in series matrix.")

    vals = re.findall(r"(\d+\.\d+)%", risk_line)
    return [float(x) for x in vals]


def _build_clinical_df(series_path: Path, threshold: float = 20.0) -> pd.DataFrame:
    """
    Build a minimal clinical dataframe:
      sample_id, tyrer_cuzick_score (float, percent), target (High/Average)
    Robust to mismatch lengths by truncating to min length and later inner-joining with counts.
    """
    lines = _read_lines(series_path)
    gsm_list = _extract_gsm_list(lines)
    risk_vals = _extract_tyrer_cuzick_risk_percent(lines)

    n = min(len(gsm_list), len(risk_vals))
    gsm_list = gsm_list[:n]
    risk_vals = risk_vals[:n]

    clinical_df = pd.DataFrame(
        {
            "sample_id": gsm_list,
            "tyrer_cuzick_score": risk_vals,
        }
    )
    clinical_df["target"] = np.where(clinical_df["tyrer_cuzick_score"] >= threshold, "High", "Average")

    print(f"[Part1] GSM in series_matrix: {len(_extract_gsm_list(lines))}")
    print(f"[Part1] risk values parsed: {len(_extract_tyrer_cuzick_risk_percent(lines))}")
    print(f"[Part1] using paired rows: {n}")
    print(f"[Part1] High={int((clinical_df['target']=='High').sum())}, "
          f"Average={int((clinical_df['target']=='Average').sum())}")

    return clinical_df


def _load_counts_df(counts_path: Path) -> pd.DataFrame:
    """
    Load counts tsv:
      rows: GeneID
      cols: GSM...
    Convert to:
      rows: sample_id
      cols: genes
    """
    counts = pd.read_csv(counts_path, sep="\t")
    if "GeneID" not in counts.columns:
        raise ValueError("Counts file must contain 'GeneID' column.")
    counts = counts.set_index("GeneID").T
    counts.index.name = "sample_id"
    counts.reset_index(inplace=True)

    # Ensure gene columns are strings
    counts.columns = counts.columns.map(str)

    print(f"[Part1] counts loaded: samples={counts.shape[0]}, genes={counts.shape[1]-1}")
    return counts


def _load_gene_id_map(annot_path: Path) -> dict[str, str]:
    """
    Optional mapping gene_id -> gene_symbol from annotation TSV.
    """
    if not annot_path.exists():
        print(f"[Part1] annotation file not found, skip mapping: {annot_path}")
        return {}

    print(f"[Part1] loading annotations: {annot_path}")
    annot_df = pd.read_csv(
        annot_path,
        sep="\t",
        header=None,
        usecols=[0, 1],
        on_bad_lines="skip"
    )
    annot_df.columns = ["gene_id", "gene_symbol"]
    annot_df.dropna(subset=["gene_id", "gene_symbol"], inplace=True)

    m = dict(zip(annot_df["gene_id"].astype(str), annot_df["gene_symbol"].astype(str)))
    print(f"[Part1] mapping size: {len(m)}")
    return m


def _make_master_df(series_path: Path, counts_path: Path, annot_path: Path | None = None) -> pd.DataFrame:
    clinical_df = _build_clinical_df(series_path, threshold=20.0)
    counts_df = _load_counts_df(counts_path)

    # Optional: map gene_id -> symbol
    gene_id_map = _load_gene_id_map(annot_path) if annot_path is not None else {}
    if gene_id_map:
        # Rename gene columns only (exclude 'sample_id')
        gene_cols = counts_df.columns[1:]
        renamed = {c: gene_id_map.get(c, c) for c in gene_cols}
        counts_df = counts_df.rename(columns=renamed)
        print("[Part1] gene IDs mapped to symbols (note: symbols may not be unique).")

    master_df = pd.merge(clinical_df, counts_df, on="sample_id", how="inner")
    master_df = master_df.set_index("sample_id")

    print(f"[Part1] master_df: samples={master_df.shape[0]}, columns={master_df.shape[1]}")
    return master_df


# =========================
# Part 3: DESeq2
# =========================
def run_deseq2_from_master(master_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (normalized_counts_df, results_df)
    """
    if "target" not in master_df.columns:
        raise ValueError("master_df must contain 'target' column (High/Average).")

    # y: High=1, Average=0
    y = master_df["target"].map({"High": 1, "Average": 0})
    if y.isna().any():
        raise ValueError("Found NA in y mapping. Check 'target' values.")

    # Identify non-gene columns to drop (keep only gene expression counts)
    non_gene_cols = {"tyrer_cuzick_score", "target"}
    X_cont_raw = master_df.drop(columns=[c for c in non_gene_cols if c in master_df.columns])

    # Ensure numeric & integer
    X_cont_raw = X_cont_raw.apply(pd.to_numeric, errors="coerce")
    if X_cont_raw.isna().any().any():
        # If mapping to symbols created duplicate columns, some ops can behave oddly; still, we should fail fast.
        raise ValueError("Non-numeric values detected in counts matrix after conversion.")

    X_cont_raw = X_cont_raw.round().astype(int)

    # metadata aligned to counts
    samples_df = pd.DataFrame({"condition": y.astype(int).values}, index=X_cont_raw.index)

    print("\n--- Running DESeq2 Differential Expression Analysis ---")
    print(f"[Part3] counts matrix: samples={X_cont_raw.shape[0]}, genes={X_cont_raw.shape[1]}")
    print(f"[Part3] High={int(samples_df['condition'].sum())}, "
          f"Average={int((samples_df['condition']==0).sum())}")

    dds = DeseqDataSet(
        counts=X_cont_raw,
        metadata=samples_df,
        design_factors="condition"
    )
    dds.deseq2()

    normalized_counts_df = pd.DataFrame(
        dds.layers["normed_counts"],
        index=X_cont_raw.index,
        columns=X_cont_raw.columns
    )

    stat_res = DeseqStats(dds, contrast=["condition", 1, 0])
    stat_res.summary()
    results_df = stat_res.results_df

    return normalized_counts_df, results_df


def main():
    # -------- Part 1 --------
    master_df = _make_master_df(SERIES_PATH, COUNTS_PATH)
    master_df.to_csv(MASTER_OUT)
    print(f"[Part1] saved master_df -> {MASTER_OUT}")

    # -------- Part 3 --------
    normalized_counts_df, results_df = run_deseq2_from_master(master_df)
    normalized_counts_df.to_csv(NORMED_OUT)
    results_df.to_csv(DE_OUT)

    print(f"\n[Part3] saved normalized counts -> {NORMED_OUT}")
    print(f"[Part3] saved DE results         -> {DE_OUT}")

    sig = results_df[(results_df["padj"] < 0.05) & (results_df["log2FoldChange"].abs() > 1)].copy()
    sig.to_csv(OUT_DIR / "deseq2_significant_genes.csv")
    print("saved significant genes ->", OUT_DIR / "deseq2_significant_genes.csv")

    # Quick sanity summary
    if "padj" in results_df.columns:
        print("\n[Check] min pvalue:", float(results_df["pvalue"].min()))
        print("[Check] min padj  :", float(results_df["padj"].min()))
        sig = results_df.dropna(subset=["padj", "log2FoldChange"])
        n_sig = int(((sig["padj"] < 0.05) & (sig["log2FoldChange"].abs() > 1.0)).sum())
        print("[Check] significant genes (padj<0.05 & |log2FC|>1):", n_sig)


if __name__ == "__main__":
    main()

    