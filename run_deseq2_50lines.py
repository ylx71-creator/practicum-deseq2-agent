import re
from pathlib import Path
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# --- Paths ---
BASE = Path(__file__).resolve().parent
DATA, OUT = BASE / "data", BASE / "output"
OUT.mkdir(exist_ok=True)
SERIES = DATA / "GSE164641_series_matrix.txt"
COUNTS = DATA / "GSE164641_raw_counts_GRCh38.p13_NCBI.tsv"
ANNOT = DATA / "Human.GRCh38.p13.annot.tsv"

def load_map(annot_fp: Path) -> dict:
    """gene_id (col0) -> gene_symbol (col1); tab-separated."""
    if not annot_fp.exists(): return {}
    a = pd.read_csv(annot_fp, sep="\t", header=None, usecols=[0, 1], on_bad_lines="skip")
    a.columns = ["gene_id", "gene_symbol"]
    a = a.dropna()
    return dict(zip(a["gene_id"].astype(str), a["gene_symbol"].astype(str)))

def build_master():
    # 1) parse GSM + Tyrer-Cuzick risk%
    lines = SERIES.read_text(encoding="utf-8", errors="ignore").splitlines()
    gsm = re.findall(r"(GSM\d+)", next(l for l in lines if l.startswith("!Sample_geo_accession")))
    risk_line = next(l for l in lines if "life time risk" in l.lower() and "tyrer-cuzick" in l.lower())
    risk = [float(x) for x in re.findall(r"(\d+\.\d+)%", risk_line)]
    n = min(len(gsm), len(risk))
    clin = pd.DataFrame({"sample_id": gsm[:n], "tyrer_cuzick_score": risk[:n]})
    clin["target"] = np.where(clin["tyrer_cuzick_score"] >= 20.0, "High", "Average")

    # 2) load counts: GeneID x GSM... -> sample_id x genes
    cnt = pd.read_csv(COUNTS, sep="\t").set_index("GeneID").T
    cnt.index.name = "sample_id"
    cnt = cnt.reset_index()
    cnt.columns = cnt.columns.map(str)

    master = clin.merge(cnt, on="sample_id", how="inner").set_index("sample_id")
    print(f"[Part1] GSM={len(gsm)}, risk={len(risk)}, used={n}, merged={master.shape[0]} samples, genes={cnt.shape[1]-1}")
    return master

def main():
    gene_map = load_map(ANNOT)
    master = build_master()
    master.to_csv(OUT / "GSE164641_master_dataframe.csv")

    # DESeq2: High(1) vs Average(0)
    y = master["target"].map({"High": 1, "Average": 0}).astype(int)
    X = master.drop(columns=["tyrer_cuzick_score", "target"]).apply(pd.to_numeric, errors="raise").round().astype(int)
    meta = pd.DataFrame({"condition": y.values}, index=X.index)

    dds = DeseqDataSet(counts=X, metadata=meta, design_factors="condition")
    dds.deseq2()
    pd.DataFrame(dds.layers["normed_counts"], index=X.index, columns=X.columns).to_csv(OUT / "deseq2_normalized_counts.csv")

    stat = DeseqStats(dds, contrast=["condition", 1, 0])
    stat.summary()              
    res = stat.results_df
    if gene_map: res["gene_symbol"] = res.index.astype(str).map(gene_map)
    res.to_csv(OUT / "deseq2_statistical_results.csv")

    sig = res.dropna(subset=["padj", "log2FoldChange"])
    sig = sig[(sig["padj"] < 0.05) & (sig["log2FoldChange"].abs() > 1)]
    sig.to_csv(OUT / "deseq2_significant_genes.csv")
    print(f"[Done] sig genes={sig.shape[0]} | symbol_mapped={'yes' if gene_map else 'no'}")

if __name__ == "__main__":
    main()