import re
from pathlib import Path
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def run_pipeline(
    base_dir: Path,
    threshold: float = 20.0,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
) -> dict:
   
    data, out = base_dir/"data", base_dir/"output"
    out.mkdir(exist_ok=True)
    series = (data/"GSE164641_series_matrix.txt").read_text(encoding="utf-8", errors="ignore").splitlines()

    gsm = re.findall(r"(GSM\d+)", next(l for l in series if l.startswith("!Sample_geo_accession")))
    risk_line = next(l for l in series if "life time risk" in l.lower() and "tyrer-cuzick" in l.lower())
    risk = [float(x) for x in re.findall(r"(\d+\.\d+)%", risk_line)]
    n = min(len(gsm), len(risk))
    clin = pd.DataFrame({"sample_id": gsm[:n], "score": risk[:n]})
    clin["target"] = np.where(clin["score"] >= threshold, "High", "Average")

    cnt = pd.read_csv(data/"GSE164641_raw_counts_GRCh38.p13_NCBI.tsv", sep="\t").set_index("GeneID").T
    cnt.index.name = "sample_id"; cnt = cnt.reset_index(); cnt.columns = cnt.columns.map(str)
    master = clin.merge(cnt, on="sample_id", how="inner").set_index("sample_id")

    gene_map = {}
    ap = data/"Human.GRCh38.p13.annot.tsv"
    if ap.exists():
        a = pd.read_csv(ap, sep="\t", dtype=str)
        if "GeneID" in a.columns and "Symbol" in a.columns:
            gene_map = dict(zip(a["GeneID"], a["Symbol"]))

    y = master["target"].map({"High": 1, "Average": 0}).astype(int)
    X = master.drop(columns=["score", "target"]).apply(pd.to_numeric, errors="raise").round().astype(int)
    dds = DeseqDataSet(counts=X, metadata=pd.DataFrame({"condition": y.values}, index=X.index), design_factors="condition")
    dds.deseq2()

    norm = pd.DataFrame(dds.layers["normed_counts"], index=X.index, columns=X.columns)
    norm.to_csv(out/"deseq2_normalized_counts.csv")

    stat = DeseqStats(dds, contrast=["condition", 1, 0]); stat.summary()
    res = stat.results_df.copy()
    res.to_csv(out/"deseq2_statistical_results.csv")


    sig = res.dropna(subset=["padj","log2FoldChange"])
    sig = sig[(sig["padj"] < padj_cutoff) & (sig["log2FoldChange"].abs() > lfc_cutoff)].copy()
    sig.insert(0, "gene_symbol", sig.index.map(lambda gid: gene_map.get(str(gid), "NA")))
    sig.to_csv(out/"deseq2_significant_genes.csv")
    
    return {"n_used": int(master.shape[0]), "n_genes": int(X.shape[1]), "n_sig": int(sig.shape[0]),
            "outputs": {"normed": str(out/"deseq2_normalized_counts.csv"),
                        "de": str(out/"deseq2_statistical_results.csv"),
                        "sig": str(out/"deseq2_significant_genes.csv")}}

if __name__ == "__main__":
    s = run_pipeline(Path(__file__).resolve().parent)
    print("[Done]", s)

    