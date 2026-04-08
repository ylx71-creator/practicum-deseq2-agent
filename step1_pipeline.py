import re
from pathlib import Path
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def run_pipeline(
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
) -> dict:

    base = Path(__file__).resolve().parent
    data_dir, out_dir = base / "data", base / "output"
    out_dir.mkdir(exist_ok=True)

    # --- metadata: sample_id + risk category ---
    lines = (data_dir / "GSE164641_series_matrix.txt").read_text(
        encoding="utf-8", errors="ignore"
    ).splitlines()

    gsm = re.findall(r"GSM\d+", next(l for l in lines if l.startswith("!Sample_geo_accession")))
    cat_line = next(l for l in lines if "risk category" in l.lower())
    cats = re.findall(r"risk category: ([A-Za-z]+)", cat_line)

    n = min(len(gsm), len(cats))
    clin = pd.DataFrame({"sample_id": gsm[:n], "target": cats[:n]})

    # --- counts ---
    cnt = pd.read_csv(
        data_dir / "GSE164641_raw_counts_GRCh38.p13_NCBI.tsv",
        sep="\t"
    ).set_index("GeneID").T
    cnt.index.name = "sample_id"
    cnt = cnt.reset_index()
    cnt.columns = cnt.columns.map(str)

    # --- merge ---
    master = clin.merge(cnt, on="sample_id", how="inner").set_index("sample_id")

    # --- annotation ---
    gene_map = {}
    annot_path = data_dir / "Human.GRCh38.p13.annot.tsv"
    if annot_path.exists():
        annot = pd.read_csv(annot_path, sep="\t", dtype=str)
        if "GeneID" in annot.columns and "Symbol" in annot.columns:
            annot["GeneID"] = annot["GeneID"].str.split(".").str[0]
            gene_map = dict(zip(annot["GeneID"], annot["Symbol"]))

    # --- DESeq2 input ---
    condition = master["target"].astype(str)
    X = master.drop(columns=["target"]).apply(pd.to_numeric, errors="raise").round().astype(int)
    meta = pd.DataFrame({"condition": condition.values}, index=X.index)

    print("Samples in metadata:", len(clin))
    print("Samples after merge:", master.shape[0])
    print("Group counts:\n", condition.value_counts())

    # --- DESeq2 ---
    dds = DeseqDataSet(counts=X, metadata=meta, design_factors="condition")
    dds.deseq2()

    norm = pd.DataFrame(dds.layers["normed_counts"], index=X.index, columns=X.columns)
    norm.to_csv(out_dir / "deseq2_normalized_counts.csv")

    stat = DeseqStats(dds, contrast=["condition", "High", "Average"])
    stat.summary()
    res = stat.results_df.copy()
    res["GeneID"] = res.index.astype(str).str.split(".").str[0]
    res.insert(0, "gene_symbol", res["GeneID"].map(lambda gid: gene_map.get(gid, "NA")))
    res.to_csv(out_dir / "deseq2_statistical_results.csv", index=False)

    sig = res.dropna(subset=["padj", "log2FoldChange"])
    sig = sig[(sig["padj"] < padj_cutoff) & (sig["log2FoldChange"].abs() > lfc_cutoff)].copy()
    sig.to_csv(out_dir / "deseq2_significant_genes.csv", index=False)

    return {
        "n_used": int(master.shape[0]),
        "n_genes": int(X.shape[1]),
        "n_sig": int(sig.shape[0]),
        "outputs": {
            "normed": str(out_dir / "deseq2_normalized_counts.csv"),
            "de": str(out_dir / "deseq2_statistical_results.csv"),
            "sig": str(out_dir / "deseq2_significant_genes.csv"),
        },
    }

if __name__ == "__main__":
    run_pipeline()