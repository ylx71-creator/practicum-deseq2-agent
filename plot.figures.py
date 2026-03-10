import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# =========================
# Paths
# =========================
BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "output"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

DE_PATH = OUT_DIR / "deseq2_statistical_results.csv"
NORM_PATH = OUT_DIR / "deseq2_normalized_counts.csv"
MASTER_PATH = OUT_DIR / "GSE164641_master_dataframe.csv"


# =========================
# Utilities
# =========================
def _savefig(name: str, dpi: int = 300):
    fp = FIG_DIR / name
    plt.savefig(fp, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"[saved] {fp}")


def _load_master_labels(master_path: Path) -> pd.Series:
    """
    Return y (0/1) aligned to sample_id index from master_df:
      High -> 1, Average -> 0
    """
    master = pd.read_csv(master_path, index_col=0)
    if "target" not in master.columns:
        raise ValueError("GSE164641_master_dataframe.csv must contain 'target' column.")
    y = master["target"].map({"High": 1, "Average": 0})
    if y.isna().any():
        raise ValueError("Found NA when mapping target to 0/1. Check 'target' values.")
    return y.astype(int)


# =========================
# 1) Volcano plot
# =========================
def volcano_plot(de_path: Path):
    de = pd.read_csv(de_path, index_col=0)

    need_cols = {"log2FoldChange", "pvalue", "padj"}
    if not need_cols.issubset(set(de.columns)):
        raise ValueError(f"DE results missing columns. Need: {need_cols}")

    df = de[["log2FoldChange", "pvalue", "padj"]].copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    # -log10(pvalue)
    df["neglog10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))

    # significance rule consistent with your code
    sig = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1.0)

    plt.figure(figsize=(9, 7))
    plt.scatter(df.loc[~sig, "log2FoldChange"], df.loc[~sig, "neglog10p"], s=8, alpha=0.5)
    plt.scatter(df.loc[sig, "log2FoldChange"], df.loc[sig, "neglog10p"], s=10, alpha=0.9)

    # threshold lines
    plt.axvline(1.0, linestyle="--", linewidth=1)
    plt.axvline(-1.0, linestyle="--", linewidth=1)
    plt.axhline(-np.log10(0.05), linestyle="--", linewidth=1)

    plt.title("Volcano Plot (High vs Average)")
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(p-value)")
    plt.tight_layout()
    _savefig("volcano_plot.png")


# =========================
# 2) PCA plot (all genes or subset)
# =========================
def pca_plot(norm_counts: pd.DataFrame, y: pd.Series, out_name: str, title: str):
    """
    norm_counts: rows=samples, cols=genes (numeric)
    y: 0/1 indexed by sample_id
    """
    # align
    common = norm_counts.index.intersection(y.index)
    X = norm_counts.loc[common]
    yy = y.loc[common]

    # log1p for stability
    X_log = np.log1p(X)

    # scale then PCA
    X_scaled = StandardScaler().fit_transform(X_log.values)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X_scaled)

    var = pca.explained_variance_ratio_
    pc_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=common)
    pc_df["group"] = yy.map({1: "High", 0: "Average"}).values

    plt.figure(figsize=(9, 7))
    for g in ["High", "Average"]:
        sub = pc_df[pc_df["group"] == g]
        plt.scatter(sub["PC1"], sub["PC2"], s=25, alpha=0.8, label=g)

    plt.title(title)
    plt.xlabel(f"PC1 ({var[0]*100:.1f}% var)")
    plt.ylabel(f"PC2 ({var[1]*100:.1f}% var)")
    plt.legend()
    plt.tight_layout()
    _savefig(out_name)


# =========================
# 3) Heatmap (DEG subset)
# =========================
def heatmap_deg(norm_counts: pd.DataFrame, y: pd.Series, deg_genes: list[str]):
    """
    Make a clustered heatmap (simple + dependency-light):
    - subset to deg_genes
    - z-score per gene
    - order samples by group
    """
    common = norm_counts.index.intersection(y.index)
    X = norm_counts.loc[common]
    yy = y.loc[common].map({1: "High", 0: "Average"})

    # keep only genes that exist
    genes = [g for g in deg_genes if g in X.columns]
    if len(genes) == 0:
        raise ValueError("No DEG genes found in normalized counts columns.")

    Xg = np.log1p(X[genes])

    # z-score per gene (column-wise)
    Xz = (Xg - Xg.mean(axis=0)) / (Xg.std(axis=0).replace(0, np.nan))
    Xz = Xz.fillna(0.0)

    # sample order: High then Average
    order = yy.sort_values(ascending=False).index  # High first
    Xz = Xz.loc[order]
    yy = yy.loc[order]

    plt.figure(figsize=(11, 7))
    plt.imshow(Xz.values, aspect="auto")
    plt.colorbar(label="Z-score (log1p norm counts)")

    plt.title("Heatmap of Significant Genes (padj<0.05 & |log2FC|>1)")
    plt.xlabel("Genes")
    plt.ylabel("Samples (ordered by group)")

    # simple group separator line
    n_high = int((yy == "High").sum())
    plt.axhline(n_high - 0.5, linewidth=2)

    # hide dense ticks
    plt.xticks([])
    plt.yticks([])

    plt.tight_layout()
    _savefig("heatmap_deg_clustermap.png")


def main():
    # load labels
    y = _load_master_labels(MASTER_PATH)

    # load normalized counts (rows=samples, cols=genes)
    norm = pd.read_csv(NORM_PATH, index_col=0)
    norm = norm.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # load DE results for volcano + DEG list
    de = pd.read_csv(DE_PATH, index_col=0)

    # 1) Volcano
    volcano_plot(DE_PATH)

    # 2) PCA all genes
    pca_plot(
        norm_counts=norm,
        y=y,
        out_name="pca_all_genes.png",
        title="PCA of Gene Expression (All Genes)",
    )

    # DEG list from your threshold
    deg = de.dropna(subset=["padj", "log2FoldChange"])
    deg = deg[(deg["padj"] < 0.05) & (deg["log2FoldChange"].abs() > 1.0)]
    deg_genes = deg.index.astype(str).tolist()

    # 3) PCA significant genes
    if len(deg_genes) >= 2:
        pca_plot(
            norm_counts=norm[deg_genes].copy(),
            y=y,
            out_name="pca_significant_genes.png",
            title="PCA of Gene Expression (Significant Genes)",
        )
    else:
        print("[skip] PCA significant genes: not enough significant genes.")

    # 4) Heatmap (significant genes)
    if len(deg_genes) >= 2:
        heatmap_deg(norm, y, deg_genes)
    else:
        print("[skip] Heatmap: not enough significant genes.")

    print("\nDone. Check figures/ folder.")


if __name__ == "__main__":
    main()