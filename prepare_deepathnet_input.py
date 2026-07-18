#!/usr/bin/env python3
"""Prepare a DeePathNet-ready RNA feature matrix from local practicum files."""

from pathlib import Path

import pandas as pd


ANNOTATION_PATH = Path("data/Human.GRCh38.p13.annot.tsv")
COUNTS_PATH = Path("data/GSE164641_raw_counts_GRCh38.p13_NCBI.tsv")
SIGNIFICANT_GENES_PATH = Path("output/deseq2_significant_genes.csv")
OUTPUT_DIR = Path("output")
OUTPUT_PATH = OUTPUT_DIR / "deepathnet_data_file.csv"


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is missing: {path}")


def require_columns(frame: pd.DataFrame, path: Path, required_columns: list[str]) -> None:
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(missing)}"
        )


def main() -> None:
    for path in (ANNOTATION_PATH, COUNTS_PATH, SIGNIFICANT_GENES_PATH):
        require_file(path)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    annotation = pd.read_csv(
        ANNOTATION_PATH,
        sep="\t",
        dtype={"GeneID": str},
        low_memory=False,
    )
    counts = pd.read_csv(COUNTS_PATH, sep="\t", dtype={"GeneID": str})
    significant_genes = pd.read_csv(SIGNIFICANT_GENES_PATH, dtype={"gene_symbol": str})

    require_columns(annotation, ANNOTATION_PATH, ["GeneID", "Symbol"])
    require_columns(counts, COUNTS_PATH, ["GeneID"])
    require_columns(significant_genes, SIGNIFICANT_GENES_PATH, ["gene_symbol"])

    sample_columns = [column for column in counts.columns if column != "GeneID"]
    if not sample_columns:
        raise ValueError(f"{COUNTS_PATH} must contain sample columns after GeneID")

    print(f"Annotation input shape: {annotation.shape}")
    print(f"Raw counts input shape: {counts.shape}")
    print(f"Significant genes input shape: {significant_genes.shape}")
    print(f"Detected {len(sample_columns)} sample columns in raw counts")

    gene_symbols = set(significant_genes["gene_symbol"].dropna())
    if not gene_symbols:
        raise ValueError(
            f"{SIGNIFICANT_GENES_PATH} does not contain any non-missing gene_symbol values"
        )

    mapped_counts = counts.merge(
        annotation[["GeneID", "Symbol"]],
        on="GeneID",
        how="left",
    )
    print(f"Merged shape after GeneID to Symbol mapping: {mapped_counts.shape}")
    print(f"Rows with missing Symbol after mapping: {mapped_counts['Symbol'].isna().sum()}")

    filtered_counts = mapped_counts.dropna(subset=["Symbol"])
    filtered_counts = filtered_counts[filtered_counts["Symbol"].isin(gene_symbols)]
    print(f"Filtered shape before duplicate Symbol removal: {filtered_counts.shape}")

    duplicate_symbol_count = filtered_counts["Symbol"].duplicated().sum()
    filtered_counts = filtered_counts.drop_duplicates(subset=["Symbol"], keep="first")
    print(f"Duplicate Symbols removed: {duplicate_symbol_count}")
    print(f"Filtered shape after duplicate Symbol removal: {filtered_counts.shape}")

    if filtered_counts.empty:
        raise ValueError(
            "No rows remain after filtering counts to significant gene symbols"
        )

    feature_matrix = filtered_counts.set_index("Symbol")
    feature_matrix = feature_matrix.drop(columns=["GeneID"])

    deepathnet_data = feature_matrix.transpose().reset_index()
    deepathnet_data = deepathnet_data.rename(columns={"index": "Cell_line"})
    deepathnet_data = deepathnet_data.rename(
        columns={
            column: f"{column}_RNA"
            for column in deepathnet_data.columns
            if column != "Cell_line"
        }
    )
    deepathnet_data.columns.name = None

    print(f"Final DeePathNet output shape: {deepathnet_data.shape}")
    print("First few rows of final DeePathNet file:")
    print(deepathnet_data.head())

    deepathnet_data.to_csv(OUTPUT_PATH, index=False)
    print(f"Saved DeePathNet-ready feature file to: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
