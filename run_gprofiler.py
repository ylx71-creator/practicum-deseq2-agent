from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "output" / "deseq2_significant_genes.csv"
OUTPUT_DIR = BASE_DIR / "output" / "gprofiler"

INPUT_GENES_OUT = OUTPUT_DIR / "gprofiler_input_genes.csv"
RECOGNIZED_GENES_OUT = OUTPUT_DIR / "gprofiler_recognized_genes.csv"
UNRECOGNIZED_GENES_OUT = OUTPUT_DIR / "gprofiler_unrecognized_genes.csv"
GO_BP_OUT = OUTPUT_DIR / "gprofiler_GO_BP_results.csv"
ALL_SOURCES_OUT = OUTPUT_DIR / "gprofiler_all_sources_results.csv"


def load_clean_gene_list(input_file: Path) -> list[str]:
    """Load gene_symbol values and return unique, non-empty gene symbols."""
    if not input_file.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    print(f"Reading DESeq2 significant genes from: {input_file}")
    df = pd.read_csv(input_file)

    if "gene_symbol" not in df.columns:
        raise ValueError("The input CSV is missing the required 'gene_symbol' column.")

    if df.columns[0] != "gene_symbol":
        raise ValueError(
            "The first column of the input CSV must be 'gene_symbol'. "
            f"Found first column: '{df.columns[0]}'"
        )

    genes = (
        df["gene_symbol"]
        .dropna()
        .astype(str)
        .str.strip()
    )
    genes = genes[genes != ""]
    genes = genes.drop_duplicates()

    clean_genes = genes.tolist()
    if not clean_genes:
        raise ValueError("No valid genes found after removing missing values, empty strings, and duplicates.")

    return clean_genes


def run_enrichment(gp, genes: list[str], sources: list[str], output_file: Path, description: str) -> pd.DataFrame:
    print()
    print(f"Running g:Profiler enrichment analysis: {description}")
    print(f"Sources: {sources}")

    results = gp.profile(
        organism="hsapiens",
        query=genes,
        sources=sources,
    )

    results.to_csv(output_file, index=False)
    print(f"Saved results to: {output_file}")
    print(f"Number of enriched terms/pathways returned: {len(results)}")
    return results


def check_gene_recognition(gp, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    print()
    print("Checking which cleaned input genes are recognized by g:Profiler.")
    print("Running: gp.convert(organism=\"hsapiens\", query=genes)")

    converted = gp.convert(organism="hsapiens", query=genes)

    if "incoming" not in converted.columns or "converted" not in converted.columns:
        raise ValueError(
            "Unexpected g:Profiler convert() output. "
            "Expected columns named 'incoming' and 'converted'."
        )

    converted_values = converted["converted"].astype(str).str.strip()
    recognized = converted[
        converted["incoming"].notna()
        & converted["converted"].notna()
        & (converted_values != "")
        & (~converted_values.str.lower().isin({"none", "nan"}))
    ].copy()

    recognized_gene_set = set(recognized["incoming"].astype(str))
    unrecognized_genes = [gene for gene in genes if gene not in recognized_gene_set]
    unrecognized = pd.DataFrame({"gene_symbol": unrecognized_genes})

    recognized.to_csv(RECOGNIZED_GENES_OUT, index=False)
    unrecognized.to_csv(UNRECOGNIZED_GENES_OUT, index=False)

    print(f"Total cleaned input genes: {len(genes)}")
    print(f"Number of recognized genes: {len(recognized_gene_set)}")
    print(f"Number of unrecognized genes: {len(unrecognized_genes)}")
    print(f"Saved recognized genes to: {RECOGNIZED_GENES_OUT}")
    print(f"Saved unrecognized genes to: {UNRECOGNIZED_GENES_OUT}")

    return recognized, unrecognized


def main() -> None:
    print("Starting g:Profiler enrichment analysis workflow.")
    print("This script uses gene symbols from the first column of the DESeq2 significant genes CSV.")

    try:
        from gprofiler import GProfiler
    except ImportError as exc:
        raise ImportError(
            "The 'gprofiler-official' package is required. "
            "Install it with: pip install gprofiler-official"
        ) from exc

    genes = load_clean_gene_list(INPUT_FILE)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene_symbol": genes}).to_csv(INPUT_GENES_OUT, index=False)

    print(f"Created output folder: {OUTPUT_DIR}")
    print(f"Valid unique genes used for g:Profiler: {len(genes)}")
    print(f"Saved cleaned input gene list to: {INPUT_GENES_OUT}")

    gp = GProfiler(return_dataframe=True)

    check_gene_recognition(gp, genes)

    run_enrichment(
        gp=gp,
        genes=genes,
        sources=["GO:BP"],
        output_file=GO_BP_OUT,
        description="Version A - GO Biological Process only",
    )

    run_enrichment(
        gp=gp,
        genes=genes,
        sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
        output_file=ALL_SOURCES_OUT,
        description=(
            "Version B - GO Biological Process, Molecular Function, Cellular Component, "
            "KEGG pathways, and Reactome pathways"
        ),
    )

    print()
    print("Summary:")
    print("- GO:BP result = biological process enrichment only.")
    print(
        "- All sources result = GO biological process, molecular function, cellular component, "
        "KEGG pathways, and Reactome pathways."
    )
    print("- Each output CSV contains enriched terms/pathways, p-values, and matched genes.")
    print("g:Profiler enrichment analysis workflow complete.")


if __name__ == "__main__":
    main()
