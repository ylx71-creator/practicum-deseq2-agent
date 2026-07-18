#!/usr/bin/env python3
"""Prepare the DeePathNet label file from the local GEO series matrix."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pandas as pd


SERIES_MATRIX_PATH = Path("data/GSE164641_series_matrix.txt")
FEATURE_FILE_PATH = Path("output/deepathnet_data_file.csv")
OUTPUT_DIR = Path("output")
OUTPUT_PATH = OUTPUT_DIR / "deepathnet_label_file.csv"

RISK_CATEGORY_PATTERN = re.compile(r"risk category:\s*(.+)", re.IGNORECASE)


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is missing: {path}")


def clean_cell(value: str) -> str:
    return value.strip().strip('"')


def read_series_metadata(path: Path) -> tuple[list[str], list[str], int]:
    sample_ids = None
    risk_labels = None
    risk_row_number = None

    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=1):
            if not row:
                continue

            row_name = clean_cell(row[0])
            values = [clean_cell(value) for value in row[1:]]

            if row_name == "!Sample_geo_accession":
                sample_ids = values

            if any(RISK_CATEGORY_PATTERN.search(value) for value in values):
                risk_labels = values
                risk_row_number = row_number

    if sample_ids is None:
        raise ValueError(
            f"Could not find !Sample_geo_accession row in {path}"
        )

    if risk_labels is None or risk_row_number is None:
        raise ValueError(
            f"Could not find a row containing risk category labels in {path}"
        )

    if len(sample_ids) != len(risk_labels):
        raise ValueError(
            "Sample ID count does not match risk label count: "
            f"{len(sample_ids)} sample IDs vs {len(risk_labels)} labels"
        )

    return sample_ids, risk_labels, risk_row_number


def extract_risk_category(value: str) -> str | None:
    match = RISK_CATEGORY_PATTERN.search(value)
    if not match:
        return None
    return match.group(1).strip()


def main() -> None:
    require_file(SERIES_MATRIX_PATH)
    require_file(FEATURE_FILE_PATH)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    sample_ids, raw_risk_labels, risk_row_number = read_series_metadata(
        SERIES_MATRIX_PATH
    )
    print(f"Found {len(sample_ids)} samples in the series matrix")
    print(f"Risk category labels found on row: {risk_row_number}")

    labels = [
        {
            "Cell_line": sample_id,
            "Cancer_type": extract_risk_category(raw_label),
        }
        for sample_id, raw_label in zip(sample_ids, raw_risk_labels)
    ]
    label_table = pd.DataFrame(labels)
    label_table = label_table.dropna(subset=["Cancer_type"])
    print(f"Extracted {len(label_table)} risk category labels")

    feature_data = pd.read_csv(FEATURE_FILE_PATH, usecols=["Cell_line"], dtype=str)
    feature_samples = feature_data["Cell_line"].dropna().tolist()
    feature_sample_set = set(feature_samples)

    overlapping_labels = label_table[
        label_table["Cell_line"].isin(feature_sample_set)
    ].copy()
    print(
        "Samples overlapping with DeePathNet feature file: "
        f"{len(overlapping_labels)}"
    )

    if overlapping_labels.empty:
        raise ValueError(
            "No series matrix labels overlap with output/deepathnet_data_file.csv"
        )

    sample_order = {sample_id: index for index, sample_id in enumerate(feature_samples)}
    overlapping_labels["sample_order"] = overlapping_labels["Cell_line"].map(sample_order)
    output = overlapping_labels.sort_values("sample_order")
    output = output[["Cell_line", "Cancer_type"]].reset_index(drop=True)

    print(f"Final DeePathNet label output shape: {output.shape}")
    print("First few rows of final DeePathNet label file:")
    print(output.head())

    output.to_csv(OUTPUT_PATH, index=False)
    print(f"Saved DeePathNet-ready label file to: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
