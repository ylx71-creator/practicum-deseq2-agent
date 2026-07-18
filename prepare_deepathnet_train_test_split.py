#!/usr/bin/env python3
"""Create stratified DeePathNet train/test files from prepared data and labels."""

from pathlib import Path

import pandas as pd
from sklearn.model_selection import train_test_split


DATA_PATH = Path("output/deepathnet_data_file.csv")
LABEL_PATH = Path("output/deepathnet_label_file.csv")
OUTPUT_DIR = Path("output/deepathnet")

TRAIN_DATA_PATH = OUTPUT_DIR / "train_data_file.csv"
TEST_DATA_PATH = OUTPUT_DIR / "test_data_file.csv"
TRAIN_LABEL_PATH = OUTPUT_DIR / "train_label_file.csv"
TEST_LABEL_PATH = OUTPUT_DIR / "test_label_file.csv"

RANDOM_STATE = 42
TEST_SIZE = 0.20


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is missing: {path}")


def require_columns(frame: pd.DataFrame, path: Path, required_columns: list[str]) -> None:
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(missing)}"
        )


def fail_if_duplicates(frame: pd.DataFrame, path: Path, column: str) -> None:
    duplicate_count = frame[column].duplicated().sum()
    if duplicate_count:
        raise ValueError(
            f"{path} contains {duplicate_count} duplicate value(s) in {column}"
        )


def print_distribution(title: str, labels: pd.Series) -> None:
    print(title)
    print(labels.value_counts(dropna=False).to_string())


def main() -> None:
    require_file(DATA_PATH)
    require_file(LABEL_PATH)

    data = pd.read_csv(DATA_PATH, dtype={"Cell_line": str})
    labels = pd.read_csv(LABEL_PATH, dtype={"Cell_line": str, "Cancer_type": str})

    require_columns(data, DATA_PATH, ["Cell_line"])
    require_columns(labels, LABEL_PATH, ["Cell_line", "Cancer_type"])
    fail_if_duplicates(data, DATA_PATH, "Cell_line")
    fail_if_duplicates(labels, LABEL_PATH, "Cell_line")

    if labels["Cancer_type"].isna().any():
        missing_count = labels["Cancer_type"].isna().sum()
        raise ValueError(f"{LABEL_PATH} contains {missing_count} missing Cancer_type values")

    print(f"Original data file samples: {len(data)}")
    print(f"Original label file samples: {len(labels)}")

    aligned = data[["Cell_line"]].merge(labels, on="Cell_line", how="inner")
    print(f"Overlapping samples: {len(aligned)}")
    print_distribution("Original label distribution:", aligned["Cancer_type"])

    if aligned.empty:
        raise ValueError("No overlapping Cell_line values between data and label files")

    class_counts = aligned["Cancer_type"].value_counts()
    if len(class_counts) < 2:
        raise ValueError("Stratified split requires at least two Cancer_type classes")

    too_small = class_counts[class_counts < 2]
    if not too_small.empty:
        raise ValueError(
            "Stratified split requires at least two samples per class. "
            f"Too-small class counts: {too_small.to_dict()}"
        )

    train_samples, test_samples = train_test_split(
        aligned["Cell_line"],
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=aligned["Cancer_type"],
    )

    train_sample_set = set(train_samples)
    test_sample_set = set(test_samples)

    train_data = data[data["Cell_line"].isin(train_sample_set)].copy()
    test_data = data[data["Cell_line"].isin(test_sample_set)].copy()
    train_labels = labels[labels["Cell_line"].isin(train_sample_set)].copy()
    test_labels = labels[labels["Cell_line"].isin(test_sample_set)].copy()

    train_order = {sample_id: index for index, sample_id in enumerate(train_samples)}
    test_order = {sample_id: index for index, sample_id in enumerate(test_samples)}

    train_data["sample_order"] = train_data["Cell_line"].map(train_order)
    test_data["sample_order"] = test_data["Cell_line"].map(test_order)
    train_labels["sample_order"] = train_labels["Cell_line"].map(train_order)
    test_labels["sample_order"] = test_labels["Cell_line"].map(test_order)

    train_data = train_data.sort_values("sample_order").drop(columns=["sample_order"])
    test_data = test_data.sort_values("sample_order").drop(columns=["sample_order"])
    train_labels = train_labels.sort_values("sample_order").drop(columns=["sample_order"])
    test_labels = test_labels.sort_values("sample_order").drop(columns=["sample_order"])

    train_labels = train_labels[["Cell_line", "Cancer_type"]]
    test_labels = test_labels[["Cell_line", "Cancer_type"]]

    train_features_match_test = list(train_data.columns) == list(test_data.columns)
    train_ids_align = train_data["Cell_line"].tolist() == train_labels["Cell_line"].tolist()
    test_ids_align = test_data["Cell_line"].tolist() == test_labels["Cell_line"].tolist()

    if not train_features_match_test:
        raise ValueError("Train/test feature columns do not match exactly")
    if not train_ids_align:
        raise ValueError("Train data Cell_line order does not match train labels")
    if not test_ids_align:
        raise ValueError("Test data Cell_line order does not match test labels")

    print(f"Train data shape: {train_data.shape}")
    print(f"Test data shape: {test_data.shape}")
    print(f"Train label shape: {train_labels.shape}")
    print(f"Test label shape: {test_labels.shape}")
    print_distribution("Train label distribution:", train_labels["Cancer_type"])
    print_distribution("Test label distribution:", test_labels["Cancer_type"])
    print(f"Train/test feature columns match exactly: {train_features_match_test}")
    print(f"Train sample IDs align with train labels: {train_ids_align}")
    print(f"Test sample IDs align with test labels: {test_ids_align}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    train_data.to_csv(TRAIN_DATA_PATH, index=False)
    test_data.to_csv(TEST_DATA_PATH, index=False)
    train_labels.to_csv(TRAIN_LABEL_PATH, index=False)
    test_labels.to_csv(TEST_LABEL_PATH, index=False)

    print(f"Saved train data to: {TRAIN_DATA_PATH}")
    print(f"Saved test data to: {TEST_DATA_PATH}")
    print(f"Saved train labels to: {TRAIN_LABEL_PATH}")
    print(f"Saved test labels to: {TEST_LABEL_PATH}")


if __name__ == "__main__":
    main()
