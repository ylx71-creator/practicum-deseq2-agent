#!/usr/bin/env python3
"""Train baseline classifiers on the prepared DeePathNet risk-group split."""

from pathlib import Path

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


TRAIN_DATA_PATH = Path("output/deepathnet/train_data_file.csv")
TEST_DATA_PATH = Path("output/deepathnet/test_data_file.csv")
TRAIN_LABEL_PATH = Path("output/deepathnet/train_label_file.csv")
TEST_LABEL_PATH = Path("output/deepathnet/test_label_file.csv")

OUTPUT_DIR = Path("output/baseline")
METRICS_PATH = OUTPUT_DIR / "baseline_metrics.csv"

LABEL_TO_INT = {"Average": 0, "High": 1}
INT_TO_LABEL = {0: "Average", 1: "High"}
RANDOM_STATE = 42


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is missing: {path}")


def require_columns(frame: pd.DataFrame, path: Path, required_columns: list[str]) -> None:
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(missing)}"
        )


def fail_if_duplicate_cell_lines(frame: pd.DataFrame, path: Path) -> None:
    duplicate_count = frame["Cell_line"].duplicated().sum()
    if duplicate_count:
        raise ValueError(
            f"{path} contains {duplicate_count} duplicate Cell_line value(s)"
        )


def read_and_validate_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in (TRAIN_DATA_PATH, TEST_DATA_PATH, TRAIN_LABEL_PATH, TEST_LABEL_PATH):
        require_file(path)

    train_data = pd.read_csv(TRAIN_DATA_PATH, dtype={"Cell_line": str})
    test_data = pd.read_csv(TEST_DATA_PATH, dtype={"Cell_line": str})
    train_labels = pd.read_csv(
        TRAIN_LABEL_PATH,
        dtype={"Cell_line": str, "Cancer_type": str},
    )
    test_labels = pd.read_csv(
        TEST_LABEL_PATH,
        dtype={"Cell_line": str, "Cancer_type": str},
    )

    require_columns(train_data, TRAIN_DATA_PATH, ["Cell_line"])
    require_columns(test_data, TEST_DATA_PATH, ["Cell_line"])
    require_columns(train_labels, TRAIN_LABEL_PATH, ["Cell_line", "Cancer_type"])
    require_columns(test_labels, TEST_LABEL_PATH, ["Cell_line", "Cancer_type"])

    for frame, path in (
        (train_data, TRAIN_DATA_PATH),
        (test_data, TEST_DATA_PATH),
        (train_labels, TRAIN_LABEL_PATH),
        (test_labels, TEST_LABEL_PATH),
    ):
        fail_if_duplicate_cell_lines(frame, path)

    if list(train_data.columns) != list(test_data.columns):
        raise ValueError("Train and test data feature schemas do not match exactly")

    train_unknown = set(train_labels["Cancer_type"].dropna()) - set(LABEL_TO_INT)
    test_unknown = set(test_labels["Cancer_type"].dropna()) - set(LABEL_TO_INT)
    if train_unknown:
        raise ValueError(f"Unknown train Cancer_type label(s): {sorted(train_unknown)}")
    if test_unknown:
        raise ValueError(f"Unknown test Cancer_type label(s): {sorted(test_unknown)}")
    if train_labels["Cancer_type"].isna().any() or test_labels["Cancer_type"].isna().any():
        raise ValueError("Cancer_type contains missing values")

    return train_data, test_data, train_labels, test_labels


def align_data_and_labels(
    data: pd.DataFrame,
    labels: pd.DataFrame,
    split_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    aligned = data.merge(labels, on="Cell_line", how="inner", validate="one_to_one")
    if len(aligned) != len(data) or len(aligned) != len(labels):
        raise ValueError(
            f"{split_name} data/label alignment failed: "
            f"{len(data)} data rows, {len(labels)} label rows, {len(aligned)} overlaps"
        )

    feature_columns = [column for column in data.columns if column != "Cell_line"]
    aligned_data = aligned[["Cell_line"] + feature_columns].copy()
    aligned_labels = aligned[["Cell_line", "Cancer_type"]].copy()
    return aligned_data, aligned_labels


def predicted_probabilities(model, x_test: pd.DataFrame) -> tuple[list[float], list[float]]:
    if hasattr(model, "predict_proba"):
        probabilities = model.predict_proba(x_test)
        classes = list(model.classes_)
        average_index = classes.index(0)
        high_index = classes.index(1)
        return probabilities[:, average_index].tolist(), probabilities[:, high_index].tolist()

    predictions = model.predict(x_test)
    prob_high = [float(value) for value in predictions]
    prob_average = [1.0 - value for value in prob_high]
    return prob_average, prob_high


def evaluate_model(
    model_name: str,
    model,
    x_train: pd.DataFrame,
    y_train: pd.Series,
    x_test: pd.DataFrame,
    y_test: pd.Series,
    test_cell_lines: pd.Series,
    prediction_path: Path,
) -> dict[str, object]:
    model.fit(x_train, y_train)
    y_pred = model.predict(x_test)
    prob_average, prob_high = predicted_probabilities(model, x_test)

    tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0, 1]).ravel()
    predicted_counts = pd.Series(y_pred).map(INT_TO_LABEL).value_counts()
    predicted_average = int(predicted_counts.get("Average", 0))
    predicted_high = int(predicted_counts.get("High", 0))

    try:
        roc_auc = roc_auc_score(y_test, prob_high)
    except ValueError:
        roc_auc = float("nan")

    metrics = {
        "model": model_name,
        "accuracy": accuracy_score(y_test, y_pred),
        "precision": precision_score(y_test, y_pred, zero_division=0),
        "recall": recall_score(y_test, y_pred, zero_division=0),
        "f1": f1_score(y_test, y_pred, zero_division=0),
        "roc_auc": roc_auc,
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
        "predicted_Average": predicted_average,
        "predicted_High": predicted_high,
    }

    predictions = pd.DataFrame(
        {
            "Cell_line": test_cell_lines,
            "y_true": y_test.map(INT_TO_LABEL),
            "y_pred": pd.Series(y_pred).map(INT_TO_LABEL),
            "prob_Average": prob_average,
            "prob_High": prob_high,
        }
    )
    predictions.to_csv(prediction_path, index=False)

    print(f"\n{model_name} metrics:")
    for key, value in metrics.items():
        if key != "model":
            print(f"  {key}: {value}")
    if predicted_average == 0 or predicted_high == 0:
        print(f"  Warning: {model_name} predicted only one class")
    else:
        print(f"  {model_name} predicted both classes")
    print(f"  Saved predictions to: {prediction_path}")

    return metrics


def make_xgboost_model():
    try:
        from xgboost import XGBClassifier
    except ImportError:
        print("\nXGBoost is not installed; skipping XGBoost baseline.")
        return None

    return XGBClassifier(
        n_estimators=200,
        max_depth=3,
        learning_rate=0.05,
        subsample=0.9,
        colsample_bytree=0.9,
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=RANDOM_STATE,
        n_jobs=1,
        scale_pos_weight=1.0,
    )


def main() -> None:
    train_data, test_data, train_labels, test_labels = read_and_validate_inputs()
    train_data, train_labels = align_data_and_labels(train_data, train_labels, "Train")
    test_data, test_labels = align_data_and_labels(test_data, test_labels, "Test")

    feature_columns = [column for column in train_data.columns if column != "Cell_line"]
    if not feature_columns:
        raise ValueError("No feature columns found after excluding Cell_line")

    x_train = train_data[feature_columns]
    x_test = test_data[feature_columns]
    y_train = train_labels["Cancer_type"].map(LABEL_TO_INT)
    y_test = test_labels["Cancer_type"].map(LABEL_TO_INT)

    print(f"Train data shape: {train_data.shape}")
    print(f"Test data shape: {test_data.shape}")
    print(f"Feature count: {len(feature_columns)}")
    print("Train label distribution:")
    print(train_labels["Cancer_type"].value_counts().to_string())
    print("Test label distribution:")
    print(test_labels["Cancer_type"].value_counts().to_string())
    print(
        "Train Cell_line alignment: "
        f"{train_data['Cell_line'].tolist() == train_labels['Cell_line'].tolist()}"
    )
    print(
        "Test Cell_line alignment: "
        f"{test_data['Cell_line'].tolist() == test_labels['Cell_line'].tolist()}"
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    models = [
        (
            "logistic_regression",
            Pipeline(
                steps=[
                    ("standardize", StandardScaler()),
                    (
                        "classifier",
                        LogisticRegression(
                            class_weight="balanced",
                            max_iter=5000,
                            random_state=RANDOM_STATE,
                        ),
                    ),
                ]
            ),
            OUTPUT_DIR / "logistic_predictions.csv",
        ),
        (
            "random_forest",
            RandomForestClassifier(
                n_estimators=500,
                class_weight="balanced",
                random_state=RANDOM_STATE,
                n_jobs=-1,
            ),
            OUTPUT_DIR / "random_forest_predictions.csv",
        ),
    ]

    xgboost_model = make_xgboost_model()
    if xgboost_model is not None:
        models.append(
            (
                "xgboost",
                xgboost_model,
                OUTPUT_DIR / "xgboost_predictions.csv",
            )
        )

    metrics = []
    for model_name, model, prediction_path in models:
        metrics.append(
            evaluate_model(
                model_name=model_name,
                model=model,
                x_train=x_train,
                y_train=y_train,
                x_test=x_test,
                y_test=y_test,
                test_cell_lines=test_data["Cell_line"],
                prediction_path=prediction_path,
            )
        )

    metrics_table = pd.DataFrame(metrics)
    metrics_table.to_csv(METRICS_PATH, index=False)
    print(f"\nSaved baseline metrics to: {METRICS_PATH}")


if __name__ == "__main__":
    main()
