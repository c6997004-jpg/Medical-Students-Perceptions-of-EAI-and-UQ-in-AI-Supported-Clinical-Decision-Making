"""Bootstrap 95% CIs for the raw-scale mean differences reported for RQ1 and RQ2."""

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "DataA.xlsx"
OUTPUT_FILE = BASE_DIR / "RQ1_RQ2_Bootstrap_Mean_Differences.xlsx"

GROUPS = ["G1", "G2", "G3", "G4"]
STAGE_COLUMN = "Preclinical vs clinical"
PRECLINICAL_CODE = "P"
CLINICAL_CODE = "C"

RQ1_OUTCOMES = ["Reliability Score"]
RQ2_OUTCOMES = ["Reliability Score", "Overall HCT Score"]

N_RESAMPLES = 10_000
RQ1_RANDOM_SEED = 42
RQ2_RANDOM_SEED = 20260709


def load_data():
    frames = []

    for group in GROUPS:
        group_data = pd.read_excel(DATA_FILE, sheet_name=group)
        group_data["Group"] = group
        frames.append(group_data)

    data = pd.concat(frames, ignore_index=True)
    data[STAGE_COLUMN] = (
        data[STAGE_COLUMN].astype(str).str.strip().str.upper()
    )

    required = set(RQ1_OUTCOMES + RQ2_OUTCOMES + [STAGE_COLUMN])
    missing = sorted(required.difference(data.columns))

    if missing:
        raise ValueError(f"Missing expected columns: {missing}")

    return data


def values_for(data, level_column, level, outcome):
    """Return nonmissing outcome values for one comparison level."""

    return (
        pd.to_numeric(
            data.loc[data[level_column] == level, outcome],
            errors="coerce",
        )
        .dropna()
        .to_numpy(dtype=float)
    )


def bootstrap_mean_difference(x, y, seed):
    """Calculate x minus y and its percentile bootstrap 95% CI."""

    if len(x) == 0 or len(y) == 0:
        raise ValueError(
            "Both comparison levels must contain nonmissing observations."
        )

    rng = np.random.default_rng(seed)

    x_samples = x[
        rng.integers(0, len(x), size=(N_RESAMPLES, len(x)))
    ]
    y_samples = y[
        rng.integers(0, len(y), size=(N_RESAMPLES, len(y)))
    ]

    bootstrap_differences = (
        x_samples.mean(axis=1) - y_samples.mean(axis=1)
    )

    difference = x.mean() - y.mean()
    ci_low, ci_high = np.percentile(
        bootstrap_differences,
        [2.5, 97.5],
    )

    return difference, ci_low, ci_high


def comparison_row(rq, outcome, level_1, level_2, x, y, seed):
    difference, ci_low, ci_high = bootstrap_mean_difference(
        x,
        y,
        seed,
    )

    return {
        "Research question": rq,
        "Outcome": outcome,
        "Comparison": f"{level_1} - {level_2}",
        "Level 1": level_1,
        "Level 1 n": len(x),
        "Level 1 mean": x.mean(),
        "Level 2": level_2,
        "Level 2 n": len(y),
        "Level 2 mean": y.mean(),
        "Mean difference": difference,
        "95% CI lower": ci_low,
        "95% CI upper": ci_high,
        "Mean difference (95% CI)": (
            f"{difference:.2f} "
            f"(95% CI {ci_low:.2f} to {ci_high:.2f})"
        ),
    }


def rq1_group_differences(data):
    """All six pairwise group mean differences for reliability."""

    rows = []

    for comparison_number, (group_1, group_2) in enumerate(
        combinations(GROUPS, 2)
    ):
        outcome = RQ1_OUTCOMES[0]

        x = values_for(data, "Group", group_1, outcome)
        y = values_for(data, "Group", group_2, outcome)

        rows.append(
            comparison_row(
                "RQ1",
                outcome,
                group_1,
                group_2,
                x,
                y,
                RQ1_RANDOM_SEED + 1000 + comparison_number,
            )
        )

    return pd.DataFrame(rows)


def rq2_training_stage_differences(data):
    """Preclinical-minus-clinical differences pooled across G1-G4."""

    rows = []

    for outcome_number, outcome in enumerate(RQ2_OUTCOMES):
        preclinical = values_for(
            data,
            STAGE_COLUMN,
            PRECLINICAL_CODE,
            outcome,
        )

        clinical = values_for(
            data,
            STAGE_COLUMN,
            CLINICAL_CODE,
            outcome,
        )

        rows.append(
            comparison_row(
                "RQ2",
                outcome,
                "Preclinical",
                "Clinical",
                preclinical,
                clinical,
                RQ2_RANDOM_SEED + outcome_number,
            )
        )

    return pd.DataFrame(rows)


def main():
    data = load_data()

    rq1_results = rq1_group_differences(data)
    rq2_results = rq2_training_stage_differences(data)

    with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
        rq1_results.to_excel(
            writer,
            sheet_name="RQ1_Group_Differences",
            index=False,
        )

        rq2_results.to_excel(
            writer,
            sheet_name="RQ2_Stage_Differences",
            index=False,
        )

    print(f"Saved results to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()