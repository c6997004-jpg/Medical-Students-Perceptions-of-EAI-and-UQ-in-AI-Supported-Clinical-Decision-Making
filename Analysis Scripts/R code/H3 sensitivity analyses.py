"""HC3 sensitivity analyses for Overall HCT and perceived reliability."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


SCRIPT_DIR = Path(__file__).resolve().parent
GROUPS = ["G1", "G2", "G3", "G4"]
OUTCOMES = {
    "Overall HCT": "Overall HCT Score",
    "Perceived reliability": "Reliability Score",
}
FORMULA = "Outcome ~ C(Group, Treatment(reference='G1'))"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run group-only HC3 regressions for the two principal RQ1 outcomes."
    )
    parser.add_argument(
        "--input",
        default=str(SCRIPT_DIR / "DataA.xlsx"),
        help="Input workbook (default: DataA.xlsx beside this script).",
    )
    parser.add_argument(
        "--output",
        default=str(SCRIPT_DIR / "RQ1_HC3_Sensitivity.xlsx"),
        help="Output workbook (default: RQ1_HC3_Sensitivity.xlsx beside this script).",
    )
    return parser.parse_args()


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = adjusted_ranked
    return adjusted.tolist()


def load_data(input_path: Path) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(f"Input workbook not found: {input_path}")

    frames = []
    for group in GROUPS:
        frame = pd.read_excel(input_path, sheet_name=group)
        missing = set(OUTCOMES.values()).difference(frame.columns)
        if missing:
            raise ValueError(f"Sheet {group} is missing columns: {sorted(missing)}")
        selected = frame[list(OUTCOMES.values())].copy()
        selected["Group"] = group
        frames.append(selected)

    data = pd.concat(frames, ignore_index=True)
    for column in OUTCOMES.values():
        data[column] = pd.to_numeric(data[column], errors="coerce")
    return data


def fit_outcome(data: pd.DataFrame, label: str, source_column: str):
    analysis = data[["Group", source_column]].dropna().copy()
    analysis = analysis.rename(columns={source_column: "Outcome"})
    analysis["Group"] = pd.Categorical(analysis["Group"], GROUPS, ordered=True)

    model = smf.ols(FORMULA, data=analysis).fit(cov_type="HC3", use_t=True)
    group_terms = [
        term
        for term in model.params.index
        if term.startswith("C(Group") and "[T." in term
    ]
    if len(group_terms) != 3:
        raise RuntimeError(f"Expected three group coefficients; found {group_terms}")

    restriction = ", ".join(f"{term} = 0" for term in group_terms)
    test = model.f_test(restriction)
    omnibus = {
        "Outcome": label,
        "N": int(model.nobs),
        "Robust F": float(np.asarray(test.fvalue).squeeze()),
        "df numerator": float(test.df_num),
        "df denominator": float(test.df_denom),
        "P value": float(np.asarray(test.pvalue).squeeze()),
    }

    ci = model.conf_int(alpha=0.05)
    contrasts = []
    for term in group_terms:
        comparison = term.split("[T.", 1)[1].rstrip("]") + " - G1"
        contrasts.append(
            {
                "Outcome": label,
                "Contrast": comparison,
                "Mean difference": float(model.params[term]),
                "HC3 SE": float(model.bse[term]),
                "95% CI low": float(ci.loc[term, 0]),
                "95% CI high": float(ci.loc[term, 1]),
                "t": float(model.tvalues[term]),
                "P value": float(model.pvalues[term]),
            }
        )

    descriptives = (
        analysis.groupby("Group", observed=False)["Outcome"]
        .agg(N="count", Mean="mean", SD="std", Median="median")
        .reset_index()
    )
    descriptives.insert(0, "Outcome", label)

    summary = {
        "Outcome": label,
        "Source column": source_column,
        "Formula": FORMULA,
        "N": int(model.nobs),
        "R-squared": float(model.rsquared),
        "Adjusted R-squared": float(model.rsquared_adj),
        "Covariance estimator": "HC3",
        "Reference group": "G1",
    }
    return summary, omnibus, contrasts, descriptives


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    if output_path.suffix.lower() != ".xlsx":
        raise ValueError("--output must end in .xlsx")

    data = load_data(input_path)
    summaries = []
    omnibus_rows = []
    contrast_rows = []
    descriptive_frames = []

    for label, source_column in OUTCOMES.items():
        summary, omnibus, contrasts, descriptives = fit_outcome(
            data, label, source_column
        )
        summaries.append(summary)
        omnibus_rows.append(omnibus)
        contrast_rows.extend(contrasts)
        descriptive_frames.append(descriptives)

    for row, adjusted in zip(
        omnibus_rows,
        benjamini_hochberg([row["P value"] for row in omnibus_rows]),
    ):
        row["BH-adjusted P value (2 omnibus tests)"] = adjusted

    for row, adjusted in zip(
        contrast_rows,
        benjamini_hochberg([row["P value"] for row in contrast_rows]),
    ):
        row["BH-adjusted P value (6 contrasts)"] = adjusted

    summary_df = pd.DataFrame(summaries)
    omnibus_df = pd.DataFrame(omnibus_rows)
    contrasts_df = pd.DataFrame(contrast_rows)
    descriptives_df = pd.concat(descriptive_frames, ignore_index=True)
    notes_df = pd.DataFrame(
        {
            "Item": [
                "Purpose",
                "Primary analysis",
                "Sensitivity model",
                "Principal outcomes",
                "Omnibus FDR family",
                "Contrast FDR family",
                "Excluded terms",
            ],
            "Value": [
                "Targeted RQ1 model-based sensitivity analysis for Reviewer H3",
                "Kruskal-Wallis/Dunn analyses remain primary",
                "Outcome ~ Group with HC3 heteroskedasticity-robust standard errors",
                "Overall HCT and perceived reliability",
                "Two robust omnibus Group tests",
                "Six G1-referenced contrasts (three per outcome)",
                "Recruitment site, training stage, and interactions",
            ],
        }
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="Model_Summary", index=False)
        omnibus_df.to_excel(writer, sheet_name="Robust_Omnibus", index=False)
        contrasts_df.to_excel(writer, sheet_name="Robust_Contrasts", index=False)
        descriptives_df.to_excel(writer, sheet_name="Descriptives", index=False)
        notes_df.to_excel(writer, sheet_name="Notes", index=False)

    print(omnibus_df.to_string(index=False))
    print(f"\nInput: {input_path}")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
