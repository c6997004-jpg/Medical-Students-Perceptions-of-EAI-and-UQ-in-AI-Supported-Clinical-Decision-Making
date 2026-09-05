import pandas as pd
import numpy as np
from scipy.stats import kruskal
import scikit_posthocs as sp
import pingouin as pg


file = "DataA.xlsx"
groups = ["G1", "G2", "G3", "G4"]


posthoc_outcome = "Reliability Score"

posthoc_comparisons = [
    ("G1", "G2"),
    ("G1", "G3"),
    ("G1", "G4"),
    ("G2", "G3"),
    ("G2", "G4"),
    ("G3", "G4"),
]

df_list = []
for g in groups:
    temp = pd.read_excel(file, sheet_name=g)
    temp["Group"] = g
    df_list.append(temp)

df = pd.concat(df_list, ignore_index=True)


cols = [
    "Positive Affect Score",
    "Negative Affect Score",
    "Reliability Score",
    "Technical Competence Score",
    "Understandability Score",
    "Faith Score",
    "Attachment Score",
    "Overall HCT Score",
    "Usability Score",
]


# Cliff's delta is a rank-based pairwise effect size for independent groups.
# Positive values indicate that scores tend to be higher in x than in y.
def cliffs_delta(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(x) == 0 or len(y) == 0:
        return np.nan

    pairwise_differences = x[:, None] - y[None, :]
    greater = np.sum(pairwise_differences > 0)
    less = np.sum(pairwise_differences < 0)
    return (greater - less) / (len(x) * len(y))

kw_p_list = []
kw_results = {}
report_dict = {}


for c in cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")


    data = [df[df["Group"] == g][c].dropna() for g in groups]

    desc_list = []
    for g, d in zip(groups, data):
        median_val = d.median()
        q1, q3 = np.percentile(d, [25, 75])
        desc_list.append({
            "Group": g,
            "N": len(d),
            "Mean": d.mean(),
            "SD": d.std(),
            "Min": d.min(),
            "Max": d.max(),
            "Median": median_val,
            "IQR": q3 - q1,
        })
    desc_df = pd.DataFrame(desc_list)


    H, p = kruskal(*data)
    N_total = len(df[c].dropna())
    k = len(groups)  # 组数

    eta_squared_h = max(0.0, (H - k + 1) / (N_total - k))

    kw_p_list.append(p)
    kw_results[c] = {
        "H": H,
        "df": k - 1,
        "p": p,
        "eta_squared_H": eta_squared_h,
    }


    report_dict[c] = {"Descriptive": desc_df, "KW": kw_results[c]}


_, fdr_p = pg.multicomp(kw_p_list, method="fdr_bh")
kw_fdr_df = pd.DataFrame({
    "Variable": cols,
    "KW_p": kw_p_list,
    "KW_fdr_p": fdr_p,
    "eta_squared_H": [kw_results[v]["eta_squared_H"] for v in cols],
})


reliability_fdr_p = kw_fdr_df.loc[
    kw_fdr_df["Variable"] == posthoc_outcome,
    "KW_fdr_p",
].iloc[0]

if reliability_fdr_p < 0.05:
    reliability_df = df[["Group", posthoc_outcome]].dropna().copy()

    dunn_raw = sp.posthoc_dunn(
        reliability_df,
        val_col=posthoc_outcome,
        group_col="Group",
        p_adjust=None,
    )
    dunn_fdr = sp.posthoc_dunn(
        reliability_df,
        val_col=posthoc_outcome,
        group_col="Group",
        p_adjust="fdr_bh",
    )

    posthoc_rows = []
    for group1, group2 in posthoc_comparisons:
        scores1 = reliability_df.loc[
            reliability_df["Group"] == group1,
            posthoc_outcome,
        ].to_numpy()
        scores2 = reliability_df.loc[
            reliability_df["Group"] == group2,
            posthoc_outcome,
        ].to_numpy()

        posthoc_rows.append({
            "Comparison": f"{group1} - {group2}",
            "Group1": group1,
            "Group2": group2,
            "N_Group1": len(scores1),
            "N_Group2": len(scores2),
            "Cliffs_delta": cliffs_delta(scores1, scores2),
            "Dunn_p": dunn_raw.loc[group1, group2],
            "Dunn_FDR_p": dunn_fdr.loc[group1, group2],
        })

    reliability_posthoc = pd.DataFrame(posthoc_rows)
else:
    reliability_posthoc = pd.DataFrame(columns=[
        "Comparison",
        "Group1",
        "Group2",
        "N_Group1",
        "N_Group2",
        "Cliffs_delta",
        "Dunn_p",
        "Dunn_FDR_p",
    ])

analysis_notes = pd.DataFrame({
    "Item": [
        "Omnibus FDR family",
        "Post hoc gate",
        "Dunn FDR family",
        "Pairwise effect size",
        "Bootstrap confidence intervals",
    ],

})


with pd.ExcelWriter("Stat_Report.xlsx") as writer:
    kw_fdr_df.to_excel(writer, sheet_name="KruskalFDR", index=False)
    reliability_posthoc.to_excel(writer, sheet_name="Reliability_Dunn", index=False)
    analysis_notes.to_excel(writer, sheet_name="Analysis_Notes", index=False)

    # 每个变量单独 sheet
    for var in cols:
        desc_df = report_dict[var]["Descriptive"]

        startrow = 0
        desc_df.to_excel(writer, sheet_name=var, index=False, startrow=startrow)

        # KW 统计量写在最下方
        kw_stats = pd.DataFrame([report_dict[var]["KW"]])
        startrow += len(desc_df) + 3
        kw_stats.to_excel(writer, sheet_name=var, index=False, startrow=startrow)
