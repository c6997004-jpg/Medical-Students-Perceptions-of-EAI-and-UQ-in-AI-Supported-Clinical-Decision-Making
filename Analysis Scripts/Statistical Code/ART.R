# ============================================================
# ART ANOVA for RQ2
# Factors:
#   Group = AI presentation condition (G1-G4)
#   Stage = training stage (Preclinical vs Clinical)
# Outcomes:
#   9 questionnaire/domain scores
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ARTool)
library(tibble)
library(writexl)

options(contrasts = c("contr.sum", "contr.poly"))

file_path <- "DataA.xlsx"
groups <- c("G1", "G2", "G3", "G4")
stage_col <- "Preclinical vs clinical"
alpha <- 0.05

dvs <- c(
  "Positive Affect Score",
  "Negative Affect Score",
  "Reliability Score",
  "Technical Competence Score",
  "Understandability Score",
  "Faith Score",
  "Attachment Score",
  "Overall HCT Score",
  "Usability Score"
)

if (!file.exists(file_path)) {
  stop("Cannot find data file: ", file_path)
}

available_sheets <- excel_sheets(file_path)
missing_sheets <- setdiff(groups, available_sheets)
if (length(missing_sheets) > 0) {
  stop("Missing expected sheets: ", paste(missing_sheets, collapse = ", "))
}

df_all <- lapply(groups, function(g) {
  message("Reading sheet: ", g)
  tmp <- read_excel(path = file_path, sheet = g)
  tmp$Group <- g
  tmp
}) %>%
  bind_rows()

required_cols <- c(stage_col, dvs)
missing_cols <- setdiff(required_cols, names(df_all))
if (length(missing_cols) > 0) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}

df_all <- df_all %>%
  mutate(
    Group = factor(Group, levels = groups),
    Stage_raw = toupper(trimws(as.character(.data[[stage_col]]))),
    Stage_raw = case_when(
      Stage_raw %in% c("P", "PRECLINICAL", "PRE-CLINICAL") ~ "P",
      Stage_raw %in% c("C", "CLINICAL") ~ "C",
      TRUE ~ Stage_raw
    )
  ) %>%
  filter(Stage_raw %in% c("P", "C")) %>%
  mutate(
    Stage = factor(
      Stage_raw,
      levels = c("P", "C"),
      labels = c("Preclinical", "Clinical")
    )
  ) %>%
  select(-Stage_raw)

message("\nSample size by Group x Stage:")
print(table(df_all$Group, df_all$Stage))

partial_eta2_from_f <- function(f_value, df_effect, df_error) {
  (f_value * df_effect) / (f_value * df_effect + df_error)
}

run_art_for_dv <- function(dv_name, data) {
  message("\nAnalyzing: ", dv_name)
  
  dat <- data %>%
    select(Group, Stage, DV_value = all_of(dv_name)) %>%
    mutate(DV_value = as.numeric(DV_value)) %>%
    drop_na()
  
  message("  Valid N: ", nrow(dat))
  
  if (nrow(dat) == 0) {
    stop("No valid observations for outcome: ", dv_name)
  }
  
  art_model <- art(DV_value ~ Group * Stage, data = dat)
  
  as.data.frame(anova(art_model)) %>%
    rownames_to_column(var = "Effect") %>%
    mutate(
      DV = dv_name,
      partial_eta2 = partial_eta2_from_f(`F value`, Df, Df.res)
    )
}

art_results <- lapply(dvs, run_art_for_dv, data = df_all) %>%
  bind_rows()

effect_order <- c("Group", "Stage", "Group:Stage")
effect_labels <- c(
  Group = "AI presentation condition",
  Stage = "Training stage",
  `Group:Stage` = "Condition x Stage"
)

# BH/FDR is applied to the full family of omnibus tests:
# 9 outcomes x 3 effects = 27 p-values.
omnibus_results <- art_results %>%
  filter(Effect %in% effect_order) %>%
  mutate(
    Effect_Label = unname(effect_labels[Effect]),
    p_unc = `Pr(>F)`,
    p_adj = p.adjust(p_unc, method = "BH"),
    Significant = p_adj < alpha,
    DV_order = match(DV, dvs),
    Effect_order = match(Effect, effect_order)
  ) %>%
  arrange(DV_order, Effect_order) %>%
  select(
    DV,
    Effect,
    Effect_Label,
    Df,
    Df.res,
    `F value`,
    p_unc,
    p_adj,
    Significant,
    partial_eta2
  )

full_anova <- art_results %>%
  mutate(
    p_unc = `Pr(>F)`,
    partial_eta2 = partial_eta2_from_f(`F value`, Df, Df.res)
  ) %>%
  select(DV, Effect, Df, Df.res, `F value`, p_unc, partial_eta2)

sample_size <- df_all %>%
  group_by(Group, Stage) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = Stage, values_from = N, values_fill = 0)

console_table <- omnibus_results %>%
  mutate(
    `F value` = round(`F value`, 3),
    p_unc = round(p_unc, 3),
    p_adj = round(p_adj, 3),
    partial_eta2 = round(partial_eta2, 3)
  )

message("\nART ANOVA omnibus results:")
print(console_table)

notes <- tibble(
  Item = c(
    "Model",
    "Correction",
    "Effect size",
    "Significance criterion"
  ),
  Description = c(
    "DV ~ AI presentation condition * training stage, fitted with ARTool::art() for each outcome.",
    "Benjamini-Hochberg FDR correction applied across the 27 omnibus tests unless changed in the script.",
    "Partial eta squared computed from F, numerator df, and denominator df: F*df1 / (F*df1 + df2).",
    paste0("FDR-adjusted p < ", alpha)
  )
)

write_xlsx(
  list(
    ART_omnibus_FDR = omnibus_results,
    ART_full_anova_tables = full_anova,
    Sample_Size = sample_size,
    Notes = notes
  ),
  path = "PlanB_ART_Analysis_Revised.xlsx"
)

message("\nDone. Results saved to PlanB_ART_Analysis_Revised.xlsx")
