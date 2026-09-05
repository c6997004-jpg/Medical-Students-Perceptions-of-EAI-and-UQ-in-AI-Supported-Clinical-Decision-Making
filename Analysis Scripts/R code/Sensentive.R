# ============================================================
# Supplementary analyses for recruitment site
#
# Analysis A: Recruitment-site balance across randomized groups
# Analysis B1 (RQ1 sensitivity): Outcome ~ Group * Site
# Analysis B2 (RQ2 sensitivity): Outcome ~ Group * Stage * Site
# Analysis C: Descriptive outcome statistics by recruitment site
#
# Site denotes recruitment site (country of the medical school
# through which the participant was recruited), not nationality,
# citizenship, ethnicity, or individual cultural identity.
# ============================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ARTool",
  "tibble",
  "writexl",
  "emmeans"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ARTool)
  library(tibble)
  library(writexl)
  library(emmeans)
})

options(contrasts = c("contr.sum", "contr.poly"))

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

file_path <- "DataA.xlsx"
output_path <- "Supplementary_Site_Adjusted_ART_Analyses.xlsx"

groups <- c("G1", "G2", "G3", "G4")
stage_col <- "Preclinical vs clinical"
site_col <- "Recruitment site"
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
  stop(
    "Missing expected sheets: ",
    paste(missing_sheets, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Read and clean data
# ------------------------------------------------------------

df_all <- lapply(groups, function(g) {
  message("Reading sheet: ", g)
  tmp <- read_excel(path = file_path, sheet = g)
  tmp$Group <- g
  tmp
}) %>%
  bind_rows()

required_cols <- c(stage_col, site_col, dvs)
missing_cols <- setdiff(required_cols, names(df_all))

if (length(missing_cols) > 0) {
  stop(
    "Missing expected columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

df_all <- df_all %>%
  mutate(
    Group = factor(Group, levels = groups),
    
    Stage_raw = toupper(trimws(as.character(.data[[stage_col]]))),
    Stage_raw = case_when(
      Stage_raw %in% c("P", "PRECLINICAL", "PRE-CLINICAL") ~ "P",
      Stage_raw %in% c("C", "CLINICAL") ~ "C",
      TRUE ~ NA_character_
    ),
    
    Site_raw = toupper(trimws(as.character(.data[[site_col]]))),
    Site_raw = case_when(
      Site_raw %in% c("SINGAPORE", "SG") ~ "Singapore",
      Site_raw %in% c("CHINA", "CN") ~ "China",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    Stage_raw %in% c("P", "C"),
    !is.na(Site_raw)
  ) %>%
  mutate(
    Stage = factor(
      Stage_raw,
      levels = c("P", "C"),
      labels = c("Preclinical", "Clinical")
    ),
    Site = factor(
      Site_raw,
      levels = c("Singapore", "China")
    )
  ) %>%
  select(-Stage_raw, -Site_raw)

message("\nTotal valid participants after Stage/Site cleaning: ", nrow(df_all))

if (nrow(df_all) != 193) {
  warning(
    "Expected 193 valid participants, but found ",
    nrow(df_all),
    ". Check Stage and Site coding before interpreting results."
  )
}

message("\nSample size by Group x Stage:")
print(table(df_all$Group, df_all$Stage))

message("\nSample size by Group x Site:")
print(table(df_all$Group, df_all$Site))

message("\nSample size by Site x Group x Stage:")
print(xtabs(~ Site + Group + Stage, data = df_all))

# ------------------------------------------------------------
# Shared helper functions
# ------------------------------------------------------------

partial_eta2_from_f <- function(f_value, df_effect, df_error) {
  ifelse(
    is.na(f_value) | is.na(df_effect) | is.na(df_error),
    NA_real_,
    (f_value * df_effect) / (f_value * df_effect + df_error)
  )
}

clean_art_anova <- function(art_model, dv_name, valid_n) {
  as.data.frame(anova(art_model)) %>%
    rownames_to_column(var = "Effect") %>%
    mutate(
      DV = dv_name,
      Valid_N = valid_n,
      p_unc = `Pr(>F)`,
      partial_eta2 = partial_eta2_from_f(
        `F value`,
        Df,
        Df.res
      )
    ) %>%
    select(
      DV,
      Valid_N,
      Effect,
      Df,
      Df.res,
      `F value`,
      p_unc,
      partial_eta2
    )
}

# ============================================================
# Analysis A: Recruitment-site balance
# ============================================================

site_balance <- df_all %>%
  count(Group, Site, name = "N") %>%
  group_by(Group) %>%
  mutate(
    Group_Total = sum(N),
    Percent_within_Group = 100 * N / Group_Total
  ) %>%
  ungroup() %>%
  arrange(Group, Site)

site_table <- table(df_all$Group, df_all$Site)
site_chisq <- chisq.test(site_table, correct = FALSE)

cramers_v <- sqrt(
  as.numeric(site_chisq$statistic) /
    (
      sum(site_table) *
        min(nrow(site_table) - 1, ncol(site_table) - 1)
    )
)

site_balance_test <- tibble(
  Test = "Pearson chi-square test",
  Chi_square = as.numeric(site_chisq$statistic),
  Df = as.numeric(site_chisq$parameter),
  p_value = site_chisq$p.value,
  Cramers_V = cramers_v,
  Minimum_expected_count = min(site_chisq$expected),
  Total_N = sum(site_table)
)

site_expected_counts <- as.data.frame(
  as.table(site_chisq$expected),
  stringsAsFactors = FALSE,
  responseName = "Expected_N"
)

# Column names produced by as.data.frame.table() can depend on
# whether the input dimensions already have names. Assigning the
# three names by position is stable across R versions.
names(site_expected_counts) <- c("Group", "Site", "Expected_N")

site_expected_counts <- site_expected_counts %>%
  arrange(Group, Site)

message("\nAnalysis A: Site balance")
print(site_balance)
print(site_balance_test)

# ============================================================
# Analysis B1: RQ1 full-factorial ART including Site
# Model: Outcome ~ Group * Site
# ============================================================

fit_rq1_site_adjusted <- function(dv_name, data) {
  message("\nRQ1 site-adjusted ART: ", dv_name)
  
  dat <- data %>%
    select(Group, Site, DV_value = all_of(dv_name)) %>%
    mutate(DV_value = as.numeric(DV_value)) %>%
    drop_na()
  
  if (nrow(dat) == 0) {
    stop("No valid observations for outcome: ", dv_name)
  }
  
  if (nlevels(droplevels(dat$Group)) < 4) {
    stop("Fewer than four Group levels for outcome: ", dv_name)
  }
  
  if (nlevels(droplevels(dat$Site)) < 2) {
    stop("Fewer than two Site levels for outcome: ", dv_name)
  }
  
  art_model <- art(
    DV_value ~ Group * Site,
    data = dat
  )
  
  list(
    model = art_model,
    anova = clean_art_anova(art_model, dv_name, nrow(dat))
  )
}

rq1_fits <- setNames(
  lapply(dvs, fit_rq1_site_adjusted, data = df_all),
  dvs
)

rq1_site_adjusted_full <- bind_rows(
  lapply(rq1_fits, function(x) x$anova)
)

# The inferential family for RQ1 contains the nine adjusted
# omnibus Group tests, one for each outcome.
rq1_site_adjusted_group <- rq1_site_adjusted_full %>%
  filter(Effect == "Group") %>%
  mutate(
    p_adj = p.adjust(p_unc, method = "BH"),
    Significant = p_adj < alpha,
    DV_order = match(DV, dvs)
  ) %>%
  arrange(DV_order) %>%
  select(
    DV,
    Valid_N,
    Effect,
    Df,
    Df.res,
    `F value`,
    p_unc,
    p_adj,
    Significant,
    partial_eta2
  )

# For transparency, Site effects are saved separately. If they
# are interpreted inferentially, their nine p-values are treated
# as a separate family and adjusted using BH.
rq1_site_effects <- rq1_site_adjusted_full %>%
  filter(Effect == "Site") %>%
  mutate(
    p_adj = p.adjust(p_unc, method = "BH"),
    Significant = p_adj < alpha,
    DV_order = match(DV, dvs)
  ) %>%
  arrange(DV_order) %>%
  select(
    DV,
    Valid_N,
    Effect,
    Df,
    Df.res,
    `F value`,
    p_unc,
    p_adj,
    Significant,
    partial_eta2
  )

# Group-by-Site interactions are a separate nine-test family.
# A significant interaction indicates that the Group effect
# should not be summarized using marginal Group contrasts alone.
rq1_group_site_interactions <- rq1_site_adjusted_full %>%
  filter(Effect == "Group:Site") %>%
  mutate(
    p_adj = p.adjust(p_unc, method = "BH"),
    Significant = p_adj < alpha,
    DV_order = match(DV, dvs)
  ) %>%
  arrange(DV_order) %>%
  select(
    DV,
    Valid_N,
    Effect,
    Df,
    Df.res,
    `F value`,
    p_unc,
    p_adj,
    Significant,
    partial_eta2
  )

# Marginal Group contrasts are run only when the Group omnibus
# effect is FDR-significant and the Group-by-Site interaction is
# not FDR-significant. Significant interactions require separate
# simple-effects follow-up and should not be reduced to one
# marginal Group comparison.
significant_rq1_dvs <- rq1_site_adjusted_group %>%
  select(DV, Group_Significant = Significant) %>%
  left_join(
    rq1_group_site_interactions %>%
      select(DV, Group_Site_Significant = Significant),
    by = "DV"
  ) %>%
  filter(Group_Significant, !Group_Site_Significant) %>%
  pull(DV)

if (length(significant_rq1_dvs) > 0) {
  rq1_adjusted_contrasts <- bind_rows(
    lapply(significant_rq1_dvs, function(dv_name) {
      model <- rq1_fits[[dv_name]]$model
      
      group_lm <- artlm(model, "Group")
      group_emm <- emmeans(group_lm, ~ Group)
      
      as.data.frame(
        contrast(
          group_emm,
          method = "trt.vs.ctrl",
          ref = 1,
          adjust = "BH"
        )
      ) %>%
        mutate(
          DV = dv_name,
          Reference_Group = "G1",
          Adjustment = "BH within the three G2/G3/G4 vs G1 contrasts"
        ) %>%
        select(
          DV,
          Reference_Group,
          contrast,
          estimate,
          SE,
          df,
          t.ratio,
          p.value,
          Adjustment
        )
    })
  )
} else {
  rq1_adjusted_contrasts <- tibble(
    DV = character(),
    Reference_Group = character(),
    contrast = character(),
    estimate = numeric(),
    SE = numeric(),
    df = numeric(),
    t.ratio = numeric(),
    p.value = numeric(),
    Adjustment = character()
  )
}

message("\nAnalysis B1: RQ1 site-adjusted Group effects")
print(
  rq1_site_adjusted_group %>%
    mutate(
      `F value` = round(`F value`, 3),
      p_unc = round(p_unc, 3),
      p_adj = round(p_adj, 3),
      partial_eta2 = round(partial_eta2, 3)
    )
)

# ============================================================
# Analysis B2: RQ2 full-factorial ART including Site
# Model: Outcome ~ Group * Stage * Site
# The full factorial is required by ARTool when Site is included
# as a fixed factor.
# ============================================================

fit_rq2_site_adjusted <- function(dv_name, data) {
  message("\nRQ2 site-adjusted ART: ", dv_name)
  
  dat <- data %>%
    select(
      Group,
      Stage,
      Site,
      DV_value = all_of(dv_name)
    ) %>%
    mutate(DV_value = as.numeric(DV_value)) %>%
    drop_na()
  
  if (nrow(dat) == 0) {
    stop("No valid observations for outcome: ", dv_name)
  }
  
  if (nlevels(droplevels(dat$Group)) < 4) {
    stop("Fewer than four Group levels for outcome: ", dv_name)
  }
  
  if (nlevels(droplevels(dat$Stage)) < 2) {
    stop("Fewer than two Stage levels for outcome: ", dv_name)
  }
  
  if (nlevels(droplevels(dat$Site)) < 2) {
    stop("Fewer than two Site levels for outcome: ", dv_name)
  }
  
  art_model <- art(
    DV_value ~ Group * Stage * Site,
    data = dat
  )
  
  list(
    model = art_model,
    anova = clean_art_anova(art_model, dv_name, nrow(dat))
  )
}

rq2_fits <- setNames(
  lapply(dvs, fit_rq2_site_adjusted, data = df_all),
  dvs
)

rq2_site_adjusted_full <- bind_rows(
  lapply(rq2_fits, function(x) x$anova)
)

rq2_effect_order <- c("Group", "Stage", "Group:Stage")
rq2_effect_labels <- c(
  Group = "AI presentation condition",
  Stage = "Training stage",
  `Group:Stage` = "Condition x Stage"
)

# The RQ2 family remains consistent with the original script:
# 9 outcomes x 3 target effects = 27 omnibus p-values.
# The nine nuisance Site effects are not included in this family.
rq2_site_adjusted_target <- rq2_site_adjusted_full %>%
  filter(Effect %in% rq2_effect_order) %>%
  mutate(
    Effect_Label = unname(rq2_effect_labels[Effect]),
    p_adj = p.adjust(p_unc, method = "BH"),
    Significant = p_adj < alpha,
    DV_order = match(DV, dvs),
    Effect_order = match(Effect, rq2_effect_order)
  ) %>%
  arrange(DV_order, Effect_order) %>%
  select(
    DV,
    Valid_N,
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

# Saved for transparency as a separate nine-test family.
rq2_site_effects <- rq2_site_adjusted_full %>%
  filter(Effect == "Site") %>%
  mutate(
    p_adj = p.adjust(p_unc, method = "BH"),
    Significant = p_adj < alpha,
    DV_order = match(DV, dvs)
  ) %>%
  arrange(DV_order) %>%
  select(
    DV,
    Valid_N,
    Effect,
    Df,
    Df.res,
    `F value`,
    p_unc,
    p_adj,
    Significant,
    partial_eta2
  )

# Site-related interactions form a separate family:
# 9 outcomes x 3 interactions = 27 tests.
rq2_site_interaction_order <- c(
  "Group:Site",
  "Stage:Site",
  "Group:Stage:Site"
)

rq2_site_interaction_labels <- c(
  `Group:Site` = "Condition x Site",
  `Stage:Site` = "Training stage x Site",
  `Group:Stage:Site` = "Condition x Training stage x Site"
)

rq2_site_interactions <- rq2_site_adjusted_full %>%
  filter(Effect %in% rq2_site_interaction_order) %>%
  mutate(
    Effect_Label = unname(rq2_site_interaction_labels[Effect]),
    p_adj = p.adjust(p_unc, method = "BH"),
    Significant = p_adj < alpha,
    DV_order = match(DV, dvs),
    Effect_order = match(Effect, rq2_site_interaction_order)
  ) %>%
  arrange(DV_order, Effect_order) %>%
  select(
    DV,
    Valid_N,
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

message("\nAnalysis B2: RQ2 site-adjusted target effects")
print(
  rq2_site_adjusted_target %>%
    mutate(
      `F value` = round(`F value`, 3),
      p_unc = round(p_unc, 3),
      p_adj = round(p_adj, 3),
      partial_eta2 = round(partial_eta2, 3)
    )
)

# ============================================================
# Analysis C: Descriptive statistics by recruitment site
# These values are unadjusted and are not hypothesis tests.
# ============================================================

site_descriptives <- df_all %>%
  select(Site, all_of(dvs)) %>%
  pivot_longer(
    cols = all_of(dvs),
    names_to = "DV",
    values_to = "Score"
  ) %>%
  mutate(
    Score = as.numeric(Score),
    DV_order = match(DV, dvs)
  ) %>%
  group_by(DV_order, DV, Site) %>%
  summarise(
    N = sum(!is.na(Score)),
    Missing_N = sum(is.na(Score)),
    Mean = if (sum(!is.na(Score)) > 0) mean(Score, na.rm = TRUE) else NA_real_,
    SD = if (sum(!is.na(Score)) > 1) sd(Score, na.rm = TRUE) else NA_real_,
    Median = if (sum(!is.na(Score)) > 0) median(Score, na.rm = TRUE) else NA_real_,
    Q1 = if (sum(!is.na(Score)) > 0) as.numeric(quantile(Score, 0.25, na.rm = TRUE, names = FALSE)) else NA_real_,
    Q3 = if (sum(!is.na(Score)) > 0) as.numeric(quantile(Score, 0.75, na.rm = TRUE, names = FALSE)) else NA_real_,
    IQR = if (sum(!is.na(Score)) > 0) IQR(Score, na.rm = TRUE) else NA_real_,
    Minimum = if (sum(!is.na(Score)) > 0) min(Score, na.rm = TRUE) else NA_real_,
    Maximum = if (sum(!is.na(Score)) > 0) max(Score, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(DV_order, Site) %>%
  select(-DV_order)
group_site_descriptives_long <- df_all %>%
  select(Group, Site, all_of(dvs)) %>%
  pivot_longer(
    cols = all_of(dvs),
    names_to = "DV",
    values_to = "Score"
  ) %>%
  mutate(
    Score = as.numeric(Score),
    DV_order = match(DV, dvs)
  ) %>%
  group_by(DV_order, DV, Group, Site) %>%
  summarise(
    N = sum(!is.na(Score)),
    Missing_N = sum(is.na(Score)),
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    Median = median(Score, na.rm = TRUE),
    Q1 = as.numeric(
      quantile(Score, 0.25, na.rm = TRUE, names = FALSE)
    ),
    Q3 = as.numeric(
      quantile(Score, 0.75, na.rm = TRUE, names = FALSE)
    ),
    IQR = IQR(Score, na.rm = TRUE),
    Minimum = min(Score, na.rm = TRUE),
    Maximum = max(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(DV_order, Group, Site)
message("\nAnalysis C: Descriptive statistics by Site")
print(site_descriptives)

group_site_descriptives <- group_site_descriptives_long %>%
  select(
    DV_order,
    DV,
    Group,
    Site,
    N,
    Mean,
    SD,
    Median,
    Q1,
    Q3,
    IQR
  ) %>%
  pivot_wider(
    names_from = Site,
    values_from = c(
      N,
      Mean,
      SD,
      Median,
      Q1,
      Q3,
      IQR
    ),
    names_glue = "{.value}_{Site}"
  ) %>%
  mutate(
    Mean_Difference_China_minus_Singapore =
      Mean_China - Mean_Singapore,
    
    Median_Difference_China_minus_Singapore =
      Median_China - Median_Singapore
  ) %>%
  arrange(DV_order, Group) %>%
  select(-DV_order)

# ------------------------------------------------------------
# Sample-size tables
# ------------------------------------------------------------

sample_size_group_stage <- df_all %>%
  count(Group, Stage, name = "N") %>%
  pivot_wider(
    names_from = Stage,
    values_from = N,
    values_fill = 0
  )

sample_size_site_group_stage <- df_all %>%
  count(Site, Group, Stage, name = "N") %>%
  arrange(Site, Group, Stage)

# ------------------------------------------------------------
# Analysis notes
# ------------------------------------------------------------

notes <- tibble(
  Item = c(
    "Recruitment-site definition",
    "Analysis A",
    "Analysis B1 model",
    "Analysis B1 correction",
    "Analysis B1 contrasts",
    "Analysis B2 model",
    "Analysis B2 correction",
    "Site effects",
    "Site interactions",
    "Effect size",
    "Analysis C",
    "Significance criterion"
  ),
  Description = c(
    paste0(
      "Site denotes the country of the medical school through which a participant was recruited; ",
      "it does not denote citizenship, nationality, ethnicity, or individual cultural identity."
    ),
    paste0(
      "Group-by-Site balance was assessed using a Pearson chi-square test. ",
      "Counts, within-Group percentages, minimum expected count, and Cramer's V are reported."
    ),
    paste0(
      "For each outcome, RQ1 sensitivity analysis used the ARTool-required full factorial: ",
      "ARTool::art(DV ~ Group * Site)."
    ),
    paste0(
      "Benjamini-Hochberg FDR correction was applied separately across the nine omnibus Group tests ",
      "and across the nine Group x Site interaction tests."
    ),
    paste0(
      "For outcomes with an FDR-significant Group omnibus test and no FDR-significant Group x Site interaction, ",
      "G2/G3/G4 were compared with G1 using emmeans contrasts on ARTool::artlm(); ",
      "BH correction was applied within the three contrasts for each outcome."
    ),
    paste0(
      "For each outcome, RQ2 sensitivity analysis used ",
      "the ARTool-required full factorial: ARTool::art(DV ~ Group * Stage * Site)."
    ),
    paste0(
      "Benjamini-Hochberg FDR correction was applied across the original RQ2 target family: ",
      "nine outcomes x three effects (Group, Stage, Group x Stage) = 27 tests."
    ),
    paste0(
      "Site effects are saved separately. If interpreted inferentially, ",
      "their nine p-values are treated as a separate BH-adjusted family."
    ),
    paste0(
      "RQ1 Group x Site interactions were BH-adjusted across nine outcomes. ",
      "RQ2 site-related interactions (Group x Site, Stage x Site, and Group x Stage x Site) ",
      "were treated as a separate 27-test BH-adjusted family."
    ),
    paste0(
      "Partial eta squared was calculated from the ART ANOVA F statistic as ",
      "F*df1/(F*df1+df2) and should be interpreted as an effect size on the aligned-rank analysis scale."
    ),
    paste0(
      "Site-stratified outcome values are descriptive and unadjusted. ",
      "No Singapore-versus-China hypothesis tests were conducted in Analysis C."
    ),
    paste0("FDR-adjusted p < ", alpha)
  )
)

# ------------------------------------------------------------
# Export all outputs
# ------------------------------------------------------------

write_xlsx(
  list(
    Site_Balance = site_balance,
    Site_Balance_Test = site_balance_test,
    Site_Expected_Counts = site_expected_counts,
    RQ1_Site_Adjusted_ART = rq1_site_adjusted_group,
    RQ1_Group_Site_Int = rq1_group_site_interactions,
    RQ1_Adjusted_Contrasts = rq1_adjusted_contrasts,
    RQ1_Site_Effects = rq1_site_effects,
    RQ1_ART_Full = rq1_site_adjusted_full,
    RQ2_Site_Adjusted_ART = rq2_site_adjusted_target,
    RQ2_Site_Interactions = rq2_site_interactions,
    RQ2_Site_Effects = rq2_site_effects,
    RQ2_ART_Full = rq2_site_adjusted_full,
    Site_Descriptives = site_descriptives,
    N_Group_Stage = sample_size_group_stage,
    N_Site_Group_Stage = sample_size_site_group_stage,
    Group_Site_Desc = group_site_descriptives,
    Group_Site_Desc_Long = group_site_descriptives_long,
    Notes = notes
  ),
  path = output_path
)

message("\nDone. Results saved to: ", output_path)
