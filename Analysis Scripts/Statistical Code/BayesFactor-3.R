# ==============================================================================
# Bayes factor analysis with prior sensitivity analysis
# ==============================================================================

library(readxl)
library(dplyr)
library(BayesFactor)
library(openxlsx)

# ------------------------------------------------------------------------------
# 1. Data loading
# ------------------------------------------------------------------------------
file_path <- "DataA.xlsx"
groups <- c("G1", "G2", "G3", "G4")

df_list <- list()

for (g in groups) {
  temp <- read_excel(file_path, sheet = g)
  temp$Group <- g
  df_list[[length(df_list) + 1]] <- temp
}

df <- bind_rows(df_list)
df$Group <- as.factor(df$Group)

# ------------------------------------------------------------------------------
# 2. Variables
# ------------------------------------------------------------------------------
target_cols <- c(
  "Positive Affect Score",
  "Negative Affect Score",
  "Technical Competence Score",
  "Understandability Score",
  "Faith Score",
  "Attachment Score",
  "Overall HCT Score",
  "Usability Score"
)

# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

classify_bf <- function(bf10) {
  if (is.na(bf10)) {
    return("No Data / Failed")
  } else if (bf10 > 3) {
    return("Evidence favoring H1")
  } else if (bf10 < 1/3) {
    return("Evidence favoring H0")
  } else {
    return("Data insensitive / inconclusive")
  }
}

calculate_bayes <- function(data, var_name, prior_scale) {
  
  clean_data <- data %>%
    select(Group, all_of(var_name)) %>%
    na.omit()
  
  if (nrow(clean_data) == 0) {
    return(data.frame(
      Variable = var_name,
      Prior = prior_scale,
      N = 0,
      BF10 = NA,
      BF01 = NA,
      Conclusion = "No Data / Failed"
    ))
  }
  
  clean_data$Ranked_Value <- rank(clean_data[[var_name]])
  
  bf_obj <- anovaBF(
    Ranked_Value ~ Group,
    data = clean_data,
    rscaleFixed = prior_scale,
    progress = FALSE
  )
  
  bf10 <- as.numeric(as.vector(bf_obj)[1])
  bf01 <- 1 / bf10
  
  data.frame(
    Variable = var_name,
    Prior = prior_scale,
    N = nrow(clean_data),
    BF10 = bf10,
    BF01 = bf01,
    Conclusion = classify_bf(bf10)
  )
}

# ------------------------------------------------------------------------------
# 4. Primary analysis and sensitivity analysis
# ------------------------------------------------------------------------------

# Primary prior: default medium fixed-effect prior in BayesFactor
# Sensitivity priors: wider fixed-effect priors
prior_scales <- c("medium", "wide", "ultrawide")

results_all <- bind_rows(lapply(target_cols, function(var) {
  bind_rows(lapply(prior_scales, function(prior_scale) {
    calculate_bayes(df, var, prior_scale)
  }))
}))

primary_results <- results_all %>%
  filter(Prior == "medium") %>%
  mutate(
    BF10 = round(BF10, 3),
    BF01 = round(BF01, 3)
  )

sensitivity_results <- results_all %>%
  mutate(
    BF10 = round(BF10, 3),
    BF01 = round(BF01, 3)
  )

sensitivity_summary <- results_all %>%
  group_by(Variable) %>%
  summarise(
    N = first(N),
    Medium_Conclusion = Conclusion[Prior == "medium"],
    Wide_Conclusion = Conclusion[Prior == "wide"],
    Ultrawide_Conclusion = Conclusion[Prior == "ultrawide"],
    Sensitivity_Consistent = length(unique(Conclusion)) == 1,
    .groups = "drop"
  )

notes <- data.frame(
  Item = c(
    "Primary prior",
    "Sensitivity analysis",
    "Outcome transformation",
    "BF10",
    "BF01",
    "Interpretation"
  ),
  Description = c(
    "The primary analysis used rscaleFixed = 'medium' for the fixed Group effect.",
    "Prior sensitivity was assessed by repeating analyses with rscaleFixed = 'wide' and 'ultrawide'.",
    "Outcome scores were ranked before Bayesian one-way ANOVA.",
    "BF10 quantifies evidence for the Group model relative to the null model.",
    "BF01 was calculated as 1 / BF10.",
    "BF10 > 3 was interpreted as evidence favoring H1; BF10 < 1/3 as evidence favoring H0; intermediate values as inconclusive."
  )
)

# ------------------------------------------------------------------------------
# 5. Output
# ------------------------------------------------------------------------------

print("========================================================")
print("Primary Bayes factor results: rscaleFixed = 'medium'")
print("========================================================")
print(primary_results, row.names = FALSE)

print("========================================================")
print("Prior sensitivity summary")
print("========================================================")
print(sensitivity_summary, row.names = FALSE)

write.xlsx(
  list(
    Primary_medium_prior = primary_results,
    Sensitivity_all_priors = sensitivity_results,
    Sensitivity_summary = sensitivity_summary,
    Notes = notes
  ),
  file = "Bayes_Results_with_Prior_Sensitivity.xlsx",
  overwrite = TRUE
)