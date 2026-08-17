# This script uses the exposure num_days instead of MET_mins to look at the main analysis. 
# To better understand effect of the missing physical activity data if any

source(here::here("scripts/00_setup.R"))
source(here("scripts/03_helpers_general.R"))
source(here("scripts/04_helpers_CXA.R"))
source(here("scripts/05_helpers_survival.R"))
source(here("scripts/06_helpers_figures_CXA.R"))


# Load processed data from previous scripts
sensitivity_dat  <- readRDS(file.path(DATA_DERIVED, "analysis_dat.Rds"))

# Keep all variables but rename by dropping the words 'raw' or 'clean'
sensitivity_dat <- sensitivity_dat %>%
  rename_with(~ gsub("_clean|_raw", "", .x))

# rename variables for future analysis
sensitivity_dat <- sensitivity_dat %>%
  rename(
    SRF = self_reported_fracture_A0,
    WC = waist_circ
  )

# have had issues with date columns so run this first

sensitivity_dat$age_A0 <- as.numeric(sensitivity_dat$age_A0)
sensitivity_dat <- sensitivity_dat %>%
  dplyr::mutate(
    date_assess_A0 = data.table::as.IDate(as.Date(date_assess_A0)),
    lost_to_fu     = data.table::as.IDate(as.Date(lost_to_fu))
  )


# Create middle aged cohort 

sensitivity_dat <- sensitivity_dat %>%
  filter(between(age_A0, 40, 65))

# Create new variable called num_days_total
sensitivity_dat <- sensitivity_dat %>%
  mutate(
    num_days_total = if_else(
      !is.na(num_day_walk) &
        !is.na(num_days_mod) &
        !is.na(num_days_vig),
      num_day_walk + num_days_mod + num_days_vig,
      NA_real_
    )
  )

# Analysis is restricted to 40-65 (in previous script) but age bands are in 10 years
# Check 60-69 should not have anyone over 65, relabel as 60-65 for subsequent figures/tables

# Check the range of ages in the 60-69 group
sensitivity_dat %>%
  filter(agegp_A0 == "60-69") %>%
  summarise(
    min_age = min(age_A0, na.rm = TRUE),
    max_age = max(age_A0, na.rm = TRUE),
    n = n()
  )

sensitivity_dat <- sensitivity_dat %>%
  mutate(agegp_A0 = forcats::fct_recode(agegp_A0,
                                        "60-65" = "60-69"))


## Create a complete case analysis data set but do not exclude for PA

## Exclusions 

# Before working with the data we want to exclude those with missing data in covariates and pre-existing fracture
# We will record how many participants are excluded at each of the steps (e.g. for a flow diagram):

tab_exc <- data.frame("Exclusion" = "Starting cohort (age 40-65)",
                      "Number_excluded" = NA, 
                      "Number_remaining" = nrow(sensitivity_dat)
)


# First, we exclude participants with missing fracture data

# Store number before exclusion
nb <- nrow(sensitivity_dat)

# Exclude rows with NA in SRF
sensitivity_dat <- sensitivity_dat[!is.na(sensitivity_dat$SRF), ]

# Add to exclusion table
tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing fracture data",
    Number_excluded = nb - nrow(sensitivity_dat),
    Number_remaining = nrow(sensitivity_dat)
  )
)


# ---- Missing ethnicity ----
nb <- nrow(sensitivity_dat)

sensitivity_dat <- sensitivity_dat[!is.na(sensitivity_dat$ethnicity_derived), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing Ethnicity data",
    Number_excluded = nb - nrow(sensitivity_dat),
    Number_remaining = nrow(sensitivity_dat)
  )
)

# ---- Missing deprivation ----
nb <- nrow(sensitivity_dat)

sensitivity_dat <- sensitivity_dat[!is.na(sensitivity_dat$tdi), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing deprivation data",
    Number_excluded = nb - nrow(sensitivity_dat),
    Number_remaining = nrow(sensitivity_dat)
  )
)

# ---- Missing education ----
nb <- nrow(sensitivity_dat)

sensitivity_dat <- sensitivity_dat[!is.na(sensitivity_dat$education_level), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing qualification data",
    Number_excluded = nb - nrow(sensitivity_dat),
    Number_remaining = nrow(sensitivity_dat)
  )
)

# ---- Missing anthropometrics ----
nb <- nrow(sensitivity_dat)

sensitivity_dat <- sensitivity_dat[
  !is.na(sensitivity_dat$weight) &
    !is.na(sensitivity_dat$height) &
    !is.na(sensitivity_dat$WC) &
    !is.na(sensitivity_dat$BMI),
]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing anthropometric data",
    Number_excluded = nb - nrow(sensitivity_dat),
    Number_remaining = nrow(sensitivity_dat)
  )
)

# ---- Missing PA  data ----
nb <- nrow(sensitivity_dat)

sensitivity_dat <- sensitivity_dat[
  !is.na(sensitivity_dat$num_days_total),
]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing activity data",
    Number_excluded = nb - nrow(sensitivity_dat),
    Number_remaining = nrow(sensitivity_dat)
  )
)

View(tab_exc)

### Now PA exclusions are only 8% 

# Build the PNG file in FIGURES_DIR
png(
  filename = file.path(FIGURES_DIR, "sensitivity_Exclusion_flow.png"),
  width = 1200,
  height = 800,
  res = 150
)

# Render the table as a figure
gridExtra::grid.table(tab_exc)

# Close the device
dev.off()

### Conduct main analysis of quintiles using this dataset in order to compare magnitude and direction of ORs

# =========================================================
# Load and prepare dataset
# =========================================================


sensitivity_dat <- sensitivity_dat %>%
  dplyr::select(
    eid,
    ethnicity = ethnicity_derived,
    date_assess_A0,
    assessment_country,
    age_A0,
    agegp_A0,
    age_over_50,
    sex,
    education_level,
    tdi,
    WC,
    weight,
    height,
    BMI,
    lost_to_fu,
    approx_dob_derived,
    SRF,
    Wrist,
    Hip,
   num_days_total,
  ) %>%
  dplyr::mutate(
    SRF = factor(SRF, levels = c("No", "Yes")),
    Hip = factor(Hip, levels = c(0, 1), labels = c("No", "Yes")),
    Wrist = factor(Wrist, levels = c(0, 1), labels = c("No", "Yes"))
  )


# =========================================================
# Analysis settings
# =========================================================

n_bins <- 5

# Primary quintile analysis
quintile_exposures_primary <- c("num_days_total")
quintile_outcomes_primary <- c("SRF")

# Age-sex combined stratified analysis
quintile_exposures_age_sex <- c("num_days_total")
quintile_outcomes_age_sex <- c("SRF")

model_specs_sex <- list(
  "Unadjusted" = character(0),
  "Minimal adjustment" = c("age_A0", "weight", "height"),
  "Fully adjusted" = c("age_A0", "ethnicity", "weight", "height", "tdi", "education_level")
)

model_specs_age <- list(
  "Unadjusted" = character(0),
  "Minimal adjustment" = c("sex", "weight", "height"),
  "Fully adjusted" = c("sex", "ethnicity", "weight", "height", "tdi", "education_level")
)

model_specs_age_sex <- list(
  "Unadjusted" = character(0),
  "Minimal adjustment" = c("weight", "height"),
  "Fully adjusted" = c("ethnicity", "weight", "height", "tdi", "education_level")
)

model_levels <- c("Unadjusted", "Minimal adjustment", "Fully adjusted")

# =========================================================
# 3. Quintile analysis stratified by sex and age group
# =========================================================

sex_levels <- unique(as.character(sensitivity_dat$sex))
age_levels <- setdiff(unique(as.character(sensitivity_dat$agegp_A0)), "70-79")

strata_age_sex <- expand.grid(
  sex = sex_levels,
  agegp_A0 = age_levels,
  stringsAsFactors = FALSE
)

sensitivity_quintiles_age_sex <- run_quintile_models_stratified(
  data = sensitivity_dat,
  strata_df = strata_age_sex,
  outcomes = quintile_outcomes_age_sex,
  exposures = quintile_exposures_age_sex,
  model_specs = model_specs_age_sex,
  n_bins = n_bins
) %>%
  format_quintile_results(model_levels = model_levels, digits = 5)

saveRDS(
  sensitivity_quintiles_age_sex,
  file.path(DATA_DERIVED, "sensitivity_quintiles_age_sex.rds")
)

save_table_word(
  df = sensitivity_quintiles_age_sex,
  table_number = "Table XSensitivity",
  folder_path = TABLES_DIR,
  title = "Sensitivity Odds ratios for raw PA quintiles across three models for SRF, stratified by sex and age group"
)


# =========================================================
# Sensitivity analysis: alternative exposures
# MET_mod, MET_vig, MET_walk
# Outcome = SRF
# Stratified by sex
# =========================================================

results_quintiles_sex_exposure_sens <- readRDS(
  file.path(DATA_DERIVED, "results_quintiles_sex_exposure_sens.rds")
)

save_table_word(
  df = results_quintiles_sex_exposure_sens,
  table_number = "X_OR_Mod_Vig_walk_sex",
  folder_path = TABLES_DIR,
  title = "Odds ratios for raw MET_mod, MET_vig and MET_walk quintiles across three models for self-reported fracture, stratified by sex"
)

plot_df_sex_exposure_sens <- prep_quintile_plot_df(
  df = results_quintiles_sex_exposure_sens,
  facet_vars = c("sex"),
  exposure_filter = c("MET_mod", "MET_vig", "MET_walk"),
  outcome_filter = "SRF"
) %>%
  dplyr::mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    model = factor(model, levels = model_levels),
    exposure = factor(
      exposure,
      levels = c("MET_mod", "MET_vig", "MET_walk"),
      labels = c("Moderate", "Vigorous", "Walking")
    )
  )

p_exposure_sens <- plot_quintile_or(
  data = plot_df_sex_exposure_sens,
  facet_formula = exposure ~ sex,
  title_text = "Odds ratios for MET_mod, MET_vig and MET_walk quintiles and self-reported fracture by sex"
)

p_exposure_sens

save_quintile_plot(
  plot_obj = p_exposure_sens,
  filename = "Fig_X3_7_SRF_MET_components_sex.png",
  figures_dir = FIGURES_DIR
)

# =========================================================
# Sensitivity figures: alternative PA exposures by sex
# =========================================================

results_sex_only <- readRDS(
  file.path(DATA_DERIVED, "results_sensitivity_quintiles_sex.rds")
)

results_plot_df <- results_sex_only %>%
  dplyr::mutate(
    Quintile = sub(".*(Q[1-5])$", "\\1", term),
    exposure = factor(
      exposure,
      levels = c("MET_total", "MET_MVPA", "MET_mod", "MET_vig", "MET_walk"),
      labels = c("Total", "MVPA", "Moderate", "Vigorous", "Walking")
    ),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>%
  dplyr::filter(Quintile %in% paste0("Q", 2:5))

ref_rows <- results_plot_df %>%
  dplyr::distinct(sex, outcome, exposure) %>%
  dplyr::mutate(
    term = NA_character_,
    Quintile = "Q1",
    OR = 1,
    CI_lower = 1,
    CI_upper = 1,
    p_value = NA_real_
  )

results_plot_df <- dplyr::bind_rows(results_plot_df, ref_rows) %>%
  dplyr::mutate(
    Quintile = factor(Quintile, levels = rev(paste0("Q", 1:5)))
  )

p_female <- plot_quintile_forest(
  data = dplyr::filter(results_plot_df, sex == "Female"),
  title_text = "SRF: fully adjusted odds ratios across PA quintiles in females"
)

p_male <- plot_quintile_forest(
  data = dplyr::filter(results_plot_df, sex == "Male"),
  title_text = "SRF: fully adjusted odds ratios across PA quintiles in males"
)

save_quintile_plot(p_female, "Fig_X3_5_SRF_PA_Quintiles_Female.png", FIGURES_DIR)
save_quintile_plot(p_male, "Fig_X3_6_SRF_PA_Quintiles_Male.png", FIGURES_DIR)



# =========================================================
# MET_total distribution by quintile, sex-stratified
# =========================================================

IPAQ_dat2 <- readRDS(file.path(DATA_DERIVED, "IPAQ_dat2.rds"))

quintile_MET_summary <- IPAQ_dat2 %>%
  dplyr::filter(sex %in% c("Male", "Female")) %>%
  dplyr::group_by(sex, MET_total_q) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_MET_total = round(mean(MET_total, na.rm = TRUE), 0),
    median_MET_total = round(median(MET_total, na.rm = TRUE), 0),
    p25_MET_total = round(stats::quantile(MET_total, 0.25, na.rm = TRUE), 0),
    p75_MET_total = round(stats::quantile(MET_total, 0.75, na.rm = TRUE), 0),
    min_MET_total = round(min(MET_total, na.rm = TRUE), 0),
    max_MET_total = round(max(MET_total, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    IQR_MET_total = paste0(p25_MET_total, "–", p75_MET_total),
    threshold_MET = round(exp(7.5) - 1, 0),
    threshold_in_IQR = ifelse(
      threshold_MET >= p25_MET_total & threshold_MET <= p75_MET_total,
      "Yes",
      "No"
    )
  ) %>%
  dplyr::arrange(sex, MET_total_q)

write.csv(
  quintile_MET_summary,
  file.path(TABLES_DIR, "MET_total_quintile_summary_with_threshold.csv"),
  row.names = FALSE
)


