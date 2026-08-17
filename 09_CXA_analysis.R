# =========================================================
# Cross-sectional PA analysis
# =========================================================

library(here)
source(here::here("scripts/00_setup.R"))
source(here("scripts/03_helpers_general.R"))
source(here("scripts/04_helpers_CXA.R"))
source(here("scripts/05_helpers_survival.R"))

# =========================================================
# Load and prepare dataset
# =========================================================

IPAQ_dat2 <- readRDS(file.path(DATA_DERIVED, "complete_case_dat.Rds"))

IPAQ_dat2 <- IPAQ_dat2 %>%
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
    Leg,
    Arm,
    Ankle,
    Spine,
    Other =`Other bones`,
    duration_mod,
    num_days_mod,
    duration_vig,
    num_days_vig,
    num_day_walk,
    duration_walk,
    IPAQ,
    mins_wk_mod = cc_mins_wk_mod_trunc,
    mins_wk_vig = cc_mins_wk_vig_trunc,
    mins_wk_MVPA = cc_mins_wk_MVPA_trunc,
    mins_wk_Total = cc_mins_wk_total_trunc,
    mins_wk_walk = cc_mins_wk_walk_trunc,
    MET_mod = cc_MET_mod_trunc,
    MET_vig = cc_MET_vig_trunc,
    MET_MVPA = cc_MET_MVPA_trunc,
    MET_total = cc_MET_total_trunc,
    MET_walk = cc_MET_walk_trunc,
    MET_mod_bin = any_mod,
    MET_vig_bin = any_vig,
    MET_MVPA_bin = any_MVPA,
    who_equiv_mins = WHO_MVPA_equiv_cont,
    WHO_MVPA_equiv_cat,
    risk_binary,
    drive_fast_binary,
    falls_binary,
    low_BMD,
    bmi_reds_binary,
    smoking_binary
  ) %>%
  dplyr::mutate(
    SRF = factor(SRF, levels = c("No", "Yes")),
    Hip = factor(Hip, levels = c(0, 1), labels = c("No", "Yes")),
    Wrist = factor(Wrist, levels = c(0, 1), labels = c("No", "Yes")),
    Leg = factor(Leg , levels = c(0, 1), labels = c("No", "Yes")),
    Arm = factor(Arm, levels = c(0, 1), labels = c("No", "Yes")),
    Ankle = factor(Ankle, levels = c(0, 1), labels = c("No", "Yes")),
    Spine = factor(Spine, levels = c(0, 1), labels = c("No", "Yes")),
    Other = factor(Other, levels = c(0, 1), labels = c("No", "Yes"))
  )

# Create log-transformed continuous exposure variables once
PA_cont_vars <- c("MET_MVPA", "MET_mod", "MET_vig", "MET_total", "MET_walk")

IPAQ_dat2 <- IPAQ_dat2 %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(PA_cont_vars),
      ~ log(.x + 1),
      .names = "{.col}_log"
    )
  )

# Create PA quintiles once for downstream tables and figures
for(pa in PA_cont_vars){
  
  quint_var <- paste0(pa, "_q")
  
  IPAQ_dat2[[quint_var]] <- dplyr::ntile(
    IPAQ_dat2[[paste0(pa, "_log")]],
    5
  )
  
  IPAQ_dat2[[quint_var]] <- factor(
    IPAQ_dat2[[quint_var]],
    levels = 1:5,
    labels = paste0("Q", 1:5)
  )
}

# Analysis is restricted to 40-65 (in previous script) but age bands are in 10 years
# Check 60-69 should not have anyone over 65, relabel as 60-65 for subsequent figures/tables

# Check the range of ages in the 60-69 group
IPAQ_dat2 %>%
  filter(agegp_A0 == "60-69") %>%
  summarise(
    min_age = min(age_A0, na.rm = TRUE),
    max_age = max(age_A0, na.rm = TRUE),
    n = n()
  )

IPAQ_dat2 <- IPAQ_dat2 %>%
  mutate(agegp_A0 = forcats::fct_recode(agegp_A0,
                                        "60-65" = "60-69"))

# =========================================================
# Analysis settings
# =========================================================

# Primary analysis outcome(s)
fracture_types_primary <- c("SRF")

# Outcomes for linearity testing
fracture_types_all <- c("SRF", "Hip", "Wrist")

# Incremental adjustment models
model_specs <- list(
  "Unadjusted" = NULL,
  "Age + sex adjusted" = c("age_A0", "sex"),
  "Age + sex + ethnicity adjusted" = c("age_A0", "sex", "ethnicity"),
  "Age + sex + ethnicity + anthropometry adjusted" = c("age_A0", "sex", "ethnicity", "weight", "height"),
  "Fully adjusted" = c("age_A0", "sex", "ethnicity", "weight", "height", "tdi", "education_level")
)

# Full adjustment set for linearity checks and spline plots
full_adjustment <- c("age_A0", "sex", "ethnicity", "weight", "height", "tdi", "education_level")

saveRDS(
  IPAQ_dat2,
  file.path(DATA_DERIVED, "IPAQ_dat2.rds")
)


# =========================================================
# Run primary continuous models
# =========================================================

OR_all <- run_continuous_models(
  data = IPAQ_dat2,
  outcomes = fracture_types_primary,
  exposures = PA_cont_vars,
  model_specs = model_specs
)


# Save model results for later tables/figures
saveRDS(
  OR_all,
  file.path(DATA_DERIVED, "OR_all_continuous.rds")
)


# =========================================================
# Assess linearity
# =========================================================

linearity_results <- test_linearity(
  data = IPAQ_dat2,
  outcomes = fracture_types_all,
  exposures = PA_cont_vars,
  confounders = full_adjustment,
  spline_knots = 4
) %>%
  dplyr::mutate(
    LRT_p_value = dplyr::if_else(
      LRT_p_value < 0.001,
      "<0.001",
      as.character(round(LRT_p_value, 4))
    )
  )

# Save results
saveRDS(
  linearity_results,
  file.path(DATA_DERIVED, "linearity_results.rds")
)


# =========================================================
# Spline plots
# =========================================================

save_spline_plots(
  data = IPAQ_dat2,
  outcome = "SRF",
  exposures = PA_cont_vars,
  confounders = full_adjustment,
  figures_dir = FIGURES_DIR,
  spline_knots = 4
)

# =========================================================
# Interaction testing: sex and age
# =========================================================

interaction_results <- dplyr::bind_rows(
  lapply(
    PA_cont_vars,
    function(x) {
      test_spline_interactions(
        data = IPAQ_dat2,
        outcome = "SRF",
        exposure = x,
        confounders = full_adjustment,
        spline_knots = 4,
        age_group_var = "agegp_A0"
      )
    }
  )
) %>%
  dplyr::mutate(
    LRT_p_value = dplyr::if_else(
      LRT_p_value < 0.001,
      "<0.001",
      as.character(round(LRT_p_value, 4))
    )
  )

saveRDS(
  interaction_results,
  file.path(DATA_DERIVED, "interaction_results.rds")
)


# =========================================================
# Sex-specific spline plots for MET_total
# =========================================================

sex_spline_plots <- make_sex_specific_spline_plots(
  data = IPAQ_dat2,
  outcome = "SRF",
  exposure_log = "MET_total_log",
  confounders = c("age_A0", "ethnicity", "weight", "height", "tdi", "education_level"),
  spline_knots = 4
)

p_met_male <- sex_spline_plots$male
p_met_female <- sex_spline_plots$female

p_met_male
p_met_female


# =========================================================
# Notes
# =========================================================

# Linearity of continuous physical activity measures was assessed using
# models adjusting for confounders but without interaction terms.
# This approach tests the overall functional form of the main exposure effect.
# Final models may include interaction terms where appropriate.
#
# Results demonstrating non-linearity can motivate modelling the exposure
# using splines or presenting it categorically (e.g. quintiles) for interpretability.

###################################################
#### Main analysis ########
###################################################

# =========================================================
# Quintile analysis
# =========================================================

# Primary approach:
# - Quintiles are based on the RAW exposure
# - Logistic regression uses quintile groups with Q1 as reference
# - Progressive adjustment models
# - Stratified analyses by sex, by age group, and by sex + age group



# =========================================================
# Analysis settings
# =========================================================

n_bins <- 5

# Primary quintile analysis
quintile_exposures_primary <- c("MET_total")
quintile_outcomes_primary <- c("SRF")

# Age-sex combined stratified analysis
quintile_exposures_age_sex <- c("MET_total")
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
# 1. Quintile analysis stratified by sex
# =========================================================

sex_strata <- data.frame(
  sex = c("Male", "Female"),
  stringsAsFactors = FALSE
)

results_quintiles_sex <- run_quintile_models_stratified(
  data = IPAQ_dat2,
  strata_df = sex_strata,
  outcomes = quintile_outcomes_primary,
  exposures = quintile_exposures_primary,
  model_specs = model_specs_sex,
  n_bins = n_bins
) %>%
  format_quintile_results(model_levels = model_levels, digits = 5)

saveRDS(
  results_quintiles_sex,
  file.path(DATA_DERIVED, "results_quintiles_sex.rds")
)


# =========================================================
# 2. Quintile analysis stratified by age group
# =========================================================

age_levels <- setdiff(unique(as.character(IPAQ_dat2$agegp_A0)), "70-79")

age_strata <- data.frame(
  agegp_A0 = age_levels,
  stringsAsFactors = FALSE
)

results_quintiles_age <- run_quintile_models_stratified(
  data = IPAQ_dat2,
  strata_df = age_strata,
  outcomes = quintile_outcomes_primary,
  exposures = quintile_exposures_primary,
  model_specs = model_specs_age,
  n_bins = n_bins
) %>%
  format_quintile_results(model_levels = model_levels, digits = 5)

saveRDS(
  results_quintiles_age,
  file.path(DATA_DERIVED, "results_quintiles_age.rds")
)



# =========================================================
# 3. Quintile analysis stratified by sex and age group
# =========================================================

sex_levels <- unique(as.character(IPAQ_dat2$sex))
age_levels <- setdiff(unique(as.character(IPAQ_dat2$agegp_A0)), "70-79")

strata_age_sex <- expand.grid(
  sex = sex_levels,
  agegp_A0 = age_levels,
  stringsAsFactors = FALSE
)

results_quintiles_age_sex <- run_quintile_models_stratified(
  data = IPAQ_dat2,
  strata_df = strata_age_sex,
  outcomes = quintile_outcomes_age_sex,
  exposures = quintile_exposures_age_sex,
  model_specs = model_specs_age_sex,
  n_bins = n_bins
) %>%
  format_quintile_results(model_levels = model_levels, digits = 5)

saveRDS(
  results_quintiles_age_sex,
  file.path(DATA_DERIVED, "results_quintiles_age_sex.rds")
)


# =========================================================
# Sensitivity analysis: quintiles for alternative PA exposures
# =========================================================

sens_exposures <- c("MET_total", "MET_mod", "MET_vig", "MET_walk", "MET_MVPA")
sens_outcomes <- c("SRF")

sex_strata <- data.frame(
  sex = c("Male", "Female"),
  stringsAsFactors = FALSE
)

sens_model_specs <- list(
  "Fully adjusted" = c("age_A0", "ethnicity", "weight", "height", "tdi", "education_level")
)

results_sex_only <- run_categorical_models_stratified(
  data = IPAQ_dat2,
  strata_df = sex_strata,
  outcomes = sens_outcomes,
  exposures = sens_exposures,
  model_specs = sens_model_specs,
  cat_type = "quintile",
  n_bins = 5
) %>%
  dplyr::mutate(
    OR = round(OR, 4),
    CI_lower = round(CI_lower, 4),
    CI_upper = round(CI_upper, 4),
    p_value = round(p_value, 4)
  )

saveRDS(
  results_sex_only,
  file.path(DATA_DERIVED, "results_sensitivity_quintiles_sex.rds")
)


# =========================================================
# Mutually adjusted MET intensity models, by sex
# =========================================================

mutual_exposures <- c("MET_vig", "MET_mod", "MET_walk")
sens_outcomes <- c("SRF")

mutual_model_specs <- list(
  "Mutually adjusted" = c(
    "age_A0", "ethnicity", "weight", "height",
    "tdi", "education_level"
  )
)

run_mutual_model <- function(exposure_var) {
  
  other_activity_vars <- setdiff(mutual_exposures, exposure_var)
  
  model_specs <- list(
    "Mutually adjusted" = c(
      other_activity_vars,
      "age_A0", "ethnicity", "weight", "height",
      "tdi", "education_level"
    )
  )
  
  run_categorical_models_stratified(
    data = IPAQ_dat2,
    strata_df = sex_strata,
    outcomes = sens_outcomes,
    exposures = exposure_var,
    model_specs = model_specs,
    cat_type = "quintile",
    n_bins = 5
  )
}

results_mutual <- purrr::map_dfr(
  mutual_exposures,
  run_mutual_model
)

saveRDS(
  results_mutual,
  file.path(DATA_DERIVED, "results_mutual_quintiles_sex.rds")
)

# =========================================================
# Mutually adjusted MET intensity models, by sex and age group
# =========================================================

mutual_exposures_age_sex <- c("MET_vig", "MET_mod", "MET_walk")
sens_outcomes_age_sex <- c("SRF")

sex_age_strata <- IPAQ_dat2 %>%
  dplyr::filter(
    !is.na(sex),
    !is.na(agegp_A0),
    agegp_A0 %in% c("40-49", "50-59", "60-65")
  ) %>%
  dplyr::distinct(sex, agegp_A0) %>%
  dplyr::arrange(sex, agegp_A0)

run_mutual_model_age_sex <- function(exposure_var) {
  
  other_activity_vars <- setdiff(mutual_exposures_age_sex, exposure_var)
  
  model_specs <- list(
    "Mutually adjusted" = c(
      other_activity_vars,
      "ethnicity", "weight", "height",
      "tdi", "education_level"
    )
  )
  
  run_categorical_models_stratified(
    data = IPAQ_dat2,
    strata_df = sex_age_strata,
    outcomes = sens_outcomes_age_sex,
    exposures = exposure_var,
    model_specs = model_specs,
    cat_type = "quintile",
    n_bins = 5
  )
}

results_mutual_age_sex <- purrr::map_dfr(
  mutual_exposures_age_sex,
  run_mutual_model_age_sex
)

saveRDS(
  results_mutual_age_sex,
  file.path(DATA_DERIVED, "results_mutual_quintiles_age_sex.rds")
)


# =========================================================
# Spline sensitivity analysis by sex and fracture type
# =========================================================


spline_confounders <- c("age_A0", "ethnicity", "weight", "height", "tdi", "education_level")

spline_plot_data <- dplyr::bind_rows(
  get_spline_prediction_data(IPAQ_dat2, "SRF", "MET_total_log", spline_confounders, strata = list(sex = "Male")),
  get_spline_prediction_data(IPAQ_dat2, "SRF", "MET_total_log", spline_confounders, strata = list(sex = "Female")),
  get_spline_prediction_data(IPAQ_dat2, "Hip", "MET_total_log", spline_confounders, strata = list(sex = "Male")),
  get_spline_prediction_data(IPAQ_dat2, "Hip", "MET_total_log", spline_confounders, strata = list(sex = "Female")),
  get_spline_prediction_data(IPAQ_dat2, "Wrist", "MET_total_log", spline_confounders, strata = list(sex = "Male")),
  get_spline_prediction_data(IPAQ_dat2, "Wrist", "MET_total_log", spline_confounders, strata = list(sex = "Female"))
) %>%
  dplyr::mutate(
    outcome = factor(outcome, levels = c("SRF", "Hip", "Wrist")),
    sex = factor(sex, levels = c("Male", "Female"))
  )

saveRDS(
  spline_plot_data,
  file.path(DATA_DERIVED, "spline_plot_data_by_sex_fracture_type.rds")
)


# =========================================================
# Sensitivity analysis: quintiles for Hip and Wrist
# Stratified by sex
# =========================================================

quintile_exposures_sens <- c("MET_total")
quintile_outcomes_sens <- c("Hip", "Wrist")

sex_strata <- data.frame(
  sex = c("Male", "Female"),
  stringsAsFactors = FALSE
)

results_quintiles_sex_sens <- run_quintile_models_stratified(
  data = IPAQ_dat2,
  strata_df = sex_strata,
  outcomes = quintile_outcomes_sens,
  exposures = quintile_exposures_sens,
  model_specs = model_specs_sex,
  n_bins = n_bins
) %>%
  format_quintile_results(model_levels = model_levels, digits = 5)

saveRDS(
  results_quintiles_sex_sens,
  file.path(DATA_DERIVED, "results_quintiles_sex_hip_wrist.rds")
)


# =========================================================
# Sensitivity analysis: alternative exposures
# MET_mod, MET_vig, MET_walk
# Outcome = SRF
# Stratified by sex
# =========================================================

quintile_exposures_sens_exp <- c("MET_mod", "MET_vig", "MET_walk")
quintile_outcomes_sens_exp <- c("SRF")

sex_strata <- data.frame(
  sex = c("Male", "Female"),
  stringsAsFactors = FALSE
)

results_quintiles_sex_exposure_sens <- run_quintile_models_stratified(
  data = IPAQ_dat2,
  strata_df = sex_strata,
  outcomes = quintile_outcomes_sens_exp,
  exposures = quintile_exposures_sens_exp,
  model_specs = model_specs_sex,
  n_bins = n_bins
) %>%
  format_quintile_results(model_levels = model_levels, digits = 5)

saveRDS(
  results_quintiles_sex_exposure_sens,
  file.path(DATA_DERIVED, "results_quintiles_sex_exposure_sens.rds")
)

# =========================================================
# Activity within quintiles variables
# =========================================================

sex_levels <- c("Female", "Male")
confounders_full <- c("age_A0", "ethnicity", "weight", "height", "tdi", "education_level")

IPAQ_dat2_prop <- IPAQ_dat2 %>%
  dplyr::mutate(
    MET_total_check = MET_walk + MET_mod + MET_vig,
    
    prop_walk = dplyr::if_else(MET_total_check > 0, MET_walk / MET_total_check, NA_real_),
    prop_mod  = dplyr::if_else(MET_total_check > 0, MET_mod  / MET_total_check, NA_real_),
    prop_vig  = dplyr::if_else(MET_total_check > 0, MET_vig  / MET_total_check, NA_real_),
    
    prop_walk_10 = prop_walk * 10,
    prop_mod_10  = prop_mod * 10,
    prop_vig_10  = prop_vig * 10
  )

IPAQ_dat2_prop %>%
  dplyr::summarise(
    min_prop_sum  = min(prop_walk + prop_mod + prop_vig, na.rm = TRUE),
    max_prop_sum  = max(prop_walk + prop_mod + prop_vig, na.rm = TRUE),
    mean_prop_sum = mean(prop_walk + prop_mod + prop_vig, na.rm = TRUE)
  )

# =========================================================
# Activity within MET_total quintiles
# =========================================================

results_prop_within_q_df <- run_prop_within_quintiles_by_sex(
  data = IPAQ_dat2_prop,
  sex_levels = sex_levels,
  confounders = confounders_full
) %>%
  dplyr::mutate(
    OR = round(OR, 4),
    CI_lower = round(CI_lower, 4),
    CI_upper = round(CI_upper, 4),
    p_value = round(p_value, 4)
  )

saveRDS(
  IPAQ_dat2_prop,
  file.path(DATA_DERIVED, "IPAQ_dat2_prop.rds")
)

saveRDS(
  results_prop_within_q_df,
  file.path(DATA_DERIVED, "results_prop_quintile_activity.rds")
)

# =========================================================
# Weight/adiposity sensitivity analysis
# Fully adjusted models only for presentation
# =========================================================
# This section explores whether different specifications of weight affect the outcome

n_bins <- 5

# Primary quintile analysis
quintile_exposures_primary <- c("MET_total")
quintile_outcomes_primary <- c("SRF")


# Check variable names before running:
# ethnicity_derived, WC, BMI, weight, height should all exist in IPAQ_dat2

model_specs_weight_height <- list(
  "Fully adjusted" = c(
    "age_A0", "sex", "ethnicity",
    "weight", "height",
    "tdi", "education_level"
  )
)

model_specs_BMI <- list(
  "Fully adjusted" = c(
    "age_A0", "sex", "ethnicity",
    "BMI",
    "tdi", "education_level"
  )
)

model_specs_WC_height <- list(
  "Fully adjusted" = c(
    "age_A0", "sex", "ethnicity",
    "WC", "height",
    "tdi", "education_level"
  )
)

adiposity_model_specs <- list(
  "Weight + height" = model_specs_weight_height,
  "BMI" = model_specs_BMI,
  "Waist circumference + height" = model_specs_WC_height
)

results_weight_sensitivity <- list()
idx <- 1

for (adiposity_model in names(adiposity_model_specs)) {
  
  model_specs_current <- adiposity_model_specs[[adiposity_model]]
  
  for (outcome in quintile_outcomes_primary) {
    for (exposure in quintile_exposures_primary) {
      for (model_name in names(model_specs_current)) {
        
        res <- fit_quintile_model(
          data = IPAQ_dat2,
          exposure = exposure,
          outcome = outcome,
          confounders = model_specs_current[[model_name]],
          n_bins = n_bins
        )
        
        if (is.null(res)) next
        
        res$model <- model_name
        res$adiposity_adjustment <- adiposity_model
        
        results_weight_sensitivity[[idx]] <- res
        idx <- idx + 1
      }
    }
  }
}

results_weight_sensitivity <- dplyr::bind_rows(results_weight_sensitivity) %>%
  dplyr::mutate(
    model = factor(model, levels = "Fully adjusted"),
    adiposity_adjustment = factor(
      adiposity_adjustment,
      levels = c(
        "Weight + height",
        "BMI",
        "Waist circumference + height"
      )
    ),
    Odds_Ratio = round(Odds_Ratio, 5),
    CI_lower = round(CI_lower, 5),
    CI_upper = round(CI_upper, 5),
    p_value = round(p_value, 5),
    Global_p_value = round(Global_p_value, 5)
  ) %>%
  dplyr::select(
    adiposity_adjustment,
    model,
    outcome,
    exposure,
    term,
    Odds_Ratio,
    CI_lower,
    CI_upper,
    p_value,
    Global_p_value
  )

saveRDS(
  results_weight_sensitivity,
  file.path(DATA_DERIVED, "results_weight_sensitivity.rds")
)


# =========================================================
# Mapping to Public Health - WHO-category 
# =========================================================

who_main <- fit_categorical_model(
  data = IPAQ_dat2,
  outcome = "SRF",
  exposure = "WHO_MVPA_equiv_cat",
  confounders = c("age_A0", "sex", "ethnicity", "weight", "height", "tdi", "education_level")
) %>%
  dplyr::mutate(
    model = "Fully adjusted"
  )

saveRDS(
  who_main,
  file.path(DATA_DERIVED, "who_main.rds")
)


# =========================================================
# WHO categories by sex
# =========================================================

who_model_specs_sex <- list(
  "Unadjusted" = character(0),
  "Minimal adjustment" = c("age_A0", "weight", "height"),
  "Fully adjusted" = c("age_A0", "ethnicity", "weight", "height", "tdi", "education_level")
)

results_WHO_sex_df <- run_categorical_models_stratified(
  data = IPAQ_dat2,
  strata_df = data.frame(sex = c("Male", "Female")),
  outcomes = c("SRF"),
  exposures = c("WHO_MVPA_equiv_cat"),
  model_specs = who_model_specs_sex,
  cat_type = "precat"
)

saveRDS(
  results_WHO_sex_df,
  file.path(DATA_DERIVED, "results_WHO_sex.rds")
)


# =========================================================
# WHO categories by sex and age group
# =========================================================

strata_age_sex <- expand.grid(
  sex = c("Male", "Female"),
  agegp_A0 = c("40-49", "50-59", "60-65"),
  stringsAsFactors = FALSE
)

who_model_specs_age_sex <- list(
  "Unadjusted" = character(0),
  "Minimal adjustment" = c("weight", "height"),
  "Fully adjusted" = c("ethnicity", "weight", "height", "tdi", "education_level")
)

results_WHO_age_sex_df <- run_categorical_models_stratified(
  data = IPAQ_dat2,
  strata_df = strata_age_sex,
  outcomes = c("SRF"),
  exposures = c("WHO_MVPA_equiv_cat"),
  model_specs = who_model_specs_age_sex,
  cat_type = "precat"
) %>%
  dplyr::mutate(
    model = factor(
      model,
      levels = c("Unadjusted", "Minimal adjustment", "Fully adjusted")
    )
  )

saveRDS(
  results_WHO_age_sex_df,
  file.path(DATA_DERIVED, "results_WHO_age_sex.rds")
)





