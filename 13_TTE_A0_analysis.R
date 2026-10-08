# =========================================================
# TIME-TO-EVENT ANALYSIS
# PHYSICAL ACTIVITY AND INCIDENT HIP FRACTURE
# =========================================================
#
# This script performs analysis only.
# It saves model/result objects as .rds for use in separate
# figures and tables scripts.
#
# Analysis sequence
# -----------------
# 1. Prepare exposure variables
# 2. Pooled total-LTPA continuous models
# 3. Proportional-hazards assessment
# 4. Sex-specific total-LTPA continuous models
# 5. Sex-specific total-LTPA quintile models
# 6. Sex-specific PA-component quintile models
# 7. Fine-Gray competing-risk sensitivity
# 8. Health-marker sensitivity:
#      total / walking / moderate / vigorous PA
#      × falls / grip / walking pace / heel BMD
#
# =========================================================


# =========================================================
# 0. SETUP
# =========================================================

library(here)

source(here::here("scripts/00_setup.R"))
source(here::here("scripts/03_helpers_general.R"))
source(here::here("scripts/04_helpers_CXA.R"))
source(here::here("scripts/05_helpers_survival.R"))


# =========================================================
# 1. LOAD ANALYSIS DATA
# =========================================================

TTE_A0_analysis <- readRDS(
  file.path(
    DATA_DERIVED,
    "TTE_A0_complete_case.Rds"
  )
)


# =========================================================
# 2. PREPARE ANALYSIS VARIABLES
# =========================================================
#
# Total LTPA:
#   cc_MET_total_trunc = complete-case MET-min/week,
#   truncated for implausibly high values.
#
# Total-PA quintiles are created in the pooled analysis cohort,
# matching the original primary analysis, and then retained
# when the data are split by sex.
# =========================================================

TTE_A0_analysis <- TTE_A0_analysis %>%
  dplyr::mutate(
    log_MET_total =
      log(cc_MET_total_trunc + 1),
    
    sex_raw = factor(sex_raw),
    ethnicity_derived =
      factor(ethnicity_derived),
    education_level =
      factor(education_level)
  )

TTE_A0_analysis <- add_pa_quintile(
  data = TTE_A0_analysis,
  exposure_var = "cc_MET_total_trunc",
  quintile_var = "MET_total_quintile_f"
)


# -----------------------------
# Basic cohort checks
# -----------------------------

hip_event_counts <- TTE_A0_analysis %>%
  dplyr::count(
    sex_raw,
    event_hip_PA,
    name = "n"
  )

print(hip_event_counts)

print(
  table(
    TTE_A0_analysis$MET_total_quintile_f,
    useNA = "ifany"
  )
)


# =========================================================
# 3. POOLED CONTINUOUS TOTAL-LTPA MODELS
# =========================================================

res_hip_cont <- run_pa_whole_followup(
  data = TTE_A0_analysis,
  outcome_name = "Hip fracture",
  age_exit_var = "age_exit_hip_PA",
  event_var = "event_hip_PA",
  make_ci_tables = FALSE
)

saveRDS(
  res_hip_cont,
  file.path(
    DATA_DERIVED,
    "res_hip_cont.rds"
  )
)


# -----------------------------
# Review model comparisons
# -----------------------------

cat(
  "\n--- Pooled linear vs spline ---\n"
)

print(
  res_hip_cont$
    comparisons$
    linear_vs_spline_lrt
)

print(
  res_hip_cont$
    comparisons$
    linear_vs_spline_aic
)

cat(
  "\n--- Pooled spline vs sex interaction ---\n"
)

print(
  res_hip_cont$
    comparisons$
    sex_interaction_lrt
)

print(
  res_hip_cont$
    comparisons$
    sex_interaction_aic
)


# =========================================================
# 4. PROPORTIONAL-HAZARDS CHECK:
#    POOLED FULLY ADJUSTED SPLINE
# =========================================================
#
# Sex violated the PH assumption in the original analysis,
# motivating presentation of sex-specific models.
# =========================================================

ph_hip_full_spline_model <-
  survival::coxph(
    formula =
      res_hip_cont$
      formulas$
      fully_adjusted_spline,
    data = TTE_A0_analysis,
    x = TRUE
  )

ph_test_hip_full_spline <-
  survival::cox.zph(
    ph_hip_full_spline_model
  )

ph_hip_pooled <- list(
  full_spline =
    ph_test_hip_full_spline
)

saveRDS(
  ph_hip_pooled,
  file.path(
    DATA_DERIVED,
    "ph_hip_pooled.rds"
  )
)

print(ph_test_hip_full_spline)


# =========================================================
# 5. CREATE SEX-SPECIFIC ANALYSIS DATASETS
# =========================================================

sex_levels <- levels(
  TTE_A0_analysis$sex_raw
)

if (
  !all(
    c("Female", "Male") %in%
    sex_levels
  )
) {
  stop(
    "Expected sex_raw to contain Female and Male."
  )
}

female_dat <- TTE_A0_analysis %>%
  dplyr::filter(
    sex_raw == "Female"
  ) %>%
  droplevels()

male_dat <- TTE_A0_analysis %>%
  dplyr::filter(
    sex_raw == "Male"
  ) %>%
  droplevels()


# =========================================================
# 6. CREATE SEX-SPECIFIC PA-COMPONENT QUINTILES
# =========================================================
#
# This prevents marker missingness from redefining Q1-Q5.
# =========================================================

add_component_quintiles <- function(dat) {
  
  dat <- add_pa_quintile(
    dat,
    "cc_MET_walk_trunc",
    "MET_walk_quintile_f"
  )
  
  dat <- add_pa_quintile(
    dat,
    "cc_MET_mod_trunc",
    "MET_mod_quintile_f"
  )
  
  dat <- add_pa_quintile(
    dat,
    "cc_MET_vig_trunc",
    "MET_vig_quintile_f"
  )
  
  dat
}

female_dat <-
  add_component_quintiles(
    female_dat
  )

male_dat <-
  add_component_quintiles(
    male_dat
  )


# Save sex-specific analysis datasets.
# These contain all retained PA quintile definitions and can
# also be loaded by downstream figures/tables scripts.

saveRDS(
  female_dat,
  file.path(
    DATA_DERIVED,
    "TTE_A0_analysis_female.rds"
  )
)

saveRDS(
  male_dat,
  file.path(
    DATA_DERIVED,
    "TTE_A0_analysis_male.rds"
  )
)


sex_specific_counts <-
  dplyr::bind_rows(
    
    female_dat %>%
      dplyr::summarise(
        sex = "Female",
        n = dplyr::n(),
        hip_events =
          sum(
            event_hip_PA,
            na.rm = TRUE
          )
      ),
    
    male_dat %>%
      dplyr::summarise(
        sex = "Male",
        n = dplyr::n(),
        hip_events =
          sum(
            event_hip_PA,
            na.rm = TRUE
          )
      )
  )

print(sex_specific_counts)


# =========================================================
# 7. SEX-SPECIFIC CONTINUOUS TOTAL-LTPA MODELS
# =========================================================

res_hip_female <-
  run_pa_whole_followup_sex_specific(
    data = female_dat,
    outcome_name =
      "Hip fracture - Female",
    age_exit_var =
      "age_exit_hip_PA",
    event_var =
      "event_hip_PA",
    make_ci_tables = FALSE
  )

res_hip_male <-
  run_pa_whole_followup_sex_specific(
    data = male_dat,
    outcome_name =
      "Hip fracture - Male",
    age_exit_var =
      "age_exit_hip_PA",
    event_var =
      "event_hip_PA",
    make_ci_tables = FALSE
  )

saveRDS(
  res_hip_female,
  file.path(
    DATA_DERIVED,
    "res_hip_female.rds"
  )
)

saveRDS(
  res_hip_male,
  file.path(
    DATA_DERIVED,
    "res_hip_male.rds"
  )
)

# =========================================================
# 7A. CREATE SPLINE PREDICTION DATA FOR FIGURE
# =========================================================
#
# Fully adjusted restricted cubic spline models.
# HRs are expressed relative to the median log_MET_total.
# Prediction range is restricted to the 1st–99th percentiles.
#
# Confidence intervals are calculated for the contrast:
#   LP(x) - LP(reference)
# using the model variance-covariance matrix.
# =========================================================


make_cox_spline_plot_data <- function(
    model,
    data,
    sex_label
) {
  
  # -------------------------------------------------------
  # Exposure range for plotting
  # -------------------------------------------------------
  
  x_seq <- seq(
    stats::quantile(
      data$log_MET_total,
      0.01,
      na.rm = TRUE
    ),
    stats::quantile(
      data$log_MET_total,
      0.99,
      na.rm = TRUE
    ),
    length.out = 100
  )
  
  
  # -------------------------------------------------------
  # Reference covariate profile
  # Exposure reference = median log_MET_total
  # -------------------------------------------------------
  
  ref_data <- data[1, , drop = FALSE]
  
  ref_data$log_MET_total <- stats::median(
    data$log_MET_total,
    na.rm = TRUE
  )
  
  ref_data$ethnicity_derived <- names(
    sort(
      table(data$ethnicity_derived),
      decreasing = TRUE
    )
  )[1]
  
  ref_data$height_clean <- stats::median(
    data$height_clean,
    na.rm = TRUE
  )
  
  ref_data$weight_clean <- stats::median(
    data$weight_clean,
    na.rm = TRUE
  )
  
  ref_data$tdi_raw <- stats::median(
    data$tdi_raw,
    na.rm = TRUE
  )
  
  ref_data$education_level <- names(
    sort(
      table(data$education_level),
      decreasing = TRUE
    )
  )[1]
  
  
  # Preserve factor levels
  
  ref_data$ethnicity_derived <- factor(
    ref_data$ethnicity_derived,
    levels = levels(data$ethnicity_derived)
  )
  
  ref_data$education_level <- factor(
    ref_data$education_level,
    levels = levels(data$education_level),
    ordered = is.ordered(data$education_level)
  )
  
  
  # -------------------------------------------------------
  # Prediction dataset
  # -------------------------------------------------------
  
  new_data <- ref_data[
    rep(1, length(x_seq)),
    ,
    drop = FALSE
  ]
  
  new_data$log_MET_total <- x_seq
  
  
  # -------------------------------------------------------
  # Obtain model design matrices
  # -------------------------------------------------------
  
  X_new <- stats::predict(
    model,
    newdata = new_data,
    type = "terms"
  )
  
  X_ref <- stats::predict(
    model,
    newdata = ref_data,
    type = "terms"
  )
  
  
  # -------------------------------------------------------
  # Use the Cox model X matrix directly to ensure that
  # spline basis terms exactly match the fitted model
  # -------------------------------------------------------
  
  beta <- stats::coef(model)
  V <- stats::vcov(model)
  
  X_new_full <- stats::model.matrix(
    stats::delete.response(stats::terms(model)),
    data = new_data
  )
  
  X_ref_full <- stats::model.matrix(
    stats::delete.response(stats::terms(model)),
    data = ref_data
  )
  
  # Keep only columns corresponding to fitted coefficients
  
  X_new_full <- X_new_full[
    ,
    names(beta),
    drop = FALSE
  ]
  
  X_ref_full <- X_ref_full[
    ,
    names(beta),
    drop = FALSE
  ]
  
  
  # -------------------------------------------------------
  # Contrast each exposure value with median reference
  # -------------------------------------------------------
  
  X_diff <- sweep(
    X_new_full,
    2,
    X_ref_full[1, ],
    FUN = "-"
  )
  
  log_HR <- as.vector(
    X_diff %*% beta
  )
  
  var_log_HR <- rowSums(
    (X_diff %*% V) * X_diff
  )
  
  se_log_HR <- sqrt(
    pmax(var_log_HR, 0)
  )
  
  
  # -------------------------------------------------------
  # Convert to HR and 95% CI
  # -------------------------------------------------------
  
  tibble::tibble(
    Sex = sex_label,
    log_MET_total = x_seq,
    HR = exp(log_HR),
    HR_low = exp(
      log_HR - 1.96 * se_log_HR
    ),
    HR_high = exp(
      log_HR + 1.96 * se_log_HR
    )
  )
}


# =========================================================
# Create female and male prediction datasets
# =========================================================

hip_spline_plot_data <- list(
  
  female = make_cox_spline_plot_data(
    model =
      res_hip_female$models$fully_adjusted_spline,
    data =
      female_dat,
    sex_label =
      "Female"
  ),
  
  male = make_cox_spline_plot_data(
    model =
      res_hip_male$models$fully_adjusted_spline,
    data =
      male_dat,
    sex_label =
      "Male"
  )
)


# =========================================================
# Save for figures script
# =========================================================

saveRDS(
  hip_spline_plot_data,
  file.path(
    DATA_DERIVED,
    "hip_spline_plot_data.rds"
  )
)


# -----------------------------
# Review non-linearity
# -----------------------------

cat(
  "\n--- Female linear vs spline ---\n"
)

print(
  res_hip_female$
    comparisons$
    linear_vs_spline_lrt
)

print(
  res_hip_female$
    comparisons$
    linear_vs_spline_aic
)

cat(
  "\n--- Male linear vs spline ---\n"
)

print(
  res_hip_male$
    comparisons$
    linear_vs_spline_lrt
)

print(
  res_hip_male$
    comparisons$
    linear_vs_spline_aic
)


# =========================================================
# 8. SEX-SPECIFIC PH CHECKS
# =========================================================

female_ph_spline_model <-
  survival::coxph(
    formula =
      res_hip_female$
      formulas$
      fully_adjusted_spline,
    data = female_dat,
    x = TRUE
  )

male_ph_spline_model <-
  survival::coxph(
    formula =
      res_hip_male$
      formulas$
      fully_adjusted_spline,
    data = male_dat,
    x = TRUE
  )

ph_test_female_spline <-
  survival::cox.zph(
    female_ph_spline_model
  )

ph_test_male_spline <-
  survival::cox.zph(
    male_ph_spline_model
  )

ph_hip_sex_specific <- list(
  female_spline =
    ph_test_female_spline,
  male_spline =
    ph_test_male_spline
)

saveRDS(
  ph_hip_sex_specific,
  file.path(
    DATA_DERIVED,
    "ph_hip_sex_specific.rds"
  )
)


# =========================================================
# 9. SEX-SPECIFIC TOTAL-LTPA QUINTILE MODELS
# =========================================================

res_hip_quint_female <-
  run_pa_quintile_analysis_sex_specific(
    data = female_dat,
    outcome_name =
      "Hip fracture - Female",
    age_exit_var =
      "age_exit_hip_PA",
    event_var =
      "event_hip_PA",
    quintile_var =
      "MET_total_quintile_f",
    make_ci_tables = TRUE
  )

res_hip_quint_male <-
  run_pa_quintile_analysis_sex_specific(
    data = male_dat,
    outcome_name =
      "Hip fracture - Male",
    age_exit_var =
      "age_exit_hip_PA",
    event_var =
      "event_hip_PA",
    quintile_var =
      "MET_total_quintile_f",
    make_ci_tables = TRUE
  )

saveRDS(
  res_hip_quint_female,
  file.path(
    DATA_DERIVED,
    "res_hip_quint_female.rds"
  )
)

saveRDS(
  res_hip_quint_male,
  file.path(
    DATA_DERIVED,
    "res_hip_quint_male.rds"
  )
)


# =========================================================
# 10. SEX-SPECIFIC PA-COMPONENT QUINTILE ANALYSES
# =========================================================
#
# Primary component models:
#   walking
#   moderate
#   vigorous
#
# Mutually adjusted versions additionally include the other
# two PA components as continuous MET-min/week covariates.
# =========================================================

component_specs <- list(
  
  "Total PA" = list(
    quintile_var =
      "MET_total_quintile_f",
    adjust_for =
      character(0)
  ),
  
  "Walking PA" = list(
    quintile_var =
      "MET_walk_quintile_f",
    adjust_for =
      character(0)
  ),
  
  "Moderate PA" = list(
    quintile_var =
      "MET_mod_quintile_f",
    adjust_for =
      character(0)
  ),
  
  "Vigorous PA" = list(
    quintile_var =
      "MET_vig_quintile_f",
    adjust_for =
      character(0)
  ),
  
  "Walking PA mutually adjusted" =
    list(
      quintile_var =
        "MET_walk_quintile_f",
      adjust_for = c(
        "cc_MET_mod_trunc",
        "cc_MET_vig_trunc"
      )
    ),
  
  "Moderate PA mutually adjusted" =
    list(
      quintile_var =
        "MET_mod_quintile_f",
      adjust_for = c(
        "cc_MET_walk_trunc",
        "cc_MET_vig_trunc"
      )
    ),
  
  "Vigorous PA mutually adjusted" =
    list(
      quintile_var =
        "MET_vig_quintile_f",
      adjust_for = c(
        "cc_MET_walk_trunc",
        "cc_MET_mod_trunc"
      )
    )
)


res_hip_intensity_female <-
  run_pa_component_quintiles(
    data = female_dat,
    sex_label = "Female",
    outcome_name =
      "Hip fracture",
    age_exit_var =
      "age_exit_hip_PA",
    event_var =
      "event_hip_PA",
    exposure_specs =
      component_specs
  )

res_hip_intensity_male <-
  run_pa_component_quintiles(
    data = male_dat,
    sex_label = "Male",
    outcome_name =
      "Hip fracture",
    age_exit_var =
      "age_exit_hip_PA",
    event_var =
      "event_hip_PA",
    exposure_specs =
      component_specs
  )

res_hip_intensity <- list(
  female =
    res_hip_intensity_female,
  male =
    res_hip_intensity_male,
  combined =
    dplyr::bind_rows(
      res_hip_intensity_female,
      res_hip_intensity_male
    )
)

saveRDS(
  res_hip_intensity,
  file.path(
    DATA_DERIVED,
    "res_hip_intensity_quintile.rds"
  )
)


# =========================================================
# 11. COMPETING-RISK SENSITIVITY:
#     FINE-GRAY WITH DEATH AS COMPETING EVENT
# =========================================================
#
# Event coding:
#   0 = censored
#   1 = hip fracture
#   2 = death before hip fracture
#
# Fine-Gray uses follow-up time from baseline, as required by
# cmprsk::crr(), rather than age as the time scale.
# =========================================================

TTE_A0_analysis <- TTE_A0_analysis %>%
  dplyr::mutate(
    
    date_of_death =
      as.Date(date_of_death),
    
    date_admin_cens =
      as.Date(date_admin_cens),
    
    lost_to_fu_raw =
      as.Date(lost_to_fu_raw),
    
    date_first_hip_A0 =
      as.Date(date_first_hip_A0),
    
    date_assess_A0_raw =
      as.Date(date_assess_A0_raw),
    
    date_cens_no_death =
      pmin(
        lost_to_fu_raw,
        date_admin_cens,
        na.rm = TRUE
      ),
    
    event_hip_fg =
      dplyr::case_when(
        
        !is.na(date_first_hip_A0) &
          date_first_hip_A0 <=
          date_cens_no_death &
          (
            is.na(date_of_death) |
              date_first_hip_A0 <=
              date_of_death
          ) ~ 1L,
        
        !is.na(date_of_death) &
          date_of_death <=
          date_cens_no_death &
          (
            is.na(date_first_hip_A0) |
              date_of_death <
              date_first_hip_A0
          ) ~ 2L,
        
        TRUE ~ 0L
      ),
    
    date_exit_hip_fg =
      dplyr::case_when(
        event_hip_fg == 1L ~
          date_first_hip_A0,
        event_hip_fg == 2L ~
          date_of_death,
        TRUE ~
          date_cens_no_death
      ),
    
    fg_time_hip =
      as.numeric(
        date_exit_hip_fg -
          date_assess_A0_raw
      ) / 365.25
  )


if (
  any(
    TTE_A0_analysis$fg_time_hip < 0,
    na.rm = TRUE
  )
) {
  warning(
    "Negative Fine-Gray follow-up times identified."
  )
}


female_fg_dat <-
  TTE_A0_analysis %>%
  dplyr::filter(
    sex_raw == "Female"
  ) %>%
  droplevels()

male_fg_dat <-
  TTE_A0_analysis %>%
  dplyr::filter(
    sex_raw == "Male"
  ) %>%
  droplevels()


fg_hip_female <-
  run_crr_sex_specific(
    data = female_fg_dat,
    sex_label = "Female",
    time_var = "fg_time_hip",
    event_var = "event_hip_fg"
  )

fg_hip_male <-
  run_crr_sex_specific(
    data = male_fg_dat,
    sex_label = "Male",
    time_var = "fg_time_hip",
    event_var = "event_hip_fg"
  )


finegray_counts <-
  dplyr::bind_rows(
    
    female_fg_dat %>%
      dplyr::count(
        event_hip_fg,
        name = "n"
      ) %>%
      dplyr::mutate(
        sex = "Female"
      ),
    
    male_fg_dat %>%
      dplyr::count(
        event_hip_fg,
        name = "n"
      ) %>%
      dplyr::mutate(
        sex = "Male"
      )
  ) %>%
  dplyr::select(
    sex,
    event_hip_fg,
    n
  )


fg_hip_results <- list(
  
  female =
    fg_hip_female,
  
  male =
    fg_hip_male,
  
  female_data =
    female_fg_dat,
  
  male_data =
    male_fg_dat,
  
  event_counts =
    finegray_counts
)

saveRDS(
  fg_hip_results,
  file.path(
    DATA_DERIVED,
    "fg_hip_results.rds"
  )
)


# =========================================================
# 12. HEALTH-PROFILE SENSITIVITY ANALYSES
# =========================================================
#
# Question:
# Does additional adjustment for markers of falls /
# musculoskeletal function attenuate the association between
# PA and incident hip fracture?
#
# Exposures:
#   * Total PA
#   * Walking PA
#   * Moderate PA
#   * Vigorous PA
#
# Additional markers:
#   * Self-reported falls
#   * Grip strength
#   * Self-reported walking pace
#   * Heel BMD
#
# Each base vs marker-adjusted comparison:
#   * uses the same marker-complete sample
#   * uses PA quintiles defined before marker restriction
#   * is run separately in women and men
# =========================================================


# -----------------------------
# 12A. Derive health markers
# -----------------------------

TTE_A0_health_profile <-
  dplyr::bind_rows(
    female_dat,
    male_dat
  ) %>%
  dplyr::mutate(
    
    grip_strength_mean =
      rowMeans(
        dplyr::across(
          c(
            grip_strength_L_raw,
            grip_strength_R_raw
          )
        ),
        na.rm = TRUE
      ),
    
    grip_strength_mean =
      dplyr::if_else(
        is.nan(
          grip_strength_mean
        ),
        NA_real_,
        grip_strength_mean
      ),
    
    heel_BMD =
      heel_BMD_clean,
    
    falls_binary =
      dplyr::case_when(
        
        self_reported_falls_clean ==
          "No falls" ~
          "No falls",
        
        self_reported_falls_clean %in%
          c(
            "Only one fall",
            "More than one fall"
          ) ~
          "One or more falls",
        
        TRUE ~
          NA_character_
      ),
    
    falls_binary =
      factor(
        falls_binary,
        levels = c(
          "No falls",
          "One or more falls"
        )
      )
  )


saveRDS(
  TTE_A0_health_profile,
  file.path(
    DATA_DERIVED,
    "TTE_A0_health_profile.rds"
  )
)


# -----------------------------
# 12B. Exposure and marker specifications
# -----------------------------

health_pa_exposures <-
  tibble::tribble(
    ~exposure,      ~quintile_var,
    "Total PA",     "MET_total_quintile_f",
    "Walking PA",   "MET_walk_quintile_f",
    "Moderate PA",  "MET_mod_quintile_f",
    "Vigorous PA",  "MET_vig_quintile_f"
  )


health_markers <-
  tibble::tribble(
    ~marker,                 ~marker_var,
    "Falls",                 "falls_binary",
    "Grip strength",         "grip_strength_mean",
    "Walking pace",          "usual_walking_pace_raw",
    "Heel BMD",              "heel_BMD"
  )


# -----------------------------
# 12C. Run full exposure × marker matrix
# -----------------------------

res_hip_health_markers <-
  run_all_pa_marker_sensitivities(
    data =
      TTE_A0_health_profile,
    exposures =
      health_pa_exposures,
    markers =
      health_markers,
    sexes =
      c("Female", "Male"),
    base_confounders =
      default_tte_confounders(FALSE),
    age_entry_var =
      "age_entry_PA",
    age_exit_var =
      "age_exit_hip_PA",
    event_var =
      "event_hip_PA"
  )


# -----------------------------
# 12D. Review sample/event counts
# -----------------------------

print(
  res_hip_health_markers$
    sample_summary,
  n = Inf
)


# -----------------------------
# 12E. Save one structured result object
# -----------------------------

saveRDS(
  res_hip_health_markers,
  file.path(
    DATA_DERIVED,
    "res_hip_health_marker_sensitivity.rds"
  )
)


# =========================================================
# ANALYSIS COMPLETE
# =========================================================
#
# Saved objects required by downstream figures/tables:
#
# Saved objects required by downstream figures/tables:
#
# res_hip_cont.rds
# ph_hip_pooled.rds
# TTE_A0_analysis_female.rds
# TTE_A0_analysis_male.rds
# res_hip_female.rds
# res_hip_male.rds
# hip_spline_plot_data.rds
# ph_hip_sex_specific.rds
# res_hip_quint_female.rds
# res_hip_quint_male.rds
# res_hip_intensity_quintile.rds
# fg_hip_results.rds
# TTE_A0_health_profile.rds
# res_hip_health_marker_sensitivity.rds
#
# No tables or figures are generated in this script.
# =========================================================