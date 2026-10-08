# =========================================================
# SURVIVAL ANALYSIS HELPERS
# Physical activity and fracture
# =========================================================
#
# Design principles
# -----------------
# * This file contains FUNCTIONS ONLY.
# * No datasets are created or saved here.
# * No figures/tables are written here.
# * Analysis scripts define exposures/outcomes and save .rds outputs.
# * Figures/tables scripts format and plot saved analysis outputs.
#
# Main analysis structure
# -----------------------
# 1. Continuous total LTPA: linear and restricted cubic spline Cox models
# 2. Total LTPA quintiles
# 3. PA component quintiles (walking, moderate, vigorous), including mutual adjustment
# 4. Fine-Gray competing-risk sensitivity
# 5. Health-marker sensitivity:
#       total / walking / moderate / vigorous PA
#       × falls / grip / walking pace / heel BMD
#
# =========================================================


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


# =========================================================
# 0. COMMON UTILITIES
# =========================================================

default_tte_confounders <- function(include_sex = FALSE) {
  
  covars <- c(
    "ethnicity_derived",
    "height_clean",
    "weight_clean",
    "tdi_raw",
    "education_level"
  )
  
  if (include_sex) {
    covars <- c("sex_raw", covars)
  }
  
  covars
}


check_required_vars <- function(data, vars, context = "analysis") {
  
  missing_vars <- setdiff(unique(vars), names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      "Missing variables for ", context, ": ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  invisible(TRUE)
}


# Create quintiles once and retain them for subsequent restricted samples.
#
# This is particularly useful for sensitivity analyses: membership of Q1-Q5
# should not change simply because a marker is missing.

add_pa_quintile <- function(data, exposure_var, quintile_var) {
  
  check_required_vars(
    data,
    exposure_var,
    context = paste0("creation of ", quintile_var)
  )
  
  q_num <- dplyr::ntile(data[[exposure_var]], 5)
  
  data[[quintile_var]] <- factor(
    q_num,
    levels = 1:5,
    labels = paste0("Q", 1:5)
  )
  
  data
}


# =========================================================
# 1. CONTINUOUS TOTAL-LTPA COX MODELS
# =========================================================

run_pa_whole_followup <- function(
    data,
    outcome_name,
    age_exit_var,
    event_var,
    make_ci_tables = FALSE
) {
  
  required <- c(
    "age_entry_PA",
    age_exit_var,
    event_var,
    "log_MET_total",
    default_tte_confounders(include_sex = TRUE)
  )
  
  check_required_vars(data, required, "pooled continuous Cox models")
  
  surv_txt <- paste0(
    "survival::Surv(age_entry_PA, ",
    age_exit_var, ", ",
    event_var, ")"
  )
  
  f_unadj <- stats::as.formula(
    paste0(surv_txt, " ~ log_MET_total")
  )
  
  f_min <- stats::as.formula(
    paste0(
      surv_txt,
      " ~ log_MET_total + sex_raw + ethnicity_derived + height_clean + weight_clean"
    )
  )
  
  f_full_linear <- stats::as.formula(
    paste0(
      surv_txt,
      " ~ log_MET_total + ",
      paste(default_tte_confounders(include_sex = TRUE), collapse = " + ")
    )
  )
  
  f_full_spline <- stats::as.formula(
    paste0(
      surv_txt,
      " ~ rms::rcs(log_MET_total, 4) + ",
      paste(default_tte_confounders(include_sex = TRUE), collapse = " + ")
    )
  )
  
  f_spline_sexint <- stats::as.formula(
    paste0(
      surv_txt,
      " ~ rms::rcs(log_MET_total, 4) * sex_raw + ",
      paste(default_tte_confounders(include_sex = FALSE), collapse = " + ")
    )
  )
  
  m_unadj <- survival::coxph(
    f_unadj, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  m_min <- survival::coxph(
    f_min, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  m_full_linear <- survival::coxph(
    f_full_linear, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  m_full_spline <- survival::coxph(
    f_full_spline, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  m_spline_sexint <- survival::coxph(
    f_spline_sexint, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  tidied_models <- NULL
  
  if (make_ci_tables) {
    tidied_models <- list(
      unadjusted = broom::tidy(
        m_unadj, exponentiate = TRUE, conf.int = TRUE
      ),
      minimally_adjusted = broom::tidy(
        m_min, exponentiate = TRUE, conf.int = TRUE
      ),
      fully_adjusted_linear = broom::tidy(
        m_full_linear, exponentiate = TRUE, conf.int = TRUE
      ),
      fully_adjusted_spline = broom::tidy(
        m_full_spline, exponentiate = TRUE, conf.int = TRUE
      ),
      spline_sex_interaction = broom::tidy(
        m_spline_sexint, exponentiate = TRUE, conf.int = TRUE
      )
    )
  }
  
  list(
    outcome_name = outcome_name,
    exposure_type = "continuous",
    exposure_variable = "log_MET_total",
    
    formulas = list(
      unadjusted = f_unadj,
      minimally_adjusted = f_min,
      fully_adjusted_linear = f_full_linear,
      fully_adjusted_spline = f_full_spline,
      sex_interaction = f_spline_sexint
    ),
    
    models = list(
      unadjusted = m_unadj,
      minimally_adjusted = m_min,
      fully_adjusted_linear = m_full_linear,
      fully_adjusted_spline = m_full_spline,
      sex_interaction = m_spline_sexint
    ),
    
    tidied_models = tidied_models,
    
    comparisons = list(
      linear_vs_spline_lrt = stats::anova(
        m_full_linear,
        m_full_spline,
        test = "LRT"
      ),
      linear_vs_spline_aic = stats::AIC(
        m_full_linear,
        m_full_spline
      ),
      sex_interaction_lrt = stats::anova(
        m_full_spline,
        m_spline_sexint,
        test = "LRT"
      ),
      sex_interaction_aic = stats::AIC(
        m_full_spline,
        m_spline_sexint
      )
    )
  )
}


run_pa_whole_followup_sex_specific <- function(
    data,
    outcome_name,
    age_exit_var,
    event_var,
    make_ci_tables = FALSE
) {
  
  required <- c(
    "age_entry_PA",
    age_exit_var,
    event_var,
    "log_MET_total",
    default_tte_confounders(include_sex = FALSE)
  )
  
  check_required_vars(data, required, "sex-specific continuous Cox models")
  
  surv_txt <- paste0(
    "survival::Surv(age_entry_PA, ",
    age_exit_var, ", ",
    event_var, ")"
  )
  
  f_unadj <- stats::as.formula(
    paste0(surv_txt, " ~ log_MET_total")
  )
  
  f_min <- stats::as.formula(
    paste0(
      surv_txt,
      " ~ log_MET_total + ethnicity_derived + height_clean + weight_clean"
    )
  )
  
  f_full_linear <- stats::as.formula(
    paste0(
      surv_txt,
      " ~ log_MET_total + ",
      paste(default_tte_confounders(FALSE), collapse = " + ")
    )
  )
  
  f_full_spline <- stats::as.formula(
    paste0(
      surv_txt,
      " ~ rms::rcs(log_MET_total, 4) + ",
      paste(default_tte_confounders(FALSE), collapse = " + ")
    )
  )
  
  m_unadj <- survival::coxph(
    f_unadj, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  m_min <- survival::coxph(
    f_min, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  m_full_linear <- survival::coxph(
    f_full_linear, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  m_full_spline <- survival::coxph(
    f_full_spline, data = data, model = TRUE, x = TRUE, y = TRUE
  )
  
  tidied_models <- NULL
  
  if (make_ci_tables) {
    tidied_models <- list(
      unadjusted = broom::tidy(
        m_unadj, exponentiate = TRUE, conf.int = TRUE
      ),
      minimally_adjusted = broom::tidy(
        m_min, exponentiate = TRUE, conf.int = TRUE
      ),
      fully_adjusted_linear = broom::tidy(
        m_full_linear, exponentiate = TRUE, conf.int = TRUE
      ),
      fully_adjusted_spline = broom::tidy(
        m_full_spline, exponentiate = TRUE, conf.int = TRUE
      )
    )
  }
  
  list(
    outcome_name = outcome_name,
    exposure_type = "continuous",
    exposure_variable = "log_MET_total",
    
    formulas = list(
      unadjusted = f_unadj,
      minimally_adjusted = f_min,
      fully_adjusted_linear = f_full_linear,
      fully_adjusted_spline = f_full_spline
    ),
    
    models = list(
      unadjusted = m_unadj,
      minimally_adjusted = m_min,
      fully_adjusted_linear = m_full_linear,
      fully_adjusted_spline = m_full_spline
    ),
    
    tidied_models = tidied_models,
    
    comparisons = list(
      linear_vs_spline_lrt = stats::anova(
        m_full_linear,
        m_full_spline,
        test = "LRT"
      ),
      linear_vs_spline_aic = stats::AIC(
        m_full_linear,
        m_full_spline
      )
    )
  )
}


# =========================================================
# 2. GENERIC QUINTILE COX MODEL
# =========================================================
#
# `quintile_var` allows the same Cox machinery to be used for total,
# walking, moderate and vigorous PA without duplicate helpers.
#
# Default remains MET_total_quintile_f 
# =========================================================

run_pa_quintile_analysis_sex_specific <- function(
    data,
    outcome_name,
    age_exit_var,
    event_var,
    quintile_var = "MET_total_quintile_f",
    model_specs = NULL,
    make_ci_tables = FALSE,
    age_entry_var = "age_entry_PA"
) {
  
  check_required_vars(
    data,
    c(age_entry_var, age_exit_var, event_var, quintile_var),
    "quintile Cox model"
  )
  
  if (is.null(model_specs)) {
    model_specs <- list(
      "Unadjusted" = character(0),
      
      "Minimal adjustment" = c(
        "height_clean",
        "weight_clean"
      ),
      
      "Fully adjusted" = default_tte_confounders(FALSE)
    )
  }
  
  all_model_vars <- unique(unlist(model_specs))
  
  check_required_vars(
    data,
    all_model_vars,
    "quintile Cox model covariates"
  )
  
  surv_txt <- paste0(
    "survival::Surv(",
    age_entry_var, ", ",
    age_exit_var, ", ",
    event_var, ")"
  )
  
  fit_model <- function(covars) {
    
    rhs <- c(quintile_var, covars)
    
    form <- stats::as.formula(
      paste0(
        surv_txt,
        " ~ ",
        paste(rhs, collapse = " + ")
      )
    )
    
    survival::coxph(
      form,
      data = data,
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
  }
  
  quintile_models <- lapply(
    model_specs,
    fit_model
  )
  
  quintile_event_counts <- table(
    data[[quintile_var]],
    data[[event_var]],
    useNA = "ifany"
  )
  
  tidy_one_model <- function(model, model_name) {
    
    broom::tidy(
      model,
      exponentiate = TRUE,
      conf.int = TRUE
    ) %>%
      dplyr::filter(
        startsWith(term, quintile_var)
      ) %>%
      dplyr::mutate(
        quintile = sub(
          paste0("^", quintile_var),
          "",
          term
        ),
        model = model_name,
        outcome = outcome_name
      ) %>%
      dplyr::select(
        outcome,
        model,
        quintile,
        term,
        estimate,
        conf.low,
        conf.high,
        p.value
      )
  }
  
  tidy_quintile <- NULL
  
  if (make_ci_tables) {
    tidy_quintile <- dplyr::bind_rows(
      Map(
        tidy_one_model,
        quintile_models,
        names(quintile_models)
      )
    )
  }
  
  list(
    outcome_name = outcome_name,
    quintile_var = quintile_var,
    model_specs = model_specs,
    
    models = list(
      quintile = quintile_models
    ),
    
    tidied_models = list(
      quintile = tidy_quintile
    ),
    
    quintile_event_counts = quintile_event_counts
  )
}


# =========================================================
# 3. PA COMPONENT QUINTILE ANALYSES
# =========================================================
#
# `exposure_specs` contains:
#   label
#   quintile_var
#   optional continuous variables used for mutual adjustment
#
# Quintiles must already exist in `data`. They are not regenerated here.
# =========================================================

run_pa_component_quintiles <- function(
    data,
    sex_label,
    outcome_name,
    age_exit_var,
    event_var,
    exposure_specs,
    confounders = default_tte_confounders(FALSE),
    age_entry_var = "age_entry_PA"
) {
  
  results <- vector("list", length(exposure_specs))
  
  for (i in seq_along(exposure_specs)) {
    
    exposure_name <- names(exposure_specs)[i]
    spec <- exposure_specs[[i]]
    
    quintile_var <- spec$quintile_var
    adjust_for <- spec$adjust_for %||% character(0)
    
    required <- unique(c(
      age_entry_var,
      age_exit_var,
      event_var,
      quintile_var,
      confounders,
      adjust_for
    ))
    
    check_required_vars(
      data,
      required,
      paste0("component analysis: ", exposure_name)
    )
    
    dat <- data %>%
      dplyr::select(dplyr::all_of(required)) %>%
      tidyr::drop_na() %>%
      droplevels()
    
    if (nrow(dat) == 0 || dplyr::n_distinct(dat[[event_var]]) < 2) {
      warning("Skipping ", exposure_name, ": insufficient data/events.")
      next
    }
    
    model_specs <- list(
      "Fully adjusted" = c(confounders, adjust_for)
    )
    
    fit <- run_pa_quintile_analysis_sex_specific(
      data = dat,
      outcome_name = outcome_name,
      age_exit_var = age_exit_var,
      event_var = event_var,
      quintile_var = quintile_var,
      model_specs = model_specs,
      make_ci_tables = TRUE,
      age_entry_var = age_entry_var
    )
    
    tidy <- fit$tidied_models$quintile %>%
      dplyr::transmute(
        outcome,
        sex = sex_label,
        exposure = exposure_name,
        quintile,
        HR = estimate,
        CI_lower = conf.low,
        CI_upper = conf.high,
        p.value,
        n = stats::nobs(
          fit$models$quintile[["Fully adjusted"]]
        ),
        events =
          fit$models$quintile[["Fully adjusted"]]$nevent
      )
    
    results[[i]] <- tidy
  }
  
  dplyr::bind_rows(results)
}


# =========================================================
# 4. FINE-GRAY COMPETING-RISK MODEL
# =========================================================
#
# One generic helper replaces the duplicated run_crr_sex_specific().
# =========================================================

run_crr_sex_specific <- function(
    data,
    sex_label,
    time_var,
    event_var,
    nk = 4
) {
  
  required <- c(
    time_var,
    event_var,
    "log_MET_total",
    default_tte_confounders(FALSE)
  )
  
  check_required_vars(data, required, "Fine-Gray model")
  
  dat <- data %>%
    dplyr::filter(
      !is.na(.data[[time_var]]),
      .data[[time_var]] >= 0,
      !is.na(.data[[event_var]]),
      !is.na(log_MET_total),
      !is.na(ethnicity_derived),
      !is.na(height_clean),
      !is.na(weight_clean),
      !is.na(tdi_raw),
      !is.na(education_level)
    ) %>%
    droplevels()
  
  rcs_mat <- Hmisc::rcspline.eval(
    dat$log_MET_total,
    nk = nk,
    inclx = TRUE
  )
  
  colnames(rcs_mat) <- paste0(
    "rcs_log_MET_total_",
    seq_len(ncol(rcs_mat))
  )
  
  cov_data <- dplyr::bind_cols(
    as.data.frame(rcs_mat),
    dat %>%
      dplyr::select(
        ethnicity_derived,
        height_clean,
        weight_clean,
        tdi_raw,
        education_level
      )
  )
  
  cov_mat <- stats::model.matrix(
    ~ .,
    data = cov_data
  )[, -1, drop = FALSE]
  
  fg_model <- cmprsk::crr(
    ftime = dat[[time_var]],
    fstatus = dat[[event_var]],
    cov1 = cov_mat,
    failcode = 1,
    cencode = 0
  )
  
  tidy <- data.frame(
    sex = sex_label,
    term = names(fg_model$coef),
    estimate = exp(fg_model$coef),
    conf.low = exp(
      fg_model$coef -
        1.96 * sqrt(diag(fg_model$var))
    ),
    conf.high = exp(
      fg_model$coef +
        1.96 * sqrt(diag(fg_model$var))
    ),
    p.value = 2 * stats::pnorm(
      abs(
        fg_model$coef /
          sqrt(diag(fg_model$var))
      ),
      lower.tail = FALSE
    ),
    row.names = NULL
  )
  
  list(
    sex = sex_label,
    model = fg_model,
    tidy = tidy,
    knots = attr(rcs_mat, "knots"),
    nk = nk,
    cov_names = colnames(cov_mat),
    n = nrow(dat),
    event_counts = table(
      dat[[event_var]],
      useNA = "ifany"
    )
  )
}


# =========================================================
# 5. HEALTH-MARKER SENSITIVITY
# =========================================================
#
# Each comparison:
#   * uses a pre-existing PA quintile variable
#   * restricts to one common marker-complete sample
#   * fits the fully adjusted model
#   * fits the same model + marker
#
# Thus attenuation cannot be caused by different N between the two models,
# and the marker cannot alter the Q1-Q5 cut-points.
# =========================================================

run_pa_marker_sensitivity <- function(
    data,
    sex_label,
    exposure_label,
    quintile_var,
    marker_label,
    marker_var,
    base_confounders = default_tte_confounders(FALSE),
    sex_var = "sex_raw",
    age_entry_var = "age_entry_PA",
    age_exit_var = "age_exit_hip_PA",
    event_var = "event_hip_PA"
) {
  
  required <- unique(c(
    sex_var,
    age_entry_var,
    age_exit_var,
    event_var,
    quintile_var,
    marker_var,
    base_confounders
  ))
  
  check_required_vars(
    data,
    required,
    paste0(
      "marker sensitivity: ",
      sex_label, " / ",
      exposure_label, " / ",
      marker_label
    )
  )
  
  dat <- data %>%
    dplyr::filter(
      .data[[sex_var]] == sex_label
    ) %>%
    dplyr::select(
      dplyr::all_of(required)
    ) %>%
    tidyr::drop_na() %>%
    droplevels()
  
  if (
    nrow(dat) == 0 ||
    dplyr::n_distinct(dat[[event_var]]) < 2
  ) {
    stop(
      "Insufficient data/events for ",
      sex_label, " / ",
      exposure_label, " / ",
      marker_label
    )
  }
  
  model_specs <- list(
    "Fully adjusted" = base_confounders,
    "Marker adjusted" = c(
      base_confounders,
      marker_var
    )
  )
  
  fit <- run_pa_quintile_analysis_sex_specific(
    data = dat,
    outcome_name = "Hip fracture",
    age_exit_var = age_exit_var,
    event_var = event_var,
    quintile_var = quintile_var,
    model_specs = model_specs,
    make_ci_tables = TRUE,
    age_entry_var = age_entry_var
  )
  
  m_base <- fit$models$quintile[["Fully adjusted"]]
  m_marker <- fit$models$quintile[["Marker adjusted"]]
  
  if (stats::nobs(m_base) != stats::nobs(m_marker)) {
    stop("Nested health-marker models used different sample sizes.")
  }
  
  tidy <- fit$tidied_models$quintile %>%
    dplyr::mutate(
      sex = sex_label,
      exposure = exposure_label,
      marker = marker_label,
      quintile_var = quintile_var,
      marker_var = marker_var,
      model = dplyr::recode(
        model,
        "Marker adjusted" =
          paste0("Fully adjusted + ", marker_label)
      ),
      n = stats::nobs(m_base),
      events = m_base$nevent
    ) %>%
    dplyr::select(
      sex,
      exposure,
      marker,
      quintile_var,
      marker_var,
      model,
      quintile,
      estimate,
      conf.low,
      conf.high,
      p.value,
      n,
      events
    )
  
  # Nested-model comparison is retained as secondary descriptive output.
  # HR attenuation remains the primary sensitivity comparison.
  lrt <- stats::anova(
    m_base,
    m_marker,
    test = "LRT"
  )
  
  aic <- stats::AIC(
    m_base,
    m_marker
  )
  
  model_comparison <- tibble::tibble(
    sex = sex_label,
    exposure = exposure_label,
    marker = marker_label,
    n = stats::nobs(m_base),
    events = m_base$nevent,
    AIC_fully_adjusted =
      unname(aic["m_base", "AIC"]),
    AIC_marker_adjusted =
      unname(aic["m_marker", "AIC"]),
    delta_AIC =
      AIC_marker_adjusted -
      AIC_fully_adjusted,
    LRT_chi_square =
      unname(lrt$Chisq[2]),
    LRT_df =
      unname(lrt$Df[2]),
    LRT_p_value =
      unname(lrt$`Pr(>|Chi|)`[2])
  )
  
  hr_comparison <- tidy %>%
    dplyr::select(
      sex,
      exposure,
      marker,
      quintile,
      model,
      estimate,
      conf.low,
      conf.high
    ) %>%
    tidyr::pivot_wider(
      names_from = model,
      values_from = c(
        estimate,
        conf.low,
        conf.high
      )
    )
  
  # Create stable machine-friendly names without depending on marker text.
  base_name <- "Fully adjusted"
  marker_name <- paste0(
    "Fully adjusted + ",
    marker_label
  )
  
  hr_comparison <- hr_comparison %>%
    dplyr::transmute(
      sex,
      exposure,
      marker,
      quintile,
      
      estimate_fully_adjusted =
        .data[[paste0("estimate_", base_name)]],
      
      conf.low_fully_adjusted =
        .data[[paste0("conf.low_", base_name)]],
      
      conf.high_fully_adjusted =
        .data[[paste0("conf.high_", base_name)]],
      
      estimate_marker_adjusted =
        .data[[paste0("estimate_", marker_name)]],
      
      conf.low_marker_adjusted =
        .data[[paste0("conf.low_", marker_name)]],
      
      conf.high_marker_adjusted =
        .data[[paste0("conf.high_", marker_name)]],
      
      change_log_HR =
        log(estimate_marker_adjusted) -
        log(estimate_fully_adjusted)
    )
  
  list(
    hr_results = tidy,
    
    hr_comparison = hr_comparison,
    
    sample_summary = tibble::tibble(
      sex = sex_label,
      exposure = exposure_label,
      marker = marker_label,
      n = stats::nobs(m_base),
      events = m_base$nevent
    ),
    
    model_comparison = model_comparison
  )
}


run_all_pa_marker_sensitivities <- function(
    data,
    exposures,
    markers,
    sexes = c("Female", "Male"),
    base_confounders = default_tte_confounders(FALSE),
    sex_var = "sex_raw",
    age_entry_var = "age_entry_PA",
    age_exit_var = "age_exit_hip_PA",
    event_var = "event_hip_PA"
) {
  
  required_exposure_cols <- c(
    "exposure",
    "quintile_var"
  )
  
  required_marker_cols <- c(
    "marker",
    "marker_var"
  )
  
  if (!all(required_exposure_cols %in% names(exposures))) {
    stop(
      "`exposures` must contain: ",
      paste(required_exposure_cols, collapse = ", ")
    )
  }
  
  if (!all(required_marker_cols %in% names(markers))) {
    stop(
      "`markers` must contain: ",
      paste(required_marker_cols, collapse = ", ")
    )
  }
  
  analysis_grid <- tidyr::crossing(
    sex = sexes,
    exposures,
    markers
  )
  
  analyses <- analysis_grid %>%
    dplyr::mutate(
      result = purrr::pmap(
        list(
          sex,
          exposure,
          quintile_var,
          marker,
          marker_var
        ),
        function(
    sex,
    exposure,
    quintile_var,
    marker,
    marker_var
        ) {
          
          run_pa_marker_sensitivity(
            data = data,
            sex_label = sex,
            exposure_label = exposure,
            quintile_var = quintile_var,
            marker_label = marker,
            marker_var = marker_var,
            base_confounders = base_confounders,
            sex_var = sex_var,
            age_entry_var = age_entry_var,
            age_exit_var = age_exit_var,
            event_var = event_var
          )
        }
      )
    )
  
  list(
    hr_results = purrr::map_dfr(
      analyses$result,
      "hr_results"
    ) %>%
      dplyr::arrange(
        exposure,
        marker,
        sex,
        model,
        quintile
      ),
    
    hr_comparison = purrr::map_dfr(
      analyses$result,
      "hr_comparison"
    ) %>%
      dplyr::arrange(
        exposure,
        marker,
        sex,
        quintile
      ),
    
    sample_summary = purrr::map_dfr(
      analyses$result,
      "sample_summary"
    ) %>%
      dplyr::arrange(
        exposure,
        marker,
        sex
      ),
    
    model_comparison = purrr::map_dfr(
      analyses$result,
      "model_comparison"
    ) %>%
      dplyr::arrange(
        exposure,
        marker,
        sex
      )
  )
}