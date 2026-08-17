# =========================================================
# SURVIVAL ANALYSIS HELPERS
# =========================================================

# =========================================================
# SURVIVAL ANALYSIS HELPERS
# =========================================================

# =========================================================
# 1. Continuous Cox model helpers
# =========================================================



run_pa_whole_followup <- function(data,
                                  outcome_name,
                                  age_exit_var,
                                  event_var,
                                  make_ci_tables = FALSE) {
  
  stopifnot(age_exit_var %in% names(data))
  stopifnot(event_var %in% names(data))
  
  surv_txt <- paste0("Surv(age_entry_PA, ", age_exit_var, ", ", event_var, ")")
  
  f_unadj <- as.formula(
    paste0(surv_txt, " ~ log_MET_total")
  )
  
  f_min <- as.formula(
    paste0(
      surv_txt,
      " ~ log_MET_total + sex_raw + ethnicity_derived + height_clean + weight_clean"
    )
  )
  
  f_full_linear <- as.formula(
    paste0(
      surv_txt,
      " ~ log_MET_total + sex_raw + ethnicity_derived + height_clean + weight_clean + tdi_raw + education_level"
    )
  )
  
  f_full_spline <- as.formula(
    paste0(
      surv_txt,
      " ~ rms::rcs(log_MET_total, 4) + sex_raw + ethnicity_derived + height_clean + weight_clean + tdi_raw + education_level"
    )
  )
  
  f_spline_sexint <- as.formula(
    paste0(
      surv_txt,
      " ~ rms::rcs(log_MET_total, 4) * sex_raw + ethnicity_derived + height_clean + weight_clean + tdi_raw + education_level"
    )
  )
  
  m_unadj <- survival::coxph(f_unadj, data = data, model = TRUE, x = TRUE, y = TRUE)
  m_min <- survival::coxph(f_min, data = data, model = TRUE, x = TRUE, y = TRUE)
  m_full_linear <- survival::coxph(f_full_linear, data = data, model = TRUE, x = TRUE, y = TRUE)
  m_full_spline <- survival::coxph(f_full_spline, data = data, model = TRUE, x = TRUE, y = TRUE)
  m_spline_sexint <- survival::coxph(f_spline_sexint, data = data, model = TRUE, x = TRUE, y = TRUE)
  
  linear_vs_spline_lrt <- anova(m_full_linear, m_full_spline, test = "LRT")
  linear_vs_spline_aic <- AIC(m_full_linear, m_full_spline)
  
  sex_interaction_lrt <- anova(m_full_spline, m_spline_sexint, test = "LRT")
  sex_interaction_aic <- AIC(m_full_spline, m_spline_sexint)
  
  tidied_models <- NULL
  
  if (make_ci_tables) {
    tidied_models <- list(
      unadjusted = broom::tidy(m_unadj, exponentiate = TRUE, conf.int = TRUE),
      minimally_adjusted = broom::tidy(m_min, exponentiate = TRUE, conf.int = TRUE),
      fully_adjusted_linear = broom::tidy(m_full_linear, exponentiate = TRUE, conf.int = TRUE),
      fully_adjusted_spline = broom::tidy(m_full_spline, exponentiate = TRUE, conf.int = TRUE),
      spline_sex_interaction = broom::tidy(m_spline_sexint, exponentiate = TRUE, conf.int = TRUE)
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
      linear_vs_spline_lrt = linear_vs_spline_lrt,
      linear_vs_spline_aic = linear_vs_spline_aic,
      sex_interaction_lrt = sex_interaction_lrt,
      sex_interaction_aic = sex_interaction_aic
    )
  )
}


# SEX-SPECIFIC COX MODEL HELPER


run_pa_whole_followup_sex_specific <- function(data,
                                               outcome_name,
                                               age_exit_var,
                                               event_var,
                                               make_ci_tables = FALSE) {
  
  stopifnot(age_exit_var %in% names(data))
  stopifnot(event_var %in% names(data))
  
  surv_txt <- paste0("Surv(age_entry_PA, ", age_exit_var, ", ", event_var, ")")
  
  f_unadj <- as.formula(
    paste0(surv_txt, " ~ log_MET_total")
  )
  
  f_min <- as.formula(
    paste0(
      surv_txt,
      " ~ log_MET_total + ethnicity_derived + height_clean + weight_clean"
    )
  )
  
  f_full_linear <- as.formula(
    paste0(
      surv_txt,
      " ~ log_MET_total + ethnicity_derived + height_clean + weight_clean + tdi_raw + education_level"
    )
  )
  
  f_full_spline <- as.formula(
    paste0(
      surv_txt,
      " ~ rms::rcs(log_MET_total, 4) + ethnicity_derived + height_clean + weight_clean + tdi_raw + education_level"
    )
  )
  
  m_unadj <- survival::coxph(f_unadj, data = data, model = TRUE, x = TRUE, y = TRUE)
  m_min <- survival::coxph(f_min, data = data, model = TRUE, x = TRUE, y = TRUE)
  m_full_linear <- survival::coxph(f_full_linear, data = data, model = TRUE, x = TRUE, y = TRUE)
  m_full_spline <- survival::coxph(f_full_spline, data = data, model = TRUE, x = TRUE, y = TRUE)
  
  linear_vs_spline_lrt <- anova(m_full_linear, m_full_spline, test = "LRT")
  linear_vs_spline_aic <- AIC(m_full_linear, m_full_spline)
  
  tidied_models <- NULL
  
  if (make_ci_tables) {
    tidied_models <- list(
      unadjusted = broom::tidy(m_unadj, exponentiate = TRUE, conf.int = TRUE),
      minimally_adjusted = broom::tidy(m_min, exponentiate = TRUE, conf.int = TRUE),
      fully_adjusted_linear = broom::tidy(m_full_linear, exponentiate = TRUE, conf.int = TRUE),
      fully_adjusted_spline = broom::tidy(m_full_spline, exponentiate = TRUE, conf.int = TRUE)
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
      linear_vs_spline_lrt = linear_vs_spline_lrt,
      linear_vs_spline_aic = linear_vs_spline_aic
    )
  )
}


# =========================================================
# 2. QUINTILE COX MODEL HELPERS
# =========================================================

run_pa_quintile_analysis_sex_specific <- function(data,
                                                  outcome_name,
                                                  age_exit_var,
                                                  event_var,
                                                  model_specs = NULL,
                                                  make_ci_tables = FALSE) {
  
  stopifnot(age_exit_var %in% names(data))
  stopifnot(event_var %in% names(data))
  stopifnot("MET_total_quintile_f" %in% names(data))
  stopifnot("MET_total_quintile" %in% names(data))
  
  if (is.null(model_specs)) {
    model_specs <- list(
      "Unadjusted" = character(0),
      "Minimal adjustment" = c("height_clean", "weight_clean"),
      "Fully adjusted" = c("ethnicity_derived", "height_clean", "weight_clean", "tdi_raw", "education_level")
    )
  }
  
  surv_txt <- paste0("Surv(age_entry_PA, ", age_exit_var, ", ", event_var, ")")
  
  fit_model <- function(covars, use_trend = FALSE) {
    
    exposure_term <- if (use_trend) "MET_total_quintile" else "MET_total_quintile_f"
    rhs <- c(exposure_term, covars)
    
    form <- as.formula(
      paste0(surv_txt, " ~ ", paste(rhs, collapse = " + "))
    )
    
    survival::coxph(form, data = data, model = FALSE, x = FALSE, y = FALSE)
  }
  
  quintile_models <- lapply(model_specs, fit_model, use_trend = FALSE)
  trend_models <- lapply(model_specs, fit_model, use_trend = TRUE)
  
  quintile_event_counts <- table(
    data$MET_total_quintile_f,
    data[[event_var]],
    useNA = "ifany"
  )
  
  tidy_one_model <- function(model, model_name, outcome_name) {
    broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::filter(grepl("^MET_total_quintile_f", term)) %>%
      dplyr::mutate(
        quintile = sub("^MET_total_quintile_f", "", term),
        model = model_name,
        outcome = outcome_name
      ) %>%
      dplyr::select(outcome, model, quintile, term, estimate, conf.low, conf.high, p.value)
  }
  
  tidy_one_trend <- function(model, model_name, outcome_name) {
    broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::filter(term == "MET_total_quintile") %>%
      dplyr::mutate(
        model = model_name,
        outcome = outcome_name
      )
  }
  
  tidy_quintile <- NULL
  tidy_trend <- NULL
  
  if (make_ci_tables) {
    tidy_quintile <- dplyr::bind_rows(
      Map(tidy_one_model, quintile_models, names(quintile_models),
          MoreArgs = list(outcome_name = outcome_name))
    )
    
    tidy_trend <- dplyr::bind_rows(
      Map(tidy_one_trend, trend_models, names(trend_models),
          MoreArgs = list(outcome_name = outcome_name))
    )
  }
  
  list(
    outcome_name = outcome_name,
    model_specs = model_specs,
    models = list(
      quintile = quintile_models,
      quintile_trend = trend_models
    ),
    tidied_models = list(
      quintile = tidy_quintile,
      quintile_trend = tidy_trend
    ),
    quintile_event_counts = quintile_event_counts
  )
}

plot_pa_quintiles_models <- function(tidy_quint_df, outcome_name) {
  
  plot_df <- tidy_quint_df %>%
    dplyr::mutate(
      quintile = factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5")),
      model = factor(model, levels = c("Unadjusted", "Minimal adjustment", "Fully adjusted"))
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = quintile,
      y = estimate,
      colour = model,
      group = model
    )
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.5),
      size = 2.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      position = ggplot2::position_dodge(width = 0.5),
      width = 0.15
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "PA quintile",
      y = "Hazard ratio vs Q1",
      colour = "Model",
      title = paste0("Quintile model: ", outcome_name)
    ) +
    ggplot2::theme_minimal()
}

make_quintile_results_table <- function(res_quint) {
  
  res_quint$tidied_models$quintile %>%
    dplyr::mutate(
      quintile = factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5")),
      HR_CI = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high),
      p_value = formatC(p.value, format = "f", digits = 3)
    ) %>%
    dplyr::select(quintile, model, HR_CI, p_value) %>%
    tidyr::pivot_wider(
      names_from = model,
      values_from = c(HR_CI, p_value)
    ) %>%
    dplyr::arrange(quintile)
}

# Helper for DXA with pooled sex
run_pa_quintile_analysis <- function(data,
                                     outcome_name,
                                     age_exit_var,
                                     event_var,
                                     model_specs = NULL,
                                     make_ci_tables = FALSE) {
  
  stopifnot(age_exit_var %in% names(data))
  stopifnot(event_var %in% names(data))
  stopifnot("MET_total_quintile_f" %in% names(data))
  stopifnot("MET_total_quintile" %in% names(data))
  
  if (is.null(model_specs)) {
    model_specs <- list(
      "Unadjusted" = character(0),
      
      "Minimal adjustment" = c(
        "sex_raw",
        "height_clean",
        "weight_clean"
      ),
      
      "Fully adjusted" = c(
        "sex_raw",
        "ethnicity_derived",
        "height_clean",
        "weight_clean",
        "tdi_raw",
        "education_level"
      )
    )
  }
  
  surv_txt <- paste0(
    "Surv(age_entry_PA, ",
    age_exit_var,
    ", ",
    event_var,
    ")"
  )
  
  fit_model <- function(covars, use_trend = FALSE) {
    
    exposure_term <- if (use_trend) {
      "MET_total_quintile"
    } else {
      "MET_total_quintile_f"
    }
    
    rhs <- c(exposure_term, covars)
    
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
    fit_model,
    use_trend = FALSE
  )
  
  trend_models <- lapply(
    model_specs,
    fit_model,
    use_trend = TRUE
  )
  
  quintile_event_counts <- table(
    data$MET_total_quintile_f,
    data[[event_var]],
    useNA = "ifany"
  )
  
  tidy_one_model <- function(model, model_name, outcome_name) {
    
    broom::tidy(
      model,
      exponentiate = TRUE,
      conf.int = TRUE
    ) %>%
      dplyr::filter(
        grepl("^MET_total_quintile_f", term)
      ) %>%
      dplyr::mutate(
        quintile = sub(
          "^MET_total_quintile_f",
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
  
  tidy_one_trend <- function(model, model_name, outcome_name) {
    
    broom::tidy(
      model,
      exponentiate = TRUE,
      conf.int = TRUE
    ) %>%
      dplyr::filter(
        term == "MET_total_quintile"
      ) %>%
      dplyr::mutate(
        model = model_name,
        outcome = outcome_name
      )
  }
  
  tidy_quintile <- NULL
  tidy_trend <- NULL
  
  if (make_ci_tables) {
    
    tidy_quintile <- dplyr::bind_rows(
      Map(
        tidy_one_model,
        quintile_models,
        names(quintile_models),
        MoreArgs = list(
          outcome_name = outcome_name
        )
      )
    )
    
    tidy_trend <- dplyr::bind_rows(
      Map(
        tidy_one_trend,
        trend_models,
        names(trend_models),
        MoreArgs = list(
          outcome_name = outcome_name
        )
      )
    )
  }
  
  list(
    outcome_name = outcome_name,
    model_specs = model_specs,
    models = list(
      quintile = quintile_models,
      quintile_trend = trend_models
    ),
    tidied_models = list(
      quintile = tidy_quintile,
      quintile_trend = tidy_trend
    ),
    quintile_event_counts = quintile_event_counts
  )
}

run_pa_quintile_analysis_BMD <- function(data,
                                         outcome_name,
                                         age_exit_var,
                                         event_var,
                                         make_ci_tables = FALSE) {
  
  res <- run_pa_quintile_analysis(
    data = data,
    outcome_name = outcome_name,
    age_exit_var = age_exit_var,
    event_var = event_var,
    model_specs = model_specs_FN_BMD,
    make_ci_tables = make_ci_tables
  )
  
  list(
    outcome_name = outcome_name,
    model_specs = res$model_specs,
    models = res$models,
    tidied_models = res$tidied_models,
    quintile_event_counts = res$quintile_event_counts
  )
}

run_accel_quintile_analysis <- function(data,
                                        outcome_name,
                                        age_exit_var,
                                        event_var,
                                        make_ci_tables = FALSE) {
  
  stopifnot(age_exit_var %in% names(data))
  stopifnot(event_var %in% names(data))
  
  surv_txt <- paste0("Surv(age_entry_PA, ", age_exit_var, ", ", event_var, ")")
  
  f_quintile <- as.formula(
    paste0(
      surv_txt,
      " ~ MVPA_accel_quintile_f + sex_raw + ethnicity_derived + height_clean + weight_clean + tdi_raw + education_level"
    )
  )
  
  f_quintile_trend <- as.formula(
    paste0(
      surv_txt,
      " ~ MVPA_accel_quintile + sex_raw + ethnicity_derived + height_clean + weight_clean + tdi_raw + education_level"
    )
  )
  
  f_quintile_sexint <- as.formula(
    paste0(
      surv_txt,
      " ~ MVPA_accel_quintile_f * sex_raw + ethnicity_derived + height_clean + weight_clean + tdi_raw + education_level"
    )
  )
  
  m_quintile <- survival::coxph(f_quintile, data = data, model = FALSE, x = FALSE, y = FALSE)
  m_quintile_trend <- survival::coxph(f_quintile_trend, data = data, model = FALSE, x = FALSE, y = FALSE)
  m_quintile_sexint <- survival::coxph(f_quintile_sexint, data = data, model = FALSE, x = FALSE, y = FALSE)
  
  quintile_event_counts <- table(
    data$MVPA_accel_quintile_f,
    data[[event_var]],
    useNA = "ifany"
  )
  
  quintile_sexint_lrt <- anova(m_quintile, m_quintile_sexint, test = "LRT")
  quintile_sexint_aic <- AIC(m_quintile, m_quintile_sexint)
  
  tidied_models <- NULL
  
  if (make_ci_tables) {
    tidied_models <- list(
      quintile = broom::tidy(m_quintile, exponentiate = TRUE, conf.int = TRUE),
      quintile_trend = broom::tidy(m_quintile_trend, exponentiate = TRUE, conf.int = TRUE)
    )
  }
  
  list(
    outcome_name = outcome_name,
    formulas = list(
      quintile = f_quintile,
      quintile_trend = f_quintile_trend,
      quintile_sex_interaction = f_quintile_sexint
    ),
    models = list(
      quintile = m_quintile,
      quintile_trend = m_quintile_trend,
      quintile_sex_interaction = m_quintile_sexint
    ),
    tidied_models = tidied_models,
    comparisons = list(
      sex_interaction_lrt = quintile_sexint_lrt,
      sex_interaction_aic = quintile_sexint_aic
    ),
    quintile_event_counts = quintile_event_counts
  )
}

# =========================================================
# 3. Spline PLOTTING HELPERS
# =========================================================

add_hr_from_model <- function(model, newdat, data, x_var, ref_value = NULL) {
  
  if (is.null(ref_value)) {
    ref_value <- median(data[[x_var]], na.rm = TRUE)
  }
  
  pred <- predict(model, newdata = newdat, type = "lp", se.fit = TRUE)
  
  refdat <- newdat[1, , drop = FALSE]
  refdat[[x_var]] <- ref_value
  
  ref_lp <- predict(model, newdata = refdat, type = "lp")
  
  log_hr <- as.numeric(pred$fit - ref_lp)
  
  newdat$HR <- exp(log_hr)
  newdat$HR_low <- exp(log_hr - 1.96 * pred$se.fit)
  newdat$HR_high <- exp(log_hr + 1.96 * pred$se.fit)
  
  newdat
}

make_spline_newdata <- function(data,
                                x_var,
                                n_points = 200,
                                include_sex = TRUE,
                                sex_value = NULL) {
  
  newdat <- data.frame(
    x_tmp = seq(
      min(data[[x_var]], na.rm = TRUE),
      max(data[[x_var]], na.rm = TRUE),
      length.out = n_points
    ),
    ethnicity_derived = factor(
      rep(levels(data$ethnicity_derived)[1], n_points),
      levels = levels(data$ethnicity_derived)
    ),
    height_clean = mean(data$height_clean, na.rm = TRUE),
    weight_clean = mean(data$weight_clean, na.rm = TRUE),
    tdi_raw = mean(data$tdi_raw, na.rm = TRUE),
    education_level = factor(
      rep(levels(data$education_level)[1], n_points),
      levels = levels(data$education_level)
    )
  )
  
  names(newdat)[1] <- x_var
  
  if (include_sex) {
    if (is.null(sex_value)) {
      sex_value <- levels(data$sex_raw)[1]
    }
    
    newdat$sex_raw <- factor(
      rep(sex_value, n_points),
      levels = levels(data$sex_raw)
    )
  }
  
  newdat
}



plot_pa_spline <- function(model,
                           data,
                           outcome_name,
                           x_var = "log_MET_total",
                           x_label = "log(MET total + 1)",
                           n_points = 200) {
  
  if (is.null(model)) stop("model is NULL")
  if (!inherits(model, "coxph")) stop("model is not a coxph object")
  if (!x_var %in% names(data)) stop(paste("Variable not found in data:", x_var))
  
  newdat <- make_spline_newdata(
    data = data,
    x_var = x_var,
    n_points = n_points,
    include_sex = TRUE
  )
  
  newdat <- add_hr_from_model(
    model = model,
    newdat = newdat,
    data = data,
    x_var = x_var
  )
  
  ggplot2::ggplot(newdat, ggplot2::aes(x = .data[[x_var]], y = HR)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = HR_low, ymax = HR_high), alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = x_label,
      y = paste0("Hazard ratio for ", outcome_name),
      title = paste0("Restricted cubic spline model: ", outcome_name)
    ) +
    ggplot2::theme_minimal()
}


plot_pa_spline_sex_specific_model <- function(model,
                                              data,
                                              outcome_name,
                                              x_var = "log_MET_total",
                                              x_label = "log(MET total + 1)",
                                              n_points = 200) {
  
  if (is.null(model)) stop("model is NULL")
  if (!inherits(model, "coxph")) stop("model is not a coxph object")
  if (!x_var %in% names(data)) stop(paste("Variable not found in data:", x_var))
  
  newdat <- make_spline_newdata(
    data = data,
    x_var = x_var,
    n_points = n_points,
    include_sex = FALSE
  )
  
  newdat <- add_hr_from_model(
    model = model,
    newdat = newdat,
    data = data,
    x_var = x_var
  )
  
  ggplot2::ggplot(newdat, ggplot2::aes(x = .data[[x_var]], y = HR)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = HR_low, ymax = HR_high), alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = x_label,
      y = paste0("Hazard ratio for ", outcome_name),
      title = paste0("Restricted cubic spline model: ", outcome_name)
    ) +
    ggplot2::theme_minimal()
}

# -----------------------------
# composition
# -----------------------------

run_pa_quintile_sensitivity_sex_specific <- function(data,
                                                     sex_label,
                                                     outcome_name = "Fragility fracture",
                                                     age_exit_var = "age_exit_fragility_PA",
                                                     event_var = "event_fragility_PA",
                                                     age_entry_var = "age_entry_PA",
                                                     exposure_specs = NULL,
                                                     confounders = NULL) {
  
  stopifnot(age_entry_var %in% names(data))
  stopifnot(age_exit_var %in% names(data))
  stopifnot(event_var %in% names(data))
  
  if (is.null(confounders)) {
    confounders <- c(
      "ethnicity_derived",
      "height_clean",
      "weight_clean",
      "tdi_raw",
      "education_level"
    )
  }
  
  if (is.null(exposure_specs)) {
    exposure_specs <- list(
      "Total PA" = list(
        exposure = "cc_MET_total_trunc",
        adjust_for = character(0)
      ),
      "Vigorous PA" = list(
        exposure = "cc_MET_vig_trunc",
        adjust_for = character(0)
      ),
      "Moderate PA" = list(
        exposure = "cc_MET_mod_trunc",
        adjust_for = character(0)
      ),
      "Walking PA" = list(
        exposure = "cc_MET_walk_trunc",
        adjust_for = character(0)
      ),
      "Vigorous PA mutually adjusted" = list(
        exposure = "cc_MET_vig_trunc",
        adjust_for = c("cc_MET_mod_trunc", "cc_MET_walk_trunc")
      ),
      "Moderate PA mutually adjusted" = list(
        exposure = "cc_MET_mod_trunc",
        adjust_for = c("cc_MET_vig_trunc", "cc_MET_walk_trunc")
      ),
      "Walking PA mutually adjusted" = list(
        exposure = "cc_MET_walk_trunc",
        adjust_for = c("cc_MET_vig_trunc", "cc_MET_mod_trunc")
      )
    )
  }
  
  results <- list()
  idx <- 1
  
  for (exposure_name in names(exposure_specs)) {
    
    exposure_var <- exposure_specs[[exposure_name]]$exposure
    adjust_for <- exposure_specs[[exposure_name]]$adjust_for
    
    needed_vars <- unique(c(
      age_entry_var,
      age_exit_var,
      event_var,
      exposure_var,
      adjust_for,
      confounders
    ))
    
    missing_vars <- setdiff(needed_vars, names(data))
    
    if (length(missing_vars) > 0) {
      warning(
        "Skipping ", exposure_name,
        " because these variables are missing: ",
        paste(missing_vars, collapse = ", ")
      )
      next
    }
    
    dat <- data %>%
      dplyr::select(dplyr::all_of(needed_vars)) %>%
      tidyr::drop_na()
    
    if (nrow(dat) == 0) next
    if (dplyr::n_distinct(dat[[event_var]]) < 2) next
    
    dat <- dat %>%
      dplyr::mutate(
        exposure_quintile = dplyr::ntile(.data[[exposure_var]], 5),
        exposure_quintile_f = factor(
          paste0("Q", exposure_quintile),
          levels = paste0("Q", 1:5)
        )
      )
    
    surv_txt <- paste0(
      "Surv(",
      age_entry_var, ", ",
      age_exit_var, ", ",
      event_var,
      ")"
    )
    
    rhs_quintile <- c("exposure_quintile_f", confounders, adjust_for)
    rhs_trend <- c("exposure_quintile", confounders, adjust_for)
    
    f_quintile <- stats::as.formula(
      paste0(surv_txt, " ~ ", paste(rhs_quintile, collapse = " + "))
    )
    
    f_trend <- stats::as.formula(
      paste0(surv_txt, " ~ ", paste(rhs_trend, collapse = " + "))
    )
    
    m_quintile <- survival::coxph(f_quintile, data = dat)
    m_trend <- survival::coxph(f_trend, data = dat)
    
    trend_p <- broom::tidy(m_trend) %>%
      dplyr::filter(term == "exposure_quintile") %>%
      dplyr::pull(p.value)
    
    res <- broom::tidy(
      m_quintile,
      exponentiate = TRUE,
      conf.int = TRUE
    ) %>%
      dplyr::filter(grepl("^exposure_quintile_f", term)) %>%
      dplyr::mutate(
        outcome = outcome_name,
        sex = sex_label,
        exposure = exposure_name,
        quintile = sub("^exposure_quintile_f", "", term),
        HR = estimate,
        CI_lower = conf.low,
        CI_upper = conf.high,
        p_trend = trend_p
      ) %>%
      dplyr::select(
        outcome, sex, exposure, quintile,
        HR, CI_lower, CI_upper,
        p.value, p_trend
      )
    
    results[[idx]] <- res
    idx <- idx + 1
  }
  
  dplyr::bind_rows(results)
}

# Code adjusted to filter for only the variables below, may need to change filtering if want other 
# variables to show in forest

plot_pa_quintile_sensitivity <- function(plot_df, outcome_name = "Fragility fracture") {
  
  exposure_levels <- c(
    "Total PA",
    "Walking PA mutually adjusted",
    "Moderate PA mutually adjusted",
    "Vigorous PA mutually adjusted"
  )
  
  activity_cols <- c(
    "Total PA" = "#000000",
    "Vigorous PA mutually adjusted" = "#E69F00",  # orange
    "Moderate PA mutually adjusted" = "#0072B2",  # blue
    "Walking PA mutually adjusted" = "#CC79A7"    # purple
  )
  
  plot_df <- plot_df %>%
    dplyr::filter(exposure %in% exposure_levels) %>%
    dplyr::mutate(
      quintile = factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5")),
      exposure = factor(exposure, levels = exposure_levels)
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = quintile,
      y = HR,
      colour = exposure,
      group = exposure
    )
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.6),
      size = 2.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.15,
      linewidth = 0.7
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = activity_cols) +
    ggplot2::facet_wrap(~ sex, ncol = 1) +
    ggplot2::labs(
      title = paste0("Activity type quintile sensitivity analysis: ", outcome_name),
      x = "Activity quintile",
      y = "Hazard ratio vs Q1",
      colour = "Model"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

# =========================================================
# Helper: fracture site quintile sensitivity analysis
# =========================================================

run_fracture_site_quintile_sensitivity <- function(data,
                                                   outcome_name,
                                                   age_exit_var,
                                                   event_var,
                                                   file_stub,
                                                   tables_dir = TABLES_DIR,
                                                   figures_dir = TTE_FIGURES_DIR) {
  
  cat("\n========================================\n")
  cat("Fracture site quintile sensitivity:", outcome_name, "\n")
  cat("========================================\n")
  
  female_dat_site <- data %>%
    dplyr::filter(sex_raw == levels(data$sex_raw)[1]) %>%
    droplevels()
  
  male_dat_site <- data %>%
    dplyr::filter(sex_raw == levels(data$sex_raw)[2]) %>%
    droplevels()
  
  res_site <- dplyr::bind_rows(
    run_pa_quintile_sensitivity_sex_specific(
      data = female_dat_site,
      sex_label = levels(data$sex_raw)[1],
      outcome_name = outcome_name,
      age_exit_var = age_exit_var,
      event_var = event_var
    ),
    run_pa_quintile_sensitivity_sex_specific(
      data = male_dat_site,
      sex_label = levels(data$sex_raw)[2],
      outcome_name = outcome_name,
      age_exit_var = age_exit_var,
      event_var = event_var
    )
  ) %>%
    dplyr::mutate(
      HR = round(HR, 4),
      CI_lower = round(CI_lower, 4),
      CI_upper = round(CI_upper, 4),
      p.value = round(p.value, 4),
      p_trend = round(p_trend, 4)
    )
  
  write.csv(
    res_site,
    file.path(tables_dir, paste0("Results_survival_quintile_activity_sensitivity_", file_stub, ".csv")),
    row.names = FALSE
  )
  
  table_site <- res_site %>%
    dplyr::mutate(
      HR_CI = sprintf("%.2f (%.2f to %.2f)", HR, CI_lower, CI_upper),
      p_value = formatC(p.value, format = "f", digits = 3),
      p_trend = formatC(p_trend, format = "f", digits = 3)
    ) %>%
    dplyr::select(
      outcome, sex, exposure, quintile,
      HR_CI, p_value, p_trend
    ) %>%
    dplyr::arrange(sex, exposure, quintile)
  
  write.csv(
    table_site,
    file.path(tables_dir, paste0("Table_survival_quintile_activity_sensitivity_", file_stub, ".csv")),
    row.names = FALSE
  )
  
  p_site <- plot_pa_quintile_sensitivity(
    res_site,
    outcome_name = outcome_name
  )
  
  print(p_site)
  
  ggplot2::ggsave(
    filename = file.path(figures_dir, paste0("Fig_survival_quintile_activity_sensitivity_", file_stub, ".png")),
    plot = p_site,
    width = 14,
    height = 6,
    dpi = 300
  )
  
  list(
    outcome_name = outcome_name,
    age_exit_var = age_exit_var,
    event_var = event_var,
    file_stub = file_stub,
    results = res_site,
    table = table_site,
    plot = p_site
  )
}


plot_pa_quintile_sensitivity_selected <- function(plot_df, outcome_name) {
  
  selected_exposures <- c(
    "Total PA",
    "Vigorous PA",
    "Vigorous PA adjusted for moderate",
    "Vigorous PA adjusted for walking"
  )
  
  plot_df <- plot_df %>%
    dplyr::filter(exposure %in% selected_exposures) %>%
    dplyr::mutate(
      quintile = factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5")),
      exposure = factor(exposure, levels = selected_exposures)
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = HR,
      y = quintile,
      colour = sex
    )
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = CI_lower, xmax = CI_upper),
      height = 0.15,
      linewidth = 0.7
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_grid(sex ~ exposure) +
    ggplot2::labs(
      title = paste0("Selected activity-type quintile models: ", outcome_name),
      x = "Hazard ratio vs Q1",
      y = "Activity quintile",
      colour = "Sex"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
}



# ---------------------------------------------------------
#  Fine-Gray helper
# ---------------------------------------------------------

run_crr_fragility_sex_specific <- function(data, sex_label) {
  
  data <- data %>%
    dplyr::filter(
      !is.na(fg_time_fragility),
      fg_time_fragility >= 0,
      !is.na(event_fragility_fg),
      !is.na(log_MET_total),
      !is.na(ethnicity_derived),
      !is.na(height_clean),
      !is.na(weight_clean),
      !is.na(tdi_raw),
      !is.na(education_level)
    ) %>%
    droplevels()
  
  cov_mat <- model.matrix(
    ~ splines::ns(log_MET_total, df = 3) +
      ethnicity_derived +
      height_clean +
      weight_clean +
      tdi_raw +
      education_level,
    data = data
  )[, -1, drop = FALSE]
  
  fg_model <- cmprsk::crr(
    ftime = data$fg_time_fragility,
    fstatus = data$event_fragility_fg,
    cov1 = cov_mat,
    failcode = 1,
    cencode = 0
  )
  
  tidy <- data.frame(
    sex = sex_label,
    term = names(fg_model$coef),
    estimate = exp(fg_model$coef),
    conf.low = exp(fg_model$coef - 1.96 * sqrt(diag(fg_model$var))),
    conf.high = exp(fg_model$coef + 1.96 * sqrt(diag(fg_model$var))),
    p.value = 2 * pnorm(abs(fg_model$coef / sqrt(diag(fg_model$var))), lower.tail = FALSE),
    row.names = NULL
  )
  
  list(
    sex = sex_label,
    model = fg_model,
    tidy = tidy
  )
}


#### BMD cohort helper ###

run_pa_quintile_analysis_sex_specific_BMD <- function(data,
                                                      outcome_name,
                                                      age_exit_var,
                                                      event_var,
                                                      make_ci_tables = FALSE) {
  
  stopifnot("FN_BMD_mean" %in% names(data))
  
  model_specs_BMD <- list(
    "Unadjusted" = character(0),
    "Minimal adjustment" = c("height_clean", "weight_clean"),
    "Fully adjusted + FN BMD" = c(
      "ethnicity_derived",
      "height_clean",
      "weight_clean",
      "tdi_raw",
      "education_level",
      "FN_BMD_mean"
    )
  )
  
  run_pa_quintile_analysis_sex_specific(
    data = data,
    outcome_name = outcome_name,
    age_exit_var = age_exit_var,
    event_var = event_var,
    model_specs = model_specs_BMD,
    make_ci_tables = make_ci_tables
  )
}

#### falls ####

run_pa_quintile_analysis_sex_specific_falls <- function(data,
                                                        outcome_name,
                                                        age_exit_var,
                                                        event_var,
                                                        make_ci_tables = FALSE) {
  
  data <- data %>%
    dplyr::mutate(
      falls_binary = dplyr::case_when(
        self_reported_falls_clean == "No falls" ~ "No falls",
        self_reported_falls_clean %in% c(
          "Only one fall",
          "More than one fall"
        ) ~ "One or more falls",
        TRUE ~ NA_character_
      ),
      falls_binary = factor(
        falls_binary,
        levels = c("No falls", "One or more falls")
      )
    ) %>%
    dplyr::filter(!is.na(falls_binary))
  
  model_specs_falls <- list(
    "Fully adjusted" = c(
      "ethnicity_derived",
      "height_clean",
      "weight_clean",
      "tdi_raw",
      "education_level"
    ),
    "Fully adjusted + falls" = c(
      "ethnicity_derived",
      "height_clean",
      "weight_clean",
      "tdi_raw",
      "education_level",
      "falls_binary"
    )
  )
  
  run_pa_quintile_analysis_sex_specific(
    data = data,
    outcome_name = outcome_name,
    age_exit_var = age_exit_var,
    event_var = event_var,
    model_specs = model_specs_falls,
    make_ci_tables = make_ci_tables
  )
}

### Grip helper ###

run_pa_quintile_analysis_sex_specific_grip_light <- function(data,
                                                             outcome_name,
                                                             age_exit_var,
                                                             event_var) {
  
  data <- data %>%
    dplyr::filter(!is.na(grip_strength_mean)) %>%
    droplevels()
  
  model_specs_grip <- list(
    "Fully adjusted" = c(
      "ethnicity_derived", "height_clean", "weight_clean",
      "tdi_raw", "education_level"
    ),
    "Fully adjusted + grip" = c(
      "ethnicity_derived", "height_clean", "weight_clean",
      "tdi_raw", "education_level", "grip_strength_mean"
    )
  )
  
  res <- run_pa_quintile_analysis_sex_specific(
    data = data,
    outcome_name = outcome_name,
    age_exit_var = age_exit_var,
    event_var = event_var,
    model_specs = model_specs_grip,
    make_ci_tables = TRUE
  )
  
  list(
    outcome_name = outcome_name,
    tidied_models = res$tidied_models,
    quintile_event_counts = res$quintile_event_counts
  )
}

### Heel BMD helper ###

run_pa_quintile_analysis_sex_specific_heel_BMD <- function(data,
                                                           outcome_name,
                                                           age_exit_var,
                                                           event_var,
                                                           make_ci_tables = TRUE
) {
  
  data <- data %>%
    dplyr::filter(!is.na(heel_BMD)) %>%
    droplevels()
  
  model_specs_heel_BMD <- list(
    "Fully adjusted" = c(
      "ethnicity_derived", "height_clean", "weight_clean",
      "tdi_raw", "education_level"
    ),
    "Fully adjusted + heel BMD" = c(
      "ethnicity_derived", "height_clean", "weight_clean",
      "tdi_raw", "education_level", "heel_BMD"
    )
  )
  
  res <- run_pa_quintile_analysis_sex_specific(
    data = data,
    outcome_name = outcome_name,
    age_exit_var = age_exit_var,
    event_var = event_var,
    model_specs = model_specs_heel_BMD,
    make_ci_tables = TRUE
  )
  
  list(
    outcome_name = outcome_name,
    tidied_models = res$tidied_models,
    quintile_event_counts = res$quintile_event_counts
  )
}

run_crr_sex_specific <- function(data, sex_label, time_var, event_var, nk = 4) {
  
  data <- data %>%
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
    data$log_MET_total,
    nk = nk,
    inclx = TRUE
  )
  
  colnames(rcs_mat) <- paste0("rcs_log_MET_total_", seq_len(ncol(rcs_mat)))
  
  cov_data <- dplyr::bind_cols(
    as.data.frame(rcs_mat),
    data %>%
      dplyr::select(
        ethnicity_derived,
        height_clean,
        weight_clean,
        tdi_raw,
        education_level
      )
  )
  
  cov_mat <- model.matrix(~ ., data = cov_data)[, -1, drop = FALSE]
  
  fg_model <- cmprsk::crr(
    ftime = data[[time_var]],
    fstatus = data[[event_var]],
    cov1 = cov_mat,
    failcode = 1,
    cencode = 0
  )
  
  tidy <- data.frame(
    sex = sex_label,
    term = names(fg_model$coef),
    estimate = exp(fg_model$coef),
    conf.low = exp(fg_model$coef - 1.96 * sqrt(diag(fg_model$var))),
    conf.high = exp(fg_model$coef + 1.96 * sqrt(diag(fg_model$var))),
    p.value = 2 * pnorm(
      abs(fg_model$coef / sqrt(diag(fg_model$var))),
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
    cov_names = colnames(cov_mat)
  )
}


### ALMI helper ###

run_pa_quintile_analysis_sex_specific_ALMI <- function(data,
                                                       outcome_name,
                                                       age_exit_var,
                                                       event_var,
                                                       make_ci_tables = TRUE
) {
  
  data <- data %>%
    dplyr::filter(!is.na(heel_BMD)) %>%
    droplevels()
  
  model_specs_heel_BMD <- list(
    "Fully adjusted" = c(
      "ethnicity_derived", "height_clean", "weight_clean",
      "tdi_raw", "education_level"
    ),
    "Fully adjusted + ALMI" = c(
      "ethnicity_derived", "height_clean", "weight_clean",
      "tdi_raw", "education_level", "appendicular_LM_index"
    )
  )
  
  res <- run_pa_quintile_analysis_sex_specific(
    data = data,
    outcome_name = outcome_name,
    age_exit_var = age_exit_var,
    event_var = event_var,
    model_specs = model_specs_heel_BMD,
    make_ci_tables = TRUE
  )
  
  list(
    outcome_name = outcome_name,
    tidied_models = res$tidied_models,
    quintile_event_counts = res$quintile_event_counts
  )
}



run_crr_sex_specific <- function(data, sex_label, time_var, event_var, nk = 4) {
  
  data <- data %>%
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
    data$log_MET_total,
    nk = nk,
    inclx = TRUE
  )
  
  colnames(rcs_mat) <- paste0("rcs_log_MET_total_", seq_len(ncol(rcs_mat)))
  
  cov_data <- dplyr::bind_cols(
    as.data.frame(rcs_mat),
    data %>%
      dplyr::select(
        ethnicity_derived,
        height_clean,
        weight_clean,
        tdi_raw,
        education_level
      )
  )
  
  cov_mat <- model.matrix(~ ., data = cov_data)[, -1, drop = FALSE]
  
  fg_model <- cmprsk::crr(
    ftime = data[[time_var]],
    fstatus = data[[event_var]],
    cov1 = cov_mat,
    failcode = 1,
    cencode = 0
  )
  
  tidy <- data.frame(
    sex = sex_label,
    term = names(fg_model$coef),
    estimate = exp(fg_model$coef),
    conf.low = exp(fg_model$coef - 1.96 * sqrt(diag(fg_model$var))),
    conf.high = exp(fg_model$coef + 1.96 * sqrt(diag(fg_model$var))),
    p.value = 2 * pnorm(
      abs(fg_model$coef / sqrt(diag(fg_model$var))),
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
    cov_names = colnames(cov_mat)
  )
}

# =========================================================
# Model specifications: walking pace sensitivity analysis
# Sex-specific models
# =========================================================

model_specs_walking_pace <- list(
  
  "Fully adjusted" = c(
    "ethnicity_derived",
    "height_clean",
    "weight_clean",
    "tdi_raw",
    "education_level"
  ),
  
  "Fully adjusted + walking pace" = c(
    "ethnicity_derived",
    "height_clean",
    "weight_clean",
    "tdi_raw",
    "education_level",
    "usual_walking_pace_raw"
  )
)

# =========================================================
# PA quintile sensitivity analysis: walking pace
# =========================================================

run_pa_quintile_analysis_sex_specific_walking_pace <- function(
    data,
    outcome_name,
    age_exit_var,
    event_var,
    make_ci_tables = TRUE
) {
  
  run_pa_quintile_analysis(
    data = data,
    outcome_name = outcome_name,
    age_exit_var = age_exit_var,
    event_var = event_var,
    model_specs = model_specs_walking_pace,
    make_ci_tables = make_ci_tables
  )
}

# =========================================================
# Run a sensitivity analysis separately in women and men
# =========================================================

run_sensitivity_by_sex <- function(
    data,
    analysis_fun,
    sensitivity_var,
    outcome_label,
    age_exit_var,
    event_var,
    make_ci_tables = TRUE,
    pass_make_ci_tables = TRUE
) {
  
  required_vars <- c(
    "sex_raw",
    sensitivity_var,
    age_exit_var,
    event_var
  )
  
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      "The following required variables are missing: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  if (!is.function(analysis_fun)) {
    stop("`analysis_fun` must be a function.")
  }
  
  analysis_data <- data %>%
    dplyr::filter(
      !is.na(.data[[sensitivity_var]]),
      !is.na(.data[[age_exit_var]]),
      !is.na(.data[[event_var]])
    )
  
  female_data <- analysis_data %>%
    dplyr::filter(sex_raw == "Female") %>%
    droplevels()
  
  male_data <- analysis_data %>%
    dplyr::filter(sex_raw == "Male") %>%
    droplevels()
  
  if (nrow(female_data) == 0) {
    stop("No female observations remain for: ", sensitivity_var)
  }
  
  if (nrow(male_data) == 0) {
    stop("No male observations remain for: ", sensitivity_var)
  }
  
  run_one_sex <- function(dat, sex_label) {
    
    model_args <- list(
      data = dat,
      outcome_name = paste0(outcome_label, " - ", sex_label),
      age_exit_var = age_exit_var,
      event_var = event_var
    )
    
    if (pass_make_ci_tables) {
      model_args$make_ci_tables <- make_ci_tables
    }
    
    do.call(
      analysis_fun,
      model_args
    )
  }
  
  list(
    female = run_one_sex(female_data, "Female"),
    male = run_one_sex(male_data, "Male"),
    sample_summary = tibble::tibble(
      sex = c("Female", "Male"),
      n = c(nrow(female_data), nrow(male_data)),
      events = c(
        sum(female_data[[event_var]] == 1, na.rm = TRUE),
        sum(male_data[[event_var]] == 1, na.rm = TRUE)
      )
    )
  )
}

# =========================================================
# Save sex-specific sensitivity results
# =========================================================

save_sensitivity_results <- function(
    results,
    filename,
    output_dir = DATA_DERIVED
) {
  
  saveRDS(
    results,
    file.path(output_dir, filename)
  )
  
  invisible(results)
}

# =========================================================
# Extract sensitivity-analysis results for supplementary table
# =========================================================

extract_sensitivity_model_row <- function(
    tidy_quintile,
    model_name,
    sex_label,
    marker_label,
    analysis_label,
    n,
    events
) {
  
  hr_values <- tidy_quintile %>%
    dplyr::filter(
      model == model_name,
      quintile %in% c("Q2", "Q3", "Q4", "Q5")
    ) %>%
    dplyr::mutate(
      HR_95CI = sprintf(
        "%.2f (%.2f–%.2f)",
        estimate,
        conf.low,
        conf.high
      )
    ) %>%
    dplyr::select(
      quintile,
      HR_95CI
    ) %>%
    tidyr::pivot_wider(
      names_from = quintile,
      values_from = HR_95CI
    )
  
  dplyr::bind_cols(
    tibble::tibble(
      Sex = sex_label,
      Marker = marker_label,
      Analysis = analysis_label,
      N = n,
      Events = events
    ),
    hr_values
  )
}

# =========================================================
# Extract one complete sensitivity comparison
# =========================================================

extract_sensitivity_comparison <- function(
    sensitivity_result,
    main_result_female,
    main_result_male,
    marker_label,
    additional_model_name
) {
  
  extract_for_sex <- function(
    sex_name,
    sex_label,
    main_result
  ) {
    
    sensitivity_sex <- sensitivity_result[[sex_name]]
    
    tidy_main <- main_result$tidied_models$quintile
    tidy_sensitivity <- sensitivity_sex$tidied_models$quintile
    
    # Main fully adjusted model
    main_model <- main_result$models$quintile[["Fully adjusted"]]
    
    # Marker-complete sample counts
    sample_row <- sensitivity_result$sample_summary %>%
      dplyr::filter(sex == sex_label)
    
    if (nrow(sample_row) != 1) {
      stop(
        "Could not identify one sample-summary row for ",
        sex_label,
        " and ",
        marker_label
      )
    }
    
    dplyr::bind_rows(
      
      extract_sensitivity_model_row(
        tidy_quintile = tidy_main,
        model_name = "Fully adjusted",
        sex_label = sex_label,
        marker_label = marker_label,
        analysis_label = "Main fully adjusted model",
        n = main_model$n,
        events = main_model$nevent
      ),
      
      extract_sensitivity_model_row(
        tidy_quintile = tidy_sensitivity,
        model_name = "Fully adjusted",
        sex_label = sex_label,
        marker_label = marker_label,
        analysis_label =
          "Fully adjusted, marker-complete sample",
        n = sample_row$n,
        events = sample_row$events
      ),
      
      extract_sensitivity_model_row(
        tidy_quintile = tidy_sensitivity,
        model_name = additional_model_name,
        sex_label = sex_label,
        marker_label = marker_label,
        analysis_label = additional_model_name,
        n = sample_row$n,
        events = sample_row$events
      )
    )
  }
  
  dplyr::bind_rows(
    extract_for_sex(
      sex_name = "female",
      sex_label = "Female",
      main_result = main_result_female
    ),
    extract_for_sex(
      sex_name = "male",
      sex_label = "Male",
      main_result = main_result_male
    )
  )
}

#### helper to compare Total PA and Vig with adjustment for walking speed to see different predictive ability for hip fracture using AIC
# =============================================================================
# Walking-pace sensitivity analysis helper
# =============================================================================

#' Compare Cox models with and without adjustment for walking pace
#'
#' Fits two nested Cox models on an identical complete-case sample:
#'   1. Fully adjusted model
#'   2. Fully adjusted model additionally adjusted for walking pace
#'
#' The function returns the fitted models, likelihood ratio test, AIC
#' comparison, and physical-activity quintile hazard ratios.
#'
#' @param data Analysis dataset.
#' @param sex_label Sex to analyse, matching the values in sex_var.
#' @param exposure_var Continuous PA variable from which quintiles are derived.
#' @param exposure_label Descriptive label used in output tables.
#' @param sex_var Name of the sex variable.
#' @param age_entry_var Start-time variable.
#' @param age_exit_var End-time variable.
#' @param event_var Event indicator.
#' @param confounders Variables included in the fully adjusted model.
#' @param walking_pace_var Original walking-pace variable.
#'
#' @return A list containing the analysis data, fitted models, LRT,
#'   AIC comparison and HR comparison.
#'
compare_walking_pace_models <- function(
    data,
    sex_label,
    exposure_var,
    exposure_label,
    sex_var = "sex_raw",
    age_entry_var = "age_entry_PA",
    age_exit_var = "age_exit_hip_PA",
    event_var = "event_hip_PA",
    confounders = c(
      "ethnicity_derived",
      "height_clean",
      "weight_clean",
      "tdi_raw",
      "education_level"
    ),
    walking_pace_var = "usual_walking_pace_raw"
) {
  
  # ---------------------------------------------------------------------------
  # Check required variables
  # ---------------------------------------------------------------------------
  
  required_vars <- unique(c(
    sex_var,
    age_entry_var,
    age_exit_var,
    event_var,
    exposure_var,
    confounders,
    walking_pace_var
  ))
  
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      "Variables missing from the analysis dataset: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  # ---------------------------------------------------------------------------
  # Create common complete-case dataset
  # ---------------------------------------------------------------------------
  
  analysis_dat <- data %>%
    dplyr::filter(
      .data[[sex_var]] == sex_label
    ) %>%
    dplyr::select(
      dplyr::all_of(required_vars)
    ) %>%
    dplyr::mutate(
      walking_pace = dplyr::case_when(
        .data[[walking_pace_var]] == "Slow pace" ~
          "Slow pace",
        
        .data[[walking_pace_var]] == "Steady average pace" ~
          "Steady average pace",
        
        .data[[walking_pace_var]] == "Brisk pace" ~
          "Brisk pace",
        
        TRUE ~ NA_character_
      ),
      
      walking_pace = factor(
        walking_pace,
        levels = c(
          "Slow pace",
          "Steady average pace",
          "Brisk pace"
        )
      )
    ) %>%
    tidyr::drop_na() %>%
    droplevels() %>%
    dplyr::mutate(
      exposure_quintile = dplyr::ntile(
        .data[[exposure_var]],
        5
      ),
      
      exposure_quintile_f = factor(
        paste0("Q", exposure_quintile),
        levels = paste0("Q", 1:5)
      )
    )
  
  if (nrow(analysis_dat) == 0) {
    stop("No complete observations remained for ", sex_label, ".")
  }
  
  if (dplyr::n_distinct(analysis_dat[[event_var]]) < 2) {
    stop("The analysis dataset does not contain both events and non-events.")
  }
  
  # ---------------------------------------------------------------------------
  # Construct model formulas
  # ---------------------------------------------------------------------------
  
  surv_term <- paste0(
    "survival::Surv(",
    age_entry_var, ", ",
    age_exit_var, ", ",
    event_var,
    ")"
  )
  
  base_rhs <- c(
    "exposure_quintile_f",
    confounders
  )
  
  walking_rhs <- c(
    base_rhs,
    "walking_pace"
  )
  
  formula_base <- stats::as.formula(
    paste(
      surv_term,
      "~",
      paste(base_rhs, collapse = " + ")
    )
  )
  
  formula_walking <- stats::as.formula(
    paste(
      surv_term,
      "~",
      paste(walking_rhs, collapse = " + ")
    )
  )
  
  # ---------------------------------------------------------------------------
  # Fit nested Cox models
  # ---------------------------------------------------------------------------
  
  model_base <- survival::coxph(
    formula = formula_base,
    data = analysis_dat,
    ties = "efron",
    x = TRUE,
    model = TRUE
  )
  
  model_walking <- survival::coxph(
    formula = formula_walking,
    data = analysis_dat,
    ties = "efron",
    x = TRUE,
    model = TRUE
  )
  
  if (stats::nobs(model_base) != stats::nobs(model_walking)) {
    stop(
      "The nested models used different numbers of observations."
    )
  }
  
  # ---------------------------------------------------------------------------
  # Formal model comparison
  # ---------------------------------------------------------------------------
  
  lrt <- stats::anova(
    model_base,
    model_walking,
    test = "LRT"
  )
  
  aic_results <- stats::AIC(
    model_base,
    model_walking
  )
  
  model_comparison <- tibble::tibble(
    sex = sex_label,
    exposure = exposure_label,
    N = stats::nobs(model_base),
    events = model_base$nevent,
    
    AIC_fully_adjusted =
      unname(aic_results["model_base", "AIC"]),
    
    AIC_walking_pace =
      unname(aic_results["model_walking", "AIC"]),
    
    # Negative value favours the walking-pace-adjusted model.
    delta_AIC =
      AIC_walking_pace - AIC_fully_adjusted,
    
    LRT_chi_square = unname(lrt$Chisq[2]),
    LRT_df = unname(lrt$Df[2]),
    LRT_p_value = unname(lrt$`Pr(>|Chi|)`[2])
  )
  
  # ---------------------------------------------------------------------------
  # Extract PA quintile estimates
  # ---------------------------------------------------------------------------
  
  extract_hr <- function(model, model_label) {
    
    broom::tidy(
      model,
      exponentiate = TRUE,
      conf.int = TRUE
    ) %>%
      dplyr::filter(
        grepl("^exposure_quintile_f", term)
      ) %>%
      dplyr::transmute(
        sex = sex_label,
        exposure = exposure_label,
        model = model_label,
        
        quintile = sub(
          "^exposure_quintile_f",
          "",
          term
        ),
        
        estimate = estimate,
        conf.low = conf.low,
        conf.high = conf.high,
        p.value = p.value
      )
  }
  
  hr_results <- dplyr::bind_rows(
    extract_hr(
      model_base,
      "Fully adjusted"
    ),
    
    extract_hr(
      model_walking,
      "Fully adjusted + walking pace"
    )
  ) %>%
    dplyr::mutate(
      quintile = factor(
        quintile,
        levels = c("Q2", "Q3", "Q4", "Q5")
      ),
      
      model = factor(
        model,
        levels = c(
          "Fully adjusted",
          "Fully adjusted + walking pace"
        )
      )
    )
  
  # ---------------------------------------------------------------------------
  # Side-by-side estimate comparison
  # ---------------------------------------------------------------------------
  
  hr_comparison <- hr_results %>%
    dplyr::select(
      sex,
      exposure,
      quintile,
      model,
      estimate,
      conf.low,
      conf.high,
      p.value
    ) %>%
    tidyr::pivot_wider(
      names_from = model,
      values_from = c(
        estimate,
        conf.low,
        conf.high,
        p.value
      ),
      names_glue = "{.value}_{model}"
    ) %>%
    dplyr::rename_with(
      ~ gsub(
        "Fully adjusted \\+ walking pace",
        "walking_pace",
        .x
      )
    ) %>%
    dplyr::rename_with(
      ~ gsub(
        "Fully adjusted",
        "fully_adjusted",
        .x
      )
    ) %>%
    dplyr::mutate(
      change_log_HR =
        log(estimate_walking_pace) -
        log(estimate_fully_adjusted)
    )
  
  # ---------------------------------------------------------------------------
  # Return complete result object
  # ---------------------------------------------------------------------------
  
  list(
    metadata = list(
      sex = sex_label,
      exposure = exposure_label,
      exposure_var = exposure_var,
      N = nrow(analysis_dat),
      events = sum(analysis_dat[[event_var]])
    ),
    
    models = list(
      fully_adjusted = model_base,
      walking_pace = model_walking
    ),
    
    likelihood_ratio_test = lrt,
    model_comparison = model_comparison,
    hr_results = hr_results,
    hr_comparison = hr_comparison
  )
}