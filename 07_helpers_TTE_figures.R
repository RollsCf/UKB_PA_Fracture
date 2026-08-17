# =============================================================================
# 07_helpers_figures_survival_clean.R
# Plotting and table helpers for hip fracture survival analyses
# =============================================================================

# =============================================================================
# 1. Colour schemes
# =============================================================================

activity_cols <- c(
  "Walking"  = "#E69F00",
  "Moderate" = "#CC79A7",
  "Vigorous" = "#000000"
)

model_cols <- c(
  "Unadjusted"        = "#4D4D4D",
  "Minimal adjustment" = "#56B4E9",
  "Fully adjusted"    = "#0072B2"
)


# =============================================================================
# 2. Shared plot theme
# =============================================================================

panel_theme <- ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 12),
    strip.text = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 12),
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )


# =============================================================================
# 3. Hazard-ratio table helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Total physical activity quintiles by sex
# -----------------------------------------------------------------------------

make_quintile_HR_table <- function(res_female, res_male) {
  
  female_tab <- res_female$tidied_models$quintile %>%
    dplyr::mutate(
      HR_95CI = sprintf(
        "%.2f (%.2f–%.2f)",
        estimate,
        conf.low,
        conf.high
      )
    ) %>%
    dplyr::select(
      model,
      quintile,
      Female = HR_95CI
    )
  
  male_tab <- res_male$tidied_models$quintile %>%
    dplyr::mutate(
      HR_95CI = sprintf(
        "%.2f (%.2f–%.2f)",
        estimate,
        conf.low,
        conf.high
      )
    ) %>%
    dplyr::select(
      model,
      quintile,
      Male = HR_95CI
    )
  
  dplyr::left_join(
    female_tab,
    male_tab,
    by = c("model", "quintile")
  ) %>%
    dplyr::arrange(
      factor(
        model,
        levels = c(
          "Unadjusted",
          "Minimal adjustment",
          "Fully adjusted"
        )
      ),
      quintile
    )
}


# -----------------------------------------------------------------------------
# Intensity-specific activity quintiles
# -----------------------------------------------------------------------------

make_intensity_HR_table <- function(tidy_df) {
  
  tidy_df %>%
    dplyr::mutate(
      HR_95CI = sprintf(
        "%.2f (%.2f–%.2f)",
        HR,
        CI_lower,
        CI_upper
      )
    ) %>%
    dplyr::select(
      sex,
      exposure,
      quintile,
      HR_95CI,
      p_trend
    ) %>%
    tidyr::pivot_wider(
      names_from = quintile,
      values_from = HR_95CI
    ) %>%
    dplyr::distinct() %>%
    dplyr::select(
      sex,
      exposure,
      Q2,
      Q3,
      Q4,
      Q5,
      p_trend
    )
}


# =============================================================================
# 4. Main forest-plot helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Activity-intensity quintiles
# -----------------------------------------------------------------------------

make_intensity_quintile_plot <- function(dat, plot_title) {
  
  plot_dat <- dat %>%
    dplyr::filter(
      exposure %in% c(
        "Walking PA",
        "Moderate PA",
        "Vigorous PA"
      )
    ) %>%
    dplyr::mutate(
      exposure = dplyr::recode(
        exposure,
        "Walking PA" = "Walking",
        "Moderate PA" = "Moderate",
        "Vigorous PA" = "Vigorous"
      ),
      exposure = factor(
        exposure,
        levels = c("Walking", "Moderate", "Vigorous")
      ),
      quintile = factor(
        quintile,
        levels = c("Q2", "Q3", "Q4", "Q5")
      )
    )
  
  y_lims <- range(
    c(plot_dat$CI_lower, plot_dat$CI_upper),
    na.rm = TRUE
  )
  
  y_pad <- diff(y_lims) * 0.15
  y_lims <- c(
    y_lims[1] - y_pad,
    y_lims[2] + y_pad
  )
  
  ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = quintile,
      y = HR,
      ymin = CI_lower,
      ymax = CI_upper,
      colour = exposure,
      group = exposure
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    ggplot2::geom_pointrange(
      position = ggplot2::position_dodge(width = 0.6),
      linewidth = 0.7
    ) +
    ggplot2::scale_colour_manual(
      values = activity_cols,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(ylim = y_lims) +
    ggplot2::labs(
      title = plot_title,
      x = "Activity quintile",
      y = "Hazard ratio (95% CI)",
      colour = "Activity type"
    ) +
    panel_theme
}


# -----------------------------------------------------------------------------
# Generic sensitivity comparison
# -----------------------------------------------------------------------------

make_sensitivity_plot <- function(
    tidy_df,
    plot_title,
    model_names,
    model_cols
) {
  
  plot_dat <- tidy_df %>%
    dplyr::filter(model %in% model_names) %>%
    dplyr::mutate(
      quintile = factor(
        quintile,
        levels = c("Q2", "Q3", "Q4", "Q5")
      ),
      model = factor(
        model,
        levels = model_names
      )
    )
  
  ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = quintile,
      y = estimate,
      ymin = conf.low,
      ymax = conf.high,
      colour = model,
      group = model
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    ggplot2::geom_pointrange(
      position = ggplot2::position_dodge(width = 0.6),
      linewidth = 0.7
    ) +
    ggplot2::scale_colour_manual(
      values = model_cols,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Total PA quintile",
      y = "Hazard ratio (95% CI)",
      colour = "Model"
    ) +
    panel_theme
}


# -----------------------------------------------------------------------------
# DXA femoral-neck BMD sensitivity comparison
# -----------------------------------------------------------------------------

make_dxa_forest_plot <- function(
    tidy_no_bmd,
    tidy_bmd,
    plot_title,
    adjusted_label = "Fully adjusted + femoral neck BMD"
) {
  
  plot_dat <- dplyr::bind_rows(
    tidy_no_bmd %>%
      dplyr::filter(
        model == "Fully adjusted",
        grepl("MET_total_quintile", term)
      ) %>%
      dplyr::mutate(
        model_plot = "Fully adjusted"
      ),
    
    tidy_bmd %>%
      dplyr::filter(
        model == "Fully adjusted + FN BMD",
        grepl("MET_total_quintile", term)
      ) %>%
      dplyr::mutate(
        model_plot = adjusted_label
      )
  ) %>%
    dplyr::mutate(
      quintile = factor(
        quintile,
        levels = c("Q2", "Q3", "Q4", "Q5")
      ),
      model_plot = factor(
        model_plot,
        levels = c(
          "Fully adjusted",
          adjusted_label
        )
      )
    )
  
  dxa_model_cols <- stats::setNames(
    c("#0072B2", "#D55E00"),
    c(
      "Fully adjusted",
      adjusted_label
    )
  )
  
  dodge <- ggplot2::position_dodge(width = 0.4)
  
  ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = quintile,
      y = estimate,
      colour = model_plot,
      group = model_plot
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = conf.low,
        ymax = conf.high
      ),
      position = dodge,
      width = 0.12
    ) +
    ggplot2::geom_point(
      position = dodge,
      size = 2.5
    ) +
    ggplot2::scale_colour_manual(
      values = dxa_model_cols,
      breaks = names(dxa_model_cols)
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Quintile of total physical activity",
      y = "Hazard ratio (95% CI)",
      colour = NULL
    ) +
    panel_theme
}


# =============================================================================
# 5. Cox and Fine-Gray prediction helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Cox restricted cubic spline curve
# -----------------------------------------------------------------------------

predict_cox_curve <- function(model, data, outcome, sex_label) {
  
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
  
  ref_data$ethnicity_derived <- factor(
    ref_data$ethnicity_derived,
    levels = levels(data$ethnicity_derived)
  )
  
  ref_data$education_level <- factor(
    ref_data$education_level,
    levels = levels(data$education_level),
    ordered = is.ordered(data$education_level)
  )
  
  new_data <- ref_data[rep(1, length(x_seq)), ]
  new_data$log_MET_total <- x_seq
  
  lp_new <- stats::predict(
    model,
    newdata = new_data,
    type = "lp"
  )
  
  lp_ref <- stats::predict(
    model,
    newdata = ref_data,
    type = "lp"
  )
  
  tibble::tibble(
    Outcome = outcome,
    Sex = sex_label,
    Model = "Cox",
    log_MET_total = x_seq,
    HR = exp(lp_new - as.numeric(lp_ref))
  )
}


# -----------------------------------------------------------------------------
# Fine-Gray restricted cubic spline curve
# -----------------------------------------------------------------------------

predict_fg_curve <- function(fg_obj, data, outcome, sex_label) {
  
  if (is.null(fg_obj$knots)) {
    stop(
      paste(
        "fg_obj$knots is missing.",
        "The Fine-Gray model object needs stored RCS knots."
      )
    )
  }
  
  if (is.null(fg_obj$cov_names)) {
    stop(
      paste(
        "fg_obj$cov_names is missing.",
        "The Fine-Gray model object needs stored covariate names."
      )
    )
  }
  
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
  
  ref_x <- stats::median(
    data$log_MET_total,
    na.rm = TRUE
  )
  
  knots <- fg_obj$knots
  
  make_cov <- function(x_values) {
    
    rcs_mat <- Hmisc::rcspline.eval(
      x_values,
      knots = knots,
      inclx = TRUE
    )
    
    colnames(rcs_mat) <- paste0(
      "rcs_log_MET_total_",
      seq_len(ncol(rcs_mat))
    )
    
    cov_data <- data.frame(
      rcs_mat,
      ethnicity_derived = names(
        sort(
          table(data$ethnicity_derived),
          decreasing = TRUE
        )
      )[1],
      height_clean = stats::median(
        data$height_clean,
        na.rm = TRUE
      ),
      weight_clean = stats::median(
        data$weight_clean,
        na.rm = TRUE
      ),
      tdi_raw = stats::median(
        data$tdi_raw,
        na.rm = TRUE
      ),
      education_level = names(
        sort(
          table(data$education_level),
          decreasing = TRUE
        )
      )[1]
    )
    
    cov_data$ethnicity_derived <- factor(
      cov_data$ethnicity_derived,
      levels = levels(data$ethnicity_derived)
    )
    
    cov_data$education_level <- factor(
      cov_data$education_level,
      levels = levels(data$education_level),
      ordered = is.ordered(data$education_level)
    )
    
    design_matrix <- stats::model.matrix(
      ~ .,
      data = cov_data
    )[, -1, drop = FALSE]
    
    design_matrix[, fg_obj$cov_names, drop = FALSE]
  }
  
  X_new <- make_cov(x_seq)
  X_ref <- make_cov(ref_x)
  
  lp_new <- as.vector(
    X_new %*% fg_obj$model$coef
  )
  lp_ref <- as.vector(
    X_ref %*% fg_obj$model$coef
  )
  
  tibble::tibble(
    Outcome = outcome,
    Sex = sex_label,
    Model = "Fine-Gray",
    log_MET_total = x_seq,
    HR = exp(lp_new - lp_ref)
  )
}


# =============================================================================
# 6. Supplementary table extraction helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Supplementary Table S4
# Continuous Cox models, spline LRT and sex interaction LRT
# -----------------------------------------------------------------------------

extract_main_model_results <- function(
    res_obj,
    outcome_label = "Hip fracture",
    exposure = "log_MET_total"
) {
  
  n_events <- res_obj$models$fully_adjusted_linear$nevent
  n_model <- res_obj$models$fully_adjusted_linear$n
  
  extract_hr <- function(model, model_name) {
    
    broom::tidy(
      model,
      exponentiate = TRUE,
      conf.int = TRUE
    ) %>%
      dplyr::filter(term == exposure) %>%
      dplyr::transmute(
        Outcome = outcome_label,
        N = n_model,
        Events = n_events,
        Analysis = "Continuous model",
        Result = model_name,
        Estimate = sprintf(
          "HR %.3f (95%% CI %.3f, %.3f)",
          estimate,
          conf.low,
          conf.high
        ),
        P_value = dplyr::if_else(
          p.value < 0.001,
          "<0.001",
          sprintf("%.3f", p.value)
        )
      )
  }
  
  tab_hr <- dplyr::bind_rows(
    extract_hr(
      res_obj$models$unadjusted,
      "Unadjusted"
    ),
    extract_hr(
      res_obj$models$minimally_adjusted,
      "Minimally adjusted"
    ),
    extract_hr(
      res_obj$models$fully_adjusted_linear,
      "Fully adjusted"
    )
  )
  
  lrt_spline <- res_obj$comparisons$linear_vs_spline_lrt
  
  tab_spline_lrt <- tibble::tibble(
    Outcome = outcome_label,
    N = n_model,
    Events = n_events,
    Analysis = "Non-linearity",
    Result = "Linear vs restricted cubic spline LRT",
    Estimate = sprintf(
      "Chi-square %.2f, df %s",
      lrt_spline$Chisq[2],
      lrt_spline$Df[2]
    ),
    P_value = dplyr::if_else(
      lrt_spline$`Pr(>|Chi|)`[2] < 0.001,
      "<0.001",
      sprintf(
        "%.3f",
        lrt_spline$`Pr(>|Chi|)`[2]
      )
    )
  )
  
  lrt_sex <- res_obj$comparisons$sex_interaction_lrt
  
  tab_sex_lrt <- tibble::tibble(
    Outcome = outcome_label,
    N = n_model,
    Events = n_events,
    Analysis = "Interaction",
    Result = "Sex × LTPA volume LRT",
    Estimate = sprintf(
      "Chi-square %.2f, df %s",
      lrt_sex$Chisq[2],
      lrt_sex$Df[2]
    ),
    P_value = dplyr::if_else(
      lrt_sex$`Pr(>|Chi|)`[2] < 0.001,
      "<0.001",
      sprintf(
        "%.3f",
        lrt_sex$`Pr(>|Chi|)`[2]
      )
    )
  )
  
  dplyr::bind_rows(
    tab_hr,
    tab_spline_lrt,
    tab_sex_lrt
  )
}


# -----------------------------------------------------------------------------
# Supplementary Table S5
# Proportional-hazards diagnostics for spline models
# -----------------------------------------------------------------------------

extract_ph_spline <- function(
    ph_object,
    outcome_label = "Hip fracture"
) {
  
  as.data.frame(ph_object$table) %>%
    tibble::rownames_to_column("Variable") %>%
    dplyr::transmute(
      Outcome = outcome_label,
      Model = "Fully adjusted spline",
      Variable = dplyr::recode(
        Variable,
        "rms::rcs(log_MET_total, 4)" = "Physical activity (RCS)",
        "sex_raw" = "Sex",
        "ethnicity_derived" = "Ethnicity",
        "height_clean" = "Height",
        "weight_clean" = "Weight",
        "tdi_raw" = "Townsend deprivation index",
        "education_level" = "Education",
        "GLOBAL" = "Global",
        .default = Variable
      ),
      Chi_square = round(chisq, 3),
      df = df,
      P_value = dplyr::case_when(
        is.na(p) ~ NA_character_,
        p < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p)
      )
    )
}


make_vig_walking_forest <- function(data, sex_label) {
  
  plot_dat <- data %>%
    dplyr::filter(sex == sex_label)
  
  ggplot(
    plot_dat,
    aes(
      x = quintile,
      y = HR,
      ymin = CI_lower,
      ymax = CI_upper,
      shape = model,
      group = model
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_errorbar(
      position = position_dodge(width = 0.45),
      width = 0.12,
      linewidth = 0.6
    ) +
    geom_point(
      position = position_dodge(width = 0.45),
      size = 2.7
    ) +
    scale_y_log10() +
    labs(
      title = sex_label,
      x = "Vigorous physical activity quintile",
      y = "Hazard ratio (95% CI)",
      shape = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "bottom",
      legend.justification = "center"
    )
}


# =============================================================================
# End of file
# =============================================================================