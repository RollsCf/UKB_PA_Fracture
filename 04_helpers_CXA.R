# =========================================================
#  Cross-Sectional Helper - Functions
# =========================================================


# =========================================================
# Continuous logistic regression
# =========================================================

fit_continuous_or <- function(data, outcome, exposure_log, confounders = NULL) {
  
  data[[outcome]] <- factor(data[[outcome]], levels = c("No", "Yes"))
  
  rhs <- c(exposure_log, confounders)
  formula <- reformulate(rhs, response = outcome)
  
  mod <- glm(formula, data = data, family = binomial())
  coef_tab <- summary(mod)$coefficients
  
  if (!(exposure_log %in% rownames(coef_tab))) {
    return(
      tibble::tibble(
        term = exposure_log,
        OR = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_
      )
    )
  }
  
  beta <- coef_tab[exposure_log, "Estimate"]
  se   <- coef_tab[exposure_log, "Std. Error"]
  pval <- coef_tab[exposure_log, "Pr(>|z|)"]
  
  tibble::tibble(
    term = exposure_log,
    OR = exp(beta),
    conf.low = exp(beta - 1.96 * se),
    conf.high = exp(beta + 1.96 * se),
    p.value = pval
  )
}


run_continuous_models <- function(data, outcomes, exposures, model_specs) {
  
  purrr::map_dfr(outcomes, function(outcome) {
    purrr::imap_dfr(model_specs, function(confounders, model_name) {
      purrr::map_dfr(exposures, function(exposure) {
        
        exposure_log <- paste0(exposure, "_log")
        
        fit_continuous_or(
          data = data,
          outcome = outcome,
          exposure_log = exposure_log,
          confounders = confounders
        ) %>%
          dplyr::mutate(
            Outcome = outcome,
            Exposure = exposure,
            Model = model_name
          )
      })
    })
  }) %>%
    dplyr::mutate(
      OR = round(OR, 3),
      conf.low = round(conf.low, 3),
      conf.high = round(conf.high, 3),
      p.value = round(p.value, 3)
    )
}


test_linearity <- function(data, outcomes, exposures, confounders, spline_knots = 4) {
  
  purrr::map_dfr(outcomes, function(outcome) {
    
    data[[outcome]] <- factor(data[[outcome]], levels = c("No", "Yes"))
    
    purrr::map_dfr(exposures, function(exposure) {
      
      exposure_log <- paste0(exposure, "_log")
      
      model_linear <- glm(
        reformulate(c(exposure_log, confounders), response = outcome),
        data = data,
        family = binomial()
      )
      
      model_rcs <- glm(
        as.formula(
          paste0(
            outcome, " ~ rms::rcs(", exposure_log, ", ", spline_knots, ") + ",
            paste(confounders, collapse = " + ")
          )
        ),
        data = data,
        family = binomial()
      )
      
      lrt <- anova(model_linear, model_rcs, test = "LRT")
      p_value <- lrt$`Pr(>Chi)`[2]
      
      tibble::tibble(
        Exposure = exposure,
        Outcome = outcome,
        LRT_p_value = p_value,
        Linear_Adequate = ifelse(p_value >= 0.05, "Yes", "No")
      )
    })
  })
}

# =========================================================
# Spline plotting for continuous exposure
# =========================================================


plot_spline_or <- function(data, outcome, exposure_log, confounders, spline_knots = 4) {
  
  data <- data %>%
    dplyr::select(dplyr::all_of(c(outcome, exposure_log, confounders))) %>%
    tidyr::drop_na()
  
  data[[outcome]] <- factor(data[[outcome]], levels = c("No", "Yes"))
  
  rhs <- paste0(
    "rms::rcs(", exposure_log, ", ", spline_knots, ")",
    if (length(confounders) > 0) paste0(" + ", paste(confounders, collapse = " + ")) else ""
  )
  
  formula <- as.formula(paste(outcome, "~", rhs))
  
  model <- glm(formula, data = data, family = binomial())
  
  x_seq <- seq(
    min(data[[exposure_log]], na.rm = TRUE),
    max(data[[exposure_log]], na.rm = TRUE),
    length.out = 100
  )
  
  newdata <- data.frame(x = x_seq)
  names(newdata)[1] <- exposure_log
  
  for (v in confounders) {
    
    if (is.numeric(data[[v]])) {
      
      newdata[[v]] <- mean(data[[v]], na.rm = TRUE)
      
    } else {
      
      newdata[[v]] <- factor(
        rep(levels(factor(data[[v]]))[1], nrow(newdata)),
        levels = levels(factor(data[[v]]))
      )
    }
  }
  
  ref_val <- median(data[[exposure_log]], na.rm = TRUE)
  refdata <- newdata[1, , drop = FALSE]
  refdata[[exposure_log]] <- ref_val
  
  X_new <- predict(model, newdata = newdata, type = "terms")
  
  lp_new <- predict(model, newdata = newdata, type = "link")
  lp_ref <- predict(model, newdata = refdata, type = "link")
  
  mm_new <- model.matrix(delete.response(terms(model)), newdata)
  mm_ref <- model.matrix(delete.response(terms(model)), refdata)
  
  X_diff <- sweep(mm_new, 2, mm_ref[1, ], "-")
  
  se_diff <- sqrt(rowSums((X_diff %*% vcov(model)) * X_diff))
  
  log_or <- as.numeric(lp_new - lp_ref)
  
  newdata$OR <- exp(log_or)
  newdata$lower <- exp(log_or - 1.96 * se_diff)
  newdata$upper <- exp(log_or + 1.96 * se_diff)
  
  ggplot2::ggplot(newdata, ggplot2::aes(x = .data[[exposure_log]], y = OR)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2
    ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = paste0("log(", gsub("_log$", "", exposure_log), " + 1)"),
      y = "Odds ratio (vs median)",
      title = paste(outcome, "vs", gsub("_log$", "", exposure_log), "(restricted cubic spline)")
    ) +
    ggplot2::theme_minimal()
}


save_spline_plots <- function(data, outcome, exposures, confounders, figures_dir, spline_knots = 4) {
  
  purrr::walk(exposures, function(exposure) {
    
    exposure_log <- paste0(exposure, "_log")
    
    p <- plot_spline_or(
      data = data,
      outcome = outcome,
      exposure_log = exposure_log,
      confounders = confounders,
      spline_knots = spline_knots
    )
    
    ggplot2::ggsave(
      filename = file.path(figures_dir, paste0("Fig_rcs_", outcome, "_", exposure_log, ".png")),
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
  })
}

# =========================================================
# Interaction testing: sex and age
# =========================================================

test_spline_interactions <- function(data,
                                     outcome = "SRF",
                                     exposure = "MET_total",
                                     confounders = full_adjustment,
                                     spline_knots = 4,
                                     age_group_var = "agegp_A0") {
  
  exposure_log <- paste0(exposure, "_log")
  
  dat <- data %>%
    dplyr::select(
      dplyr::all_of(c(outcome, exposure_log, confounders, age_group_var))
    ) %>%
    tidyr::drop_na()
  
  dat[[outcome]] <- factor(dat[[outcome]], levels = c("No", "Yes"))
  
  # -----------------------------
  # Sex interaction
  # -----------------------------
  
  conf_no_sex <- setdiff(confounders, "sex")
  
  f_sex_main <- stats::as.formula(
    paste0(
      outcome, " ~ rms::rcs(", exposure_log, ", ", spline_knots, ") + sex + ",
      paste(conf_no_sex, collapse = " + ")
    )
  )
  
  f_sex_int <- stats::as.formula(
    paste0(
      outcome, " ~ rms::rcs(", exposure_log, ", ", spline_knots, ") * sex + ",
      paste(conf_no_sex, collapse = " + ")
    )
  )
  
  m_sex_main <- glm(f_sex_main, data = dat, family = binomial())
  m_sex_int  <- glm(f_sex_int,  data = dat, family = binomial())
  
  sex_lrt <- anova(m_sex_main, m_sex_int, test = "LRT")
  
  # -----------------------------
  # Age-group interaction
  # -----------------------------
  
  conf_no_agegp <- setdiff(confounders, age_group_var)
  
  f_age_main <- stats::as.formula(
    paste0(
      outcome, " ~ rms::rcs(", exposure_log, ", ", spline_knots, ") + ",
      age_group_var, " + ",
      paste(conf_no_agegp, collapse = " + ")
    )
  )
  
  f_age_int <- stats::as.formula(
    paste0(
      outcome, " ~ rms::rcs(", exposure_log, ", ", spline_knots, ") * ",
      age_group_var, " + ",
      paste(conf_no_agegp, collapse = " + ")
    )
  )
  
  m_age_main <- glm(f_age_main, data = dat, family = binomial())
  m_age_int  <- glm(f_age_int,  data = dat, family = binomial())
  
  age_lrt <- anova(m_age_main, m_age_int, test = "LRT")
  
  tibble::tibble(
    outcome = outcome,
    exposure = exposure,
    interaction = c("Sex", "Age group"),
    LRT_p_value = c(
      sex_lrt$`Pr(>Chi)`[2],
      age_lrt$`Pr(>Chi)`[2]
    ),
    Evidence_of_interaction = dplyr::if_else(
      LRT_p_value < 0.05,
      "Yes",
      "No"
    )
  )
}

make_sex_specific_spline_plots <- function(data,
                                           outcome = "SRF",
                                           exposure_log = "MET_total_log",
                                           confounders,
                                           spline_knots = 4) {
  
  p_male <- plot_spline_or(
    data = data %>% dplyr::filter(sex == "Male"),
    outcome = outcome,
    exposure_log = exposure_log,
    confounders = confounders,
    spline_knots = spline_knots
  ) +
    ggplot2::labs(
      x = paste0("log(", gsub("_log$", "", exposure_log), " + 1)"),
      title = paste0("Male: ", gsub("_log$", "", exposure_log), " and ", outcome)
    )
  
  p_female <- plot_spline_or(
    data = data %>% dplyr::filter(sex == "Female"),
    outcome = outcome,
    exposure_log = exposure_log,
    confounders = confounders,
    spline_knots = spline_knots
  ) +
    ggplot2::labs(
      x = paste0("log(", gsub("_log$", "", exposure_log), " + 1)"),
      title = paste0("Female: ", gsub("_log$", "", exposure_log), " and ", outcome)
    )
  
  list(
    male = p_male,
    female = p_female
  )
}



# # =========================================================
# Quintile exposure models (Q1 as comparator)
# =========================================================

make_quintile_var <- function(data, exposure, n_bins = 5) {
  
  quint_var <- paste0(exposure, "_q")
  
  data[[quint_var]] <- NA_integer_
  ok <- !is.na(data[[exposure]])
  
  if (sum(ok) == 0) {
    return(data)
  }
  
  data[[quint_var]][ok] <- dplyr::ntile(data[[exposure]][ok], n_bins)
  
  data[[quint_var]] <- factor(
    data[[quint_var]],
    levels = 1:n_bins,
    labels = paste0("Q", 1:n_bins)
  )
  
  data[[quint_var]] <- stats::relevel(data[[quint_var]], ref = "Q1")
  
  data
}


fit_quintile_model <- function(data, exposure, outcome, confounders = NULL, n_bins = 5) {
  
  quint_var <- paste0(exposure, "_q")
  
  data <- make_quintile_var(data, exposure = exposure, n_bins = n_bins)
  data[[outcome]] <- factor(data[[outcome]], levels = c("No", "Yes"))
  
  vars_needed <- unique(c(outcome, exposure, quint_var, confounders))
  data_model <- data %>%
    dplyr::select(dplyr::all_of(vars_needed)) %>%
    tidyr::drop_na(dplyr::all_of(c(outcome, quint_var, confounders)))
  
  if (nrow(data_model) == 0) return(NULL)
  if (dplyr::n_distinct(stats::na.omit(data_model[[outcome]])) < 2) return(NULL)
  if (dplyr::n_distinct(stats::na.omit(data_model[[quint_var]])) < 2) return(NULL)
  
  form <- if (is.null(confounders) || length(confounders) == 0) {
    stats::as.formula(paste(outcome, "~", quint_var))
  } else {
    stats::as.formula(
      paste(outcome, "~", quint_var, "+", paste(confounders, collapse = " + "))
    )
  }
  
  model <- try(
    stats::glm(form, data = data_model, family = stats::binomial()),
    silent = TRUE
  )
  
  if (inherits(model, "try-error")) return(NULL)
  
  coef_tab <- summary(model)$coefficients
  quint_rows <- grep(paste0("^", quint_var), rownames(coef_tab))
  
  if (length(quint_rows) == 0) return(NULL)
  
  out_df <- tibble::tibble(
    term = rownames(coef_tab)[quint_rows],
    Odds_Ratio = exp(coef_tab[quint_rows, "Estimate"]),
    CI_lower   = exp(coef_tab[quint_rows, "Estimate"] - 1.96 * coef_tab[quint_rows, "Std. Error"]),
    CI_upper   = exp(coef_tab[quint_rows, "Estimate"] + 1.96 * coef_tab[quint_rows, "Std. Error"]),
    p_value    = coef_tab[quint_rows, "Pr(>|z|)"]
  )
  
  # Global p-value for quintile term
  anova_res <- anova(model, test = "Chisq")
  global_p <- anova_res$`Pr(>Chi)`[2]
  
  out_df <- out_df %>%
    dplyr::mutate(
      outcome = outcome,
      exposure = exposure,
      Global_p_value = c(global_p, rep(NA_real_, dplyr::n() - 1))
    )
  
  out_df
}


run_quintile_models_stratified <- function(data, strata_df, outcomes, exposures, model_specs, n_bins = 5) {
  
  results <- list()
  idx <- 1
  
  for (i in seq_len(nrow(strata_df))) {
    
    dat_stratum <- data
    
    for (v in names(strata_df)) {
      dat_stratum <- dat_stratum %>%
        dplyr::filter(as.character(.data[[v]]) == as.character(strata_df[[v]][i]))
    }
    
    dat_stratum <- droplevels(dat_stratum)
    
    if (nrow(dat_stratum) == 0) next
    
    for (outcome in outcomes) {
      for (exposure in exposures) {
        for (model_name in names(model_specs)) {
          
          res <- fit_quintile_model(
            data = dat_stratum,
            exposure = exposure,
            outcome = outcome,
            confounders = model_specs[[model_name]],
            n_bins = n_bins
          )
          
          if (is.null(res)) next
          
          for (v in names(strata_df)) {
            res[[v]] <- strata_df[[v]][i]
          }
          
          res$model <- model_name
          results[[idx]] <- res
          idx <- idx + 1
        }
      }
    }
  }
  
  dplyr::bind_rows(results)
}


format_quintile_results <- function(df, model_levels, digits = 5) {
  
  df %>%
    dplyr::mutate(
      model = factor(model, levels = model_levels),
      Odds_Ratio = round(Odds_Ratio, digits),
      CI_lower = round(CI_lower, digits),
      CI_upper = round(CI_upper, digits),
      p_value = round(p_value, digits),
      Global_p_value = round(Global_p_value, digits)
    )
}


prep_quintile_plot_df <- function(df, facet_vars, exposure_filter = NULL, outcome_filter = NULL) {
  
  plot_df <- df %>%
    dplyr::mutate(
      Quintile = sub(".*(Q[1-5])$", "\\1", term)
    ) %>%
    dplyr::filter(Quintile %in% c("Q2", "Q3", "Q4", "Q5"))
  
  if (!is.null(exposure_filter)) {
    plot_df <- plot_df %>% dplyr::filter(exposure %in% exposure_filter)
  }
  
  if (!is.null(outcome_filter)) {
    plot_df <- plot_df %>% dplyr::filter(outcome %in% outcome_filter)
  }
  
  ref_rows <- plot_df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(facet_vars, "outcome", "exposure", "model")))) %>%
    dplyr::mutate(
      term = NA_character_,
      Quintile = "Q1",
      Odds_Ratio = 1,
      CI_lower = 1,
      CI_upper = 1,
      p_value = NA_real_,
      Global_p_value = NA_real_
    )
  
  dplyr::bind_rows(plot_df, ref_rows) %>%
    dplyr::mutate(
      Quintile = factor(Quintile, levels = paste0("Q", 1:5))
    ) %>%
    dplyr::filter(
      is.finite(Odds_Ratio), Odds_Ratio > 0,
      is.finite(CI_lower),   CI_lower > 0,
      is.finite(CI_upper),   CI_upper > 0
    )
}



# =========================================================
# Categorical exposure models
# =========================================================

# =========================================================
# MODEL FUNCTION
# =========================================================

fit_categorical_model <- function(data, outcome, exposure, confounders = NULL,
                                  cat_type = c("quintile", "who", "precat"), n_bins = 5) {
  
  cat_type <- match.arg(cat_type)
  
  if (cat_type == "quintile") {
    data <- make_quintile_var(data, exposure = exposure, n_bins = n_bins)
    cat_var <- paste0(exposure, "_q")
  } else if (cat_type == "who") {
    data <- make_who_var(data, exposure = exposure)
    cat_var <- paste0(exposure, "_WHO")
  } else if (cat_type == "precat") {
    cat_var <- exposure
    data[[cat_var]] <- factor(data[[cat_var]])
  }
  
  vars_needed <- unique(c(outcome, cat_var, confounders))
  
  dat <- data %>%
    dplyr::select(dplyr::all_of(vars_needed)) %>%
    tidyr::drop_na()
  
  if (nrow(dat) == 0) return(NULL)
  
  dat[[outcome]] <- factor(dat[[outcome]], levels = c("No", "Yes"))
  
  if (dplyr::n_distinct(stats::na.omit(dat[[outcome]])) < 2) return(NULL)
  if (dplyr::n_distinct(stats::na.omit(dat[[cat_var]])) < 2) return(NULL)
  
  form <- if (is.null(confounders) || length(confounders) == 0) {
    reformulate(cat_var, response = outcome)
  } else {
    reformulate(c(cat_var, confounders), response = outcome)
  }
  
  model <- try(
    glm(form, data = dat, family = binomial()),
    silent = TRUE
  )
  
  if (inherits(model, "try-error")) return(NULL)
  
  res <- extract_rows_by_prefix(model, cat_var)
  if (is.null(res)) return(NULL)
  
  anova_res <- anova(model, test = "Chisq")
  global_p <- anova_res$`Pr(>Chi)`[2]
  
  res %>%
    dplyr::mutate(Global_p_value = c(global_p, rep(NA_real_, dplyr::n() - 1)))
}


# =========================================================
# STRATIFIED MODEL FUNCTION
# =========================================================

run_categorical_models_stratified <- function(data, strata_df, outcomes, exposures,
                                              model_specs, cat_type = c("quintile", "who", "precat"),
                                              n_bins = 5) {
  
  cat_type <- match.arg(cat_type)
  
  results <- list()
  idx <- 1
  
  for (i in seq_len(nrow(strata_df))) {
    
    dat_stratum <- data
    
    for (v in names(strata_df)) {
      dat_stratum <- dat_stratum %>%
        dplyr::filter(as.character(.data[[v]]) == as.character(strata_df[[v]][i]))
    }
    
    dat_stratum <- droplevels(dat_stratum)
    
    if (nrow(dat_stratum) == 0) next
    
    for (outcome in outcomes) {
      for (exposure in exposures) {
        for (model_name in names(model_specs)) {
          
          res <- fit_categorical_model(
            data = dat_stratum,
            outcome = outcome,
            exposure = exposure,
            confounders = model_specs[[model_name]],
            cat_type = cat_type,
            n_bins = n_bins
          )
          
          if (is.null(res)) next
          
          for (v in names(strata_df)) {
            res[[v]] <- strata_df[[v]][i]
          }
          
          res$outcome <- outcome
          res$exposure <- exposure
          res$model <- model_name
          
          results[[idx]] <- res
          idx <- idx + 1
        }
      }
    }
  }
  
  dplyr::bind_rows(results)
}

get_spline_prediction_data <- function(data, outcome, exposure, confounders, strata = NULL, spline_knots = 4) {
  
  dat <- data
  
  if (!is.null(strata)) {
    for (v in names(strata)) {
      dat <- dat %>% dplyr::filter(as.character(.data[[v]]) == as.character(strata[[v]]))
    }
  }
  
  dat <- dat %>%
    dplyr::select(dplyr::all_of(c(outcome, exposure, confounders))) %>%
    tidyr::drop_na()
  
  if (nrow(dat) == 0) return(NULL)
  
  dat[[outcome]] <- factor(dat[[outcome]], levels = c("No", "Yes"))
  
  rhs <- paste0(
    "rms::rcs(", exposure, ", ", spline_knots, ")",
    if (length(confounders) > 0) paste0(" + ", paste(confounders, collapse = " + ")) else ""
  )
  
  form <- as.formula(paste(outcome, "~", rhs))
  
  model <- try(glm(form, data = dat, family = binomial()), silent = TRUE)
  if (inherits(model, "try-error")) return(NULL)
  
  x_seq <- seq(
    min(dat[[exposure]], na.rm = TRUE),
    max(dat[[exposure]], na.rm = TRUE),
    length.out = 100
  )
  
  newdata <- data.frame(x = x_seq)
  names(newdata)[1] <- exposure
  
  for (v in confounders) {
    
    if (is.numeric(data[[v]])) {
      
      newdata[[v]] <- mean(data[[v]], na.rm = TRUE)
      
    } else {
      
      newdata[[v]] <- factor(
        rep(levels(factor(data[[v]]))[1], nrow(newdata)),
        levels = levels(factor(data[[v]]))
      )
    }
  }
  
  ref_val <- median(dat[[exposure]], na.rm = TRUE)
  refdata <- newdata[1, , drop = FALSE]
  refdata[[exposure]] <- ref_val
  
  lp_new <- predict(model, newdata = newdata, type = "link")
  lp_ref <- predict(model, newdata = refdata, type = "link")
  
  mm_new <- model.matrix(delete.response(terms(model)), newdata)
  mm_ref <- model.matrix(delete.response(terms(model)), refdata)
  
  X_diff <- sweep(mm_new, 2, mm_ref[1, ], "-")
  se_diff <- sqrt(rowSums((X_diff %*% vcov(model)) * X_diff))
  
  log_or <- as.numeric(lp_new - lp_ref)
  
  newdata$OR <- exp(log_or)
  newdata$lower <- exp(log_or - 1.96 * se_diff)
  newdata$upper <- exp(log_or + 1.96 * se_diff)
  
  newdata$outcome <- outcome
  
  if (!is.null(strata)) {
    for (v in names(strata)) newdata[[v]] <- strata[[v]]
  }
  
  newdata
}


extract_or_term <- function(model, term_name) {
  coef_tab <- summary(model)$coefficients
  
  if (!term_name %in% rownames(coef_tab)) {
    return(
      tibble::tibble(
        term = term_name,
        OR = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  beta <- coef_tab[term_name, "Estimate"]
  se   <- coef_tab[term_name, "Std. Error"]
  pval <- coef_tab[term_name, "Pr(>|z|)"]
  
  tibble::tibble(
    term = term_name,
    OR = exp(beta),
    CI_lower = exp(beta - 1.96 * se),
    CI_upper = exp(beta + 1.96 * se),
    p_value = pval
  )
}


extract_rows_by_prefix <- function(model, prefix) {
  coef_tab <- summary(model)$coefficients
  rows <- grep(paste0("^", prefix), rownames(coef_tab))
  
  if (length(rows) == 0) return(NULL)
  
  tibble::tibble(
    term = rownames(coef_tab)[rows],
    OR = exp(coef_tab[rows, "Estimate"]),
    CI_lower = exp(coef_tab[rows, "Estimate"] - 1.96 * coef_tab[rows, "Std. Error"]),
    CI_upper = exp(coef_tab[rows, "Estimate"] + 1.96 * coef_tab[rows, "Std. Error"]),
    p_value = coef_tab[rows, "Pr(>|z|)"]
  )
}


fit_single_term_model <- function(data, outcome, exposure, confounders = NULL) {
  
  vars_needed <- unique(c(outcome, exposure, confounders))
  
  dat <- data %>%
    dplyr::select(dplyr::all_of(vars_needed)) %>%
    tidyr::drop_na()
  
  if (nrow(dat) == 0) return(NULL)
  
  dat[[outcome]] <- factor(dat[[outcome]], levels = c("No", "Yes"))
  
  if (dplyr::n_distinct(stats::na.omit(dat[[outcome]])) < 2) return(NULL)
  if (dplyr::n_distinct(stats::na.omit(dat[[exposure]])) < 2) return(NULL)
  
  form <- if (is.null(confounders) || length(confounders) == 0) {
    reformulate(exposure, response = outcome)
  } else {
    reformulate(c(exposure, confounders), response = outcome)
  }
  
  model <- try(
    glm(form, data = dat, family = binomial()),
    silent = TRUE
  )
  
  if (inherits(model, "try-error")) return(NULL)
  
  extract_or_term(model, exposure)
}


run_single_term_models_stratified <- function(data, strata_df, outcomes, exposures, confounders) {
  
  results <- list()
  idx <- 1
  
  for (i in seq_len(nrow(strata_df))) {
    
    dat_stratum <- data
    
    for (v in names(strata_df)) {
      dat_stratum <- dat_stratum %>%
        dplyr::filter(as.character(.data[[v]]) == as.character(strata_df[[v]][i]))
    }
    
    dat_stratum <- droplevels(dat_stratum)
    
    if (nrow(dat_stratum) == 0) next
    
    for (outcome in outcomes) {
      for (exposure in exposures) {
        
        res <- fit_single_term_model(
          data = dat_stratum,
          outcome = outcome,
          exposure = exposure,
          confounders = confounders
        )
        
        if (is.null(res)) next
        
        for (v in names(strata_df)) {
          res[[v]] <- strata_df[[v]][i]
        }
        
        res$outcome <- outcome
        res$exposure <- exposure
        
        results[[idx]] <- res
        idx <- idx + 1
      }
    }
  }
  
  dplyr::bind_rows(results)
}


# =========================================================
# Proportion models within MET_total quintiles
# =========================================================

run_prop_within_quintiles_by_sex <- function(data, sex_levels, confounders) {
  
  results <- list()
  idx <- 1
  
  for (sex_group in sex_levels) {
    
    dat <- data %>%
      dplyr::filter(sex == sex_group) %>%
      dplyr::select(
        SRF,
        MET_total,
        prop_vig_10,
        prop_mod_10,
        prop_walk_10,
        dplyr::all_of(confounders)
      ) %>%
      tidyr::drop_na()
    
    if (nrow(dat) == 0) next
    
    dat$SRF <- factor(dat$SRF, levels = c("No", "Yes"))
    if (dplyr::n_distinct(dat$SRF) < 2) next
    
    dat <- make_quintile_var(dat, exposure = "MET_total", n_bins = 5)
    
    for (q in levels(dat$MET_total_q)) {
      
      dat_q <- dat %>%
        dplyr::filter(MET_total_q == q)
      
      if (nrow(dat_q) == 0) next
      if (dplyr::n_distinct(dat_q$SRF) < 2) next
      
      for (term_name in c("prop_vig_10", "prop_mod_10", "prop_walk_10")) {
        
        mod <- try(
          stats::glm(
            stats::reformulate(c(term_name, confounders), response = "SRF"),
            data = dat_q,
            family = stats::binomial()
          ),
          silent = TRUE
        )
        
        if (inherits(mod, "try-error")) next
        
        res <- extract_or_term(mod, term_name) %>%
          dplyr::mutate(
            sex = sex_group,
            model = paste0(term_name, " within ", q)
          )
        
        results[[idx]] <- res
        idx <- idx + 1
      }
    }
  }
  
  dplyr::bind_rows(results)
}