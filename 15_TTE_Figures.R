# To create the figures scripts 06-09 need to have been run first

source(here::here("scripts/00_setup.R"))
source(here("scripts/03_helpers_general.R"))
source(here("scripts/04_helpers_CXA.R"))
source(here("scripts/05_helpers_survival.R"))
source(here("scripts/07_helpers_TTE_figures.R"))

# =========================================================
# Figure 1: Exclusion flow table
# =========================================================

tab_exc <- readRDS(file.path(DATA_DERIVED, "TTE_tab_exclusions.rds"))

# Build the PNG file in FIGURES_DIR
png(
  filename = file.path(TTE_FIGURES_DIR, "Fig_1_TTE_A0_Exclusion_flow.png"),
  width = 1200,
  height = 800,
  res = 150
)

# Render the table as a figure
gridExtra::grid.table(tab_exc)

# Close the device
dev.off()


# =========================================================
# Table 1: Baseline characteristics by sex
# Uses complete-case analytic cohort
# =========================================================

complete_case_dat <- readRDS(
  file.path(DATA_DERIVED, "TTE_A0_complete_case.rds")
) %>%
  dplyr::mutate(
    follow_up_years = age_exit_hip_PA - age_entry_PA
  )

complete_case_dat <- complete_case_dat %>%
  dplyr::mutate(
    hip_fracture_table = factor(
      event_hip_PA,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

tbl <- complete_case_dat %>% 
  dplyr::select(
    sex,
    age_A0,
    ethnicity_derived,
    tdi_raw,
    education_level,
    weight_clean,
    height_clean,
    follow_up_years,
    hip_fracture_table,
    cc_MET_mod_trunc,
    cc_MET_vig_trunc,
    cc_MET_walk_trunc,
    cc_MET_total_trunc
  ) %>% 
  gtsummary::tbl_summary(
    by = sex,
    statistic = list(
      c(age_A0, tdi_raw, weight_clean, height_clean, follow_up_years) ~ "{mean} ({sd})",
      c(cc_MET_mod_trunc, cc_MET_vig_trunc, cc_MET_walk_trunc, cc_MET_total_trunc) ~ "{median} ({p25}, {p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      c(age_A0, tdi_raw, weight_clean, height_clean, follow_up_years) ~ 1,
      c(cc_MET_mod_trunc, cc_MET_vig_trunc, cc_MET_walk_trunc, cc_MET_total_trunc) ~ 0
    ),
    type = list(
      hip_fracture_table ~ "dichotomous"
    ),
    value = list(
      hip_fracture_table ~ "Yes"
    ),
    label = list(
      follow_up_years ~ "Follow-up (years), mean (SD)",
      age_A0 ~ "Age (years), mean (SD)",
      ethnicity_derived ~ "Ethnicity, n (%)",
      tdi_raw ~ "Townsend deprivation index, mean (SD)",
      education_level ~ "Education level, n (%)",
      weight_clean ~ "Weight (kg), mean (SD)",
      height_clean ~ "Height (cm), mean (SD)",
      hip_fracture_table ~ "Hip fracture, n (%)",
      cc_MET_mod_trunc ~ "Moderate activity (MET-min/week), median (IQR)",
      cc_MET_vig_trunc ~ "Vigorous activity (MET-min/week), median (IQR)",
      cc_MET_walk_trunc ~ "Walking (MET-min/week), median (IQR)",
      cc_MET_total_trunc ~ "Leisure time physical activity (MET-min/week), median (IQR)"
    ),
    missing_text = "Missing"
  )%>%
  gtsummary::add_overall(last = FALSE) %>%
  gtsummary::modify_header(label ~ "**Variable**") %>%
  gtsummary::modify_spanning_header(
    c("stat_0", "stat_1", "stat_2") ~ "**Sex**"
  ) %>%
  gtsummary::modify_caption(
    "**Table 1. Baseline characteristics of the analytic cohort overall and by sex**"
  ) %>%
  gtsummary::modify_footnote(
    gtsummary::all_stat_cols() ~ NA
  )

tbl_gt <- tbl %>%
  gtsummary::as_gt() %>%
  gt::tab_source_note(
    source_note = gt::md(
      "Continuous demographic and anthropometric variables and follow-up are presented as mean (SD). Physical activity variables are presented as median (IQR). Categorical variables are presented as n/N (%)."
    )
  )

tbl_gt

tbl_gt %>%
  gt::gtsave(
    filename = file.path(TTE_TABLES_DIR, "Table1.docx")
  )

# =========================================================
# Supplementary Table 1: ICD-10 codes not provided here
# =========================================================

# =========================================================
# Supplementary Table 2 and 3: Response rate and degree of missing data
# =========================================================
# =========================================================
# Missingness and activity-data supplementary tables
# Uses FULL processed dataset, not complete-case dataset
# =========================================================

# Load full processed data
full_analysis_dat <- readRDS(file.path(DATA_DERIVED, "TTE_analysis_dat.Rds"))

# Clean variable names and rename key variables
full_analysis_dat <- full_analysis_dat %>%
  dplyr::rename_with(~ gsub("_clean|_raw", "", .x)) %>%
  dplyr::rename(
    SRF = self_reported_fracture_A0,
    WC  = waist_circ
  ) %>%
  dplyr::mutate(
    age_A0 = as.numeric(age_A0),
    date_assess_A0 = data.table::as.IDate(as.Date(date_assess_A0)),
    lost_to_fu     = data.table::as.IDate(as.Date(lost_to_fu))
  )

# Count EID in full cohort prior to age restriction
full_analysis_dat %>%
  dplyr::summarise(n_unique = dplyr::n_distinct(eid))

# Create full middle-aged cohort, before complete-case exclusions
middle_aged_full_dat <- full_analysis_dat %>%
  dplyr::filter(dplyr::between(age_A0, 40, 65))

var_labels <- c(
  ethnicity_derived = "Ethnicity",
  num_days_mod = "Days/week moderate activity",
  num_day_walk = "Days/week walking",
  num_days_vig = "Days/week vigorous activity",
  cc_MET_mod_trunc = "Moderate activity volume (MET-min/week)",
  cc_MET_vig_trunc = "Vigorous activity volume (MET-min/week)",
  cc_MET_walk_trunc = "Leisure walking volume (MET-min/week)",
  cc_MET_total_trunc = "Leisure time physical activity volume (MET-min/week)",
  education_level = "Education",
  tdi = "Townsend deprivation index",
  weight = "Weight (kg)",
  height = "Height (cm)",
  BMI = "Body mass index (kg/m²)"
)

response_rate <- function(df, var) {
  df %>%
    dplyr::summarise(
      Answered = sum(!is.na(.data[[var]])),
      Missing = sum(is.na(.data[[var]])),
      Total = dplyr::n(),
      Answered_pct = round(Answered / Total * 100, 2),
      Missing_pct = round(Missing / Total * 100, 2)
    ) %>%
    dplyr::mutate(
      Variable = unname(var_labels[[var]])
    ) %>%
    dplyr::select(
      Variable,
      Answered,
      Answered_pct,
      Missing,
      Missing_pct
    )
}

Response_rates <- lapply(
  names(var_labels),
  function(v) response_rate(middle_aged_full_dat, v)
) %>%
  dplyr::bind_rows()

Response_rates

save_table_word(
  df = Response_rates,
  table_number = "S2_response_rate",
  folder_path = TTE_TABLES_DIR,
  title = "Response rate and degree of missing data"
)

##################################################################################################
# =========================================================
# Supplementary Table 3: With vs without activity data
# =========================================================

supp_activity_dat <- middle_aged_full_dat %>%
  dplyr::filter(!is.na(SRF)) %>%
  dplyr::mutate(
    activity_data = dplyr::if_else(
      !is.na(cc_MET_walk_trunc) &
        !is.na(cc_MET_mod_trunc) &
        !is.na(cc_MET_vig_trunc),
      "Activity data available",
      "Activity data missing"
    ),
    activity_data = factor(
      activity_data,
      levels = c("Activity data available", "Activity data missing")
    )
  )

get_smd <- function(var, data) {
  
  x <- data[[var]]
  g <- data$activity_data
  
  if (is.numeric(x)) {
    
    m1 <- mean(x[g == levels(g)[1]], na.rm = TRUE)
    m2 <- mean(x[g == levels(g)[2]], na.rm = TRUE)
    
    s1 <- stats::sd(x[g == levels(g)[1]], na.rm = TRUE)
    s2 <- stats::sd(x[g == levels(g)[2]], na.rm = TRUE)
    
    pooled_sd <- sqrt((s1^2 + s2^2) / 2)
    smd <- abs(m1 - m2) / pooled_sd
    
  } else {
    
    tab <- prop.table(table(x, g), margin = 2)
    
    if (ncol(tab) < 2) {
      smd <- NA_real_
    } else {
      smd <- max(abs(tab[, 1] - tab[, 2]), na.rm = TRUE)
    }
  }
  
  tibble::tibble(
    variable = var,
    SMD = round(smd, 3)
  )
}

vars_to_test <- c(
  "sex",
  "age_A0",
  "ethnicity_derived",
  "tdi",
  "education_level",
  "weight",
  "height"
)

smd_values <- dplyr::bind_rows(
  lapply(vars_to_test, get_smd, data = supp_activity_dat)
)

supp_tbl_activity <- supp_activity_dat %>%
  dplyr::select(
    activity_data,
    sex,
    age_A0,
    ethnicity_derived,
    tdi,
    education_level,
    weight,
    height
  ) %>%
  gtsummary::tbl_summary(
    by = activity_data,
    statistic = list(
      c(age_A0, tdi, weight, height) ~ "{mean} ({sd})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      c(age_A0, tdi, weight, height) ~ 1
    ),
    type = gtsummary::all_categorical() ~ "categorical",
    label = list(
      sex ~ "Sex, n (%)",
      age_A0 ~ "Age (years), mean (SD)",
      ethnicity_derived ~ "Ethnicity, n (%)",
      tdi ~ "Townsend deprivation index, mean (SD)",
      education_level ~ "Education level, n (%)",
      weight ~ "Weight (kg), mean (SD)",
      height ~ "Height (cm), mean (SD)"
    ),
    missing_text = "Missing"
  ) %>%
  gtsummary::add_overall(last = FALSE) %>%
  gtsummary::modify_table_body(
    ~ .x %>%
      dplyr::left_join(smd_values, by = "variable") %>%
      dplyr::mutate(
        SMD = dplyr::if_else(row_type == "label", as.character(SMD), "")
      )
  ) %>%
  gtsummary::modify_header(
    label ~ "**Variable**",
    SMD ~ "**SMD**"
  ) %>%
  gtsummary::modify_caption(
    "**Supplementary Table. Baseline characteristics of participants with and without physical activity data**"
  ) %>%
  gtsummary::modify_footnote(
    SMD ~ "SMD = standardised mean difference."
  )

supp_tbl_activity

supp_tbl_activity_df <- supp_tbl_activity$table_body %>%
  dplyr::select(label, stat_0, stat_1, stat_2, SMD) %>%
  dplyr::rename(
    Variable = label,
    Overall = stat_0,
    `Activity data available` = stat_1,
    `Activity data missing` = stat_2
  )

save_table_word(
  df = supp_tbl_activity_df,
  table_number = "S3_missing_data_smd",
  folder_path = TTE_TABLES_DIR,
  title = "Missing data comparison"
)

# =========================================================
# Supplementary Figure 2: Histogram of PA and log transformed PA
# =========================================================


met_raw_vars <- c(
  "cc_MET_mod_trunc",
  "cc_MET_vig_trunc",
  "cc_MET_walk_trunc",
  "cc_MET_total_trunc"
)

met_labels <- c(
  cc_MET_mod_trunc   = "Moderate",
  cc_MET_vig_trunc   = "Vigorous",
  cc_MET_walk_trunc  = "Walking",
  cc_MET_total_trunc = "Total LTPA"
)

# Create one long dataset from the raw variables
hist_df <- full_analysis_dat %>%
  dplyr::select(dplyr::all_of(met_raw_vars)) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "variable",
    values_to = "raw_value"
  ) %>%
  dplyr::mutate(
    variable = factor(
      variable,
      levels = met_raw_vars,
      labels = unname(met_labels[met_raw_vars])
    ),
    log_value = log1p(raw_value)
  )

hist_df_plot <- hist_df %>%
  dplyr::filter(
    is.finite(raw_value),
    is.finite(log_value)
  )

p_raw <- ggplot2::ggplot(
  hist_df_plot,
  ggplot2::aes(x = raw_value)
) +
  ggplot2::geom_histogram(
    bins = 50,
    colour = "white"
  ) +
  ggplot2::facet_wrap(
    ~ variable,
    scales = "free",
    ncol = 2
  ) +
  ggplot2::labs(
    title = "Before log(x + 1) transformation",
    x = "MET-min/week",
    y = "Count"
  ) +
  ggplot2::theme_bw()

p_log <- ggplot2::ggplot(
  hist_df_plot,
  ggplot2::aes(x = log_value)
) +
  ggplot2::geom_histogram(
    bins = 50,
    colour = "white"
  ) +
  ggplot2::facet_wrap(
    ~ variable,
    scales = "free",
    ncol = 2
  ) +
  ggplot2::labs(
    title = "After log(x + 1) transformation",
    x = "log(MET-min/week + 1)",
    y = "Count"
  ) +
  ggplot2::theme_bw()

p_supp_s2 <- patchwork::wrap_plots(
  p_raw,
  p_log,
  ncol = 1
)

p_supp_s2

ggplot2::ggsave(
  filename = file.path(
    TTE_FIGURES_DIR,
    "Fig_S2_MET_histograms_raw_log.png"
  ),
  plot = p_supp_s2,
  width = 12,
  height = 10,
  dpi = 300
)


# =========================================================
# Supplementary Table S4:
# Continuous Cox models, spline LRT and sex interaction LRT
# =========================================================

res_hip <- readRDS(
  file.path(
    DATA_DERIVED,
    "res_hip_cont.rds"
  )
)

tab_hip_continuous <- extract_main_model_results(
  res_obj = res_hip,
  outcome_label = "Hip fracture",
  exposure = "log_MET_total"
)

View(tab_hip_continuous)

save_table_word(
  df = tab_hip_continuous,
  table_number = "S4",
  folder_path = TTE_TABLES_DIR,
  title = paste(
    "Continuous Cox models and likelihood ratio tests",
    "for leisure-time physical activity volume and hip fracture"
  )
)

# =========================================================
# Supplementary Table S5
# Proportional hazards test for pooled hip fracture spline model
# =========================================================

ph_hip_pooled <- readRDS(
  file.path(DATA_DERIVED, "ph_hip_pooled.rds")
)


# -----------------------------
# Create table
# -----------------------------

tab_ph_hip <- extract_ph_spline(
  ph_object = ph_hip_pooled$full_spline,
  outcome_label = "Hip fracture"
)

View(tab_ph_hip)


# -----------------------------
# Save table
# -----------------------------

save_table_word(
  df = tab_ph_hip,
  table_number = "S5_ph_hip_pooled_spline_model",
  folder_path = TTE_TABLES_DIR,
  title = paste(
    "Proportional hazards tests for the pooled",
    "fully adjusted spline model of hip fracture"
  )
)


#  =========================================================
# Figure 2 Male and Female Spline Hip Fracture
# =========================================================


hip_spline_plot_data <- readRDS(
  file.path(DATA_DERIVED, "hip_spline_plot_data.rds")
)

make_spline_plot_from_data <- function(dat) {
  
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = log_MET_total,
      y = HR
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = HR_low,
        ymax = HR_high
      ),
      fill = "grey70",
      alpha = 0.6
    ) +
    ggplot2::geom_line(
      linewidth = 0.9
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    panel_theme +
    ggplot2::labs(
      x = "Leisure time physical activity (log[MET-min/week + 1])",
      y = "Hazard ratio"
    )
}

p_hip_female_spline <- make_spline_plot_from_data(
  hip_spline_plot_data$female
) +
  ggplot2::labs(title = "Female")

p_hip_male_spline <- make_spline_plot_from_data(
  hip_spline_plot_data$male
) +
  ggplot2::labs(title = "Male")

p_hip_female_spline
p_hip_male_spline


# Add panel titles
p_hip_female_spline <- p_hip_female_spline +
  ggplot2::labs(title = "Female")

p_hip_male_spline <- p_hip_male_spline +
  ggplot2::labs(title = "Male")


# -----------------------------
# Combine plots side by side
# -----------------------------

fig_2_hip_spline_sex_panel <-
  p_hip_female_spline +
  p_hip_male_spline +
  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

fig_2_hip_spline_sex_panel


# -----------------------------
# Save figure
# -----------------------------

ggplot2::ggsave(
  filename = file.path(
    TTE_FIGURES_DIR,
    "Fig_2_hip_spline_by_sex.png"
  ),
  plot = fig_2_hip_spline_sex_panel,
  width = 14,
  height = 6,
  dpi = 300
)

# =====================================================================
# Table 2: Physical activity and Musculoskeletal markers at baseline
# =====================================================================


TTE_A0_health_profile <- readRDS(
  file.path(DATA_DERIVED, "TTE_A0_health_profile.rds")
)


fmt_num <- function(x, digits = 1) {
  formatC(
    x,
    format = "f",
    digits = digits,
    big.mark = ","
  )
}

fmt_n_pct <- function(n, pct) {
  paste0(
    formatC(n, format = "d", big.mark = ","),
    " (",
    fmt_num(pct, 1),
    "%)"
  )
}


make_panel_table <- function(data, sex_label) {
  
  dat <- data %>%
    dplyr::filter(sex_raw == sex_label) %>%
    dplyr::mutate(
      total_MET_hrwk  = cc_MET_total_trunc / 60,
      total_hrwk      = cc_mins_wk_total_trunc / 60,
      walking_hrwk    = cc_mins_wk_walk_trunc / 60,
      moderate_hrwk   = cc_mins_wk_mod_trunc / 60,
      vigorous_hrwk   = cc_mins_wk_vig_trunc / 60,
      brisk_walking   = usual_walking_pace_raw == "Brisk pace",
      recurrent_falls = self_reported_falls_clean == "More than one fall"
    )
  
  
  # Continuous variables: median and IQR in one row
  cont_summary <- function(var, label, digits = 1) {
    
    dat %>%
      dplyr::group_by(MET_total_quintile_f) %>%
      dplyr::summarise(
        median = median(.data[[var]], na.rm = TRUE),
        p25 = quantile(.data[[var]], 0.25, na.rm = TRUE),
        p75 = quantile(.data[[var]], 0.75, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        Value = paste0(
          fmt_num(median, digits),
          " (",
          fmt_num(p25, digits),
          "–",
          fmt_num(p75, digits),
          ")"
        )
      ) %>%
      dplyr::select(
        MET_total_quintile_f,
        Value
      ) %>%
      tidyr::pivot_wider(
        names_from = MET_total_quintile_f,
        values_from = Value
      ) %>%
      dplyr::mutate(
        Characteristic = label,
        Statistic = "Median (IQR)",
        .before = 1
      )
  }
  
  
  # Categorical variables: n (%) in one row
  cat_summary <- function(var, label) {
    
    dat %>%
      dplyr::group_by(MET_total_quintile_f) %>%
      dplyr::summarise(
        n = sum(.data[[var]] == TRUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        pct = 100 * n / sum(n),
        Value = fmt_n_pct(n, pct)
      ) %>%
      dplyr::select(
        MET_total_quintile_f,
        Value
      ) %>%
      tidyr::pivot_wider(
        names_from = MET_total_quintile_f,
        values_from = Value
      ) %>%
      dplyr::mutate(
        Characteristic = label,
        Statistic = "n (% across quintiles)",
        .before = 1
      )
  }
  
  
  dplyr::bind_rows(
    cont_summary(
      "total_MET_hrwk",
      "LTPA, MET-hours/week"
    ),
    cont_summary(
      "total_hrwk",
      "LTPA, hours/week"
    ),
    cont_summary(
      "walking_hrwk",
      "Walking, hours/week"
    ),
    cont_summary(
      "moderate_hrwk",
      "Moderate activity, hours/week"
    ),
    cont_summary(
      "vigorous_hrwk",
      "Vigorous activity, hours/week"
    ),
    cat_summary(
      "brisk_walking",
      "Brisk walking pace"
    ),
    cat_summary(
      "recurrent_falls",
      "More than one fall"
    ),
    cont_summary(
      "grip_strength_mean",
      "Grip strength, kg"
    ),
    cont_summary(
      "heel_BMD",
      "Heel BMD",
      digits = 2
    )
  ) %>%
    dplyr::mutate(
      Panel = paste0("Panel: ", sex_label),
      .before = 1
    )
}


table_2_dat <- dplyr::bind_rows(
  make_panel_table(
    TTE_A0_health_profile,
    "Female"
  ),
  make_panel_table(
    TTE_A0_health_profile,
    "Male"
  )
)


table_2_gt <- table_2_dat %>%
  gt::gt(
    groupname_col = "Panel"
  ) %>%
  gt::tab_header(
    title = paste(
      "Table 2. Activity composition and selected",
      "musculoskeletal characteristics across",
      "physical activity quintiles by sex"
    )
  ) %>%
  gt::cols_label(
    Characteristic = "Characteristic",
    Statistic = "Statistic"
  ) %>%
  gt::tab_source_note(
    source_note = paste(
      "Values are median (interquartile range) for continuous variables",
      "and n (%) for categorical variables.",
      "MET-hours/week were calculated as MET-min/week divided by 60."
    )
  )


gt::gtsave(
  table_2_gt,
  filename = "Table_2_activity_composition_msk_by_quintile_sex.docx",
  path = TTE_TABLES_DIR
)



# =========================================================
# Supplementary table S6: PA quintiles and Hip fracture risk
# HRs by sex
# =========================================================


res_hip_quint_female <- readRDS(
  file.path(DATA_DERIVED, "res_hip_quint_female.rds")
)

res_hip_quint_male <- readRDS(
  file.path(DATA_DERIVED, "res_hip_quint_male.rds")
)

supp_table_hip_quint_HR <- make_quintile_HR_table(
  res_hip_quint_female,
  res_hip_quint_male
)

View(supp_table_hip_quint_HR)

save_table_word(
  df = supp_table_hip_quint_HR,
  table_number = "S6_hip_quintile_HR",
  folder_path = TTE_TABLES_DIR,
  title = "Hazard ratios for incident hip fracture across sex-specific quintiles of total physical activity"
)



# =========================================================
# Figure 3: Hip Quintile plots
# =========================================================

res_hip_quint_female <- readRDS(
  file.path(DATA_DERIVED, "res_hip_quint_female.rds")
)

res_hip_quint_male <- readRDS(
  file.path(DATA_DERIVED, "res_hip_quint_male.rds")
)

p_hip_quint_female <- plot_pa_quintiles_models(
  tidy_quint_df =
    res_hip_quint_female$tidied_models$quintile,
  outcome_name =
    "hip fracture - Female"
) +
  ggplot2::labs(title = "Female")

p_hip_quint_male <- plot_pa_quintiles_models(
  tidy_quint_df =
    res_hip_quint_male$tidied_models$quintile,
  outcome_name =
    "hip fracture - Male"
) +
  ggplot2::labs(title = "Male")


# Add panel titles if not already present
p_hip_quint_female <- p_hip_quint_female +
  labs(title = "Female")

p_hip_quint_male <- p_hip_quint_male +
  labs(title = "Male")


# -----------------------------
# Combine plots
# -----------------------------

fig_3_hip_quintile_forest_sex_panel <-
  p_hip_quint_female +
  p_hip_quint_male +
  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

fig_3_hip_quintile_forest_sex_panel


# -----------------------------
# Save figure
# -----------------------------

ggsave(
  filename = file.path(
    TTE_FIGURES_DIR,
    "Fig_3_hip_quintile_forest_by_sex.png"
  ),
  plot = fig_3_hip_quintile_forest_sex_panel,
  width = 14,
  height = 6,
  dpi = 300
)



# =========================================================
# Supplementary Figure S3:
# Fine-Gray versus Cox spline comparison
# =========================================================

res_hip_female <- readRDS(
  file.path(
    DATA_DERIVED,
    "res_hip_female.rds"
  )
)

res_hip_male <- readRDS(
  file.path(
    DATA_DERIVED,
    "res_hip_male.rds"
  )
)

fg_hip_results <- readRDS(
  file.path(
    DATA_DERIVED,
    "fg_hip_results.rds"
  )
)

fg_hip_female <- fg_hip_results$female
fg_hip_male   <- fg_hip_results$male

female_dat_fg <- fg_hip_results$female_data
male_dat_fg   <- fg_hip_results$male_data


# =========================================================
# Build prediction dataset
# =========================================================

curve_dat <- dplyr::bind_rows(
  
  predict_cox_curve(
    res_hip_female$models$fully_adjusted_spline,
    female_dat,
    "Hip fracture",
    "Female"
  ),
  
  predict_fg_curve(
    fg_hip_female,
    female_dat_fg,
    "Hip fracture",
    "Female"
  ),
  
  predict_cox_curve(
    res_hip_male$models$fully_adjusted_spline,
    male_dat,
    "Hip fracture",
    "Male"
  ),
  
  predict_fg_curve(
    fg_hip_male,
    male_dat_fg,
    "Hip fracture",
    "Male"
  )
)

# =========================================================
# Plot: Cox versus Fine-Gray spline comparison
# =========================================================

p_fg_vs_cox_hip <- ggplot2::ggplot(
  curve_dat,
  ggplot2::aes(
    x = log_MET_total,
    y = HR,
    linetype = Model
  )
) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = "dotted"
  ) +
  ggplot2::geom_line(
    linewidth = 1
  ) +
  ggplot2::facet_wrap(
    ~ Sex,
    ncol = 2
  ) +
  ggplot2::labs(
    x = "LTPA (log[MET-min/week + 1])",
    y = "Relative hazard",
    linetype = "Model",
    title = "Cox and Fine-Gray spline models for hip fracture",
    subtitle = "Fine-Gray models account for death as a competing event"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom",
    strip.background = ggplot2::element_rect(fill = "grey90"),
    panel.grid.minor = ggplot2::element_blank()
  )

print(p_fg_vs_cox_hip)


# =========================================================
# Save figure
# =========================================================

ggplot2::ggsave(
  filename = file.path(
    TTE_FIGURES_DIR,
    "Fig_S3_FineGray_vs_Cox_splines_hip.png"
  ),
  plot = p_fg_vs_cox_hip,
  width = 12,
  height = 5.5,
  dpi = 300
)

# =============================================================================
# Figure 4: Hip fracture intensity quintile forest plots by sex
# =============================================================================

fig4_dat <- readRDS(
  file.path(
    DATA_DERIVED,
    "res_hip_intensity_quintile.rds"
  )
)

# ---------------------------------------------------------
# Forest plots
# ---------------------------------------------------------

p_hip_intensity_female <- make_intensity_quintile_plot(
  fig4_dat$female,
  "Female"
)

p_hip_intensity_male <- make_intensity_quintile_plot(
  fig4_dat$male,
  "Male"
)

# ---------------------------------------------------------
# Supplementary table: Hazard ratios for intensity quintiles
# ---------------------------------------------------------

hip_intensity_HR <- dplyr::bind_rows(
  make_intensity_HR_table(fig4_dat$female),
  make_intensity_HR_table(fig4_dat$male)
)


# ---------------------------------------------------------
# Figure 4
# ---------------------------------------------------------

p_hip_intensity_female <- p_hip_intensity_female +
  labs(title = "Female")

p_hip_intensity_male <- p_hip_intensity_male +
  labs(title = "Male")

fig_4_hip_intensity_quintile_panel <-
  p_hip_intensity_female +
  p_hip_intensity_male +
  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

fig_4_hip_intensity_quintile_panel


# ---------------------------------------------------------
# Save figure
# ---------------------------------------------------------

ggsave(
  filename = file.path(
    TTE_FIGURES_DIR,
    "Fig_4_hip_intensity_quintile_panel.png"
  ),
  plot = fig_4_hip_intensity_quintile_panel,
  width = 14,
  height = 6,
  dpi = 300
)


# =========================================================
# Supplementary Table S7:
# Hazard ratios for hip fracture according to PA component quintiles
# =========================================================
#
# Built directly from the saved intensity-analysis results.
# p for trend is deliberately not reported because the exposure-response
# relationship was shown to be non-linear.
# =========================================================

fig4_dat <- readRDS(
  file.path(
    DATA_DERIVED,
    "res_hip_intensity_quintile.rds"
  )
)

hip_intensity_HR <- dplyr::bind_rows(
  fig4_dat$female,
  fig4_dat$male
) %>%
  dplyr::mutate(
    HR_CI = sprintf(
      "%.2f (%.2f–%.2f)",
      HR,
      CI_lower,
      CI_upper
    )
  ) %>%
  dplyr::select(
    Sex = sex,
    Exposure = exposure,
    quintile,
    HR_CI
  ) %>%
  tidyr::pivot_wider(
    names_from = quintile,
    values_from = HR_CI
  ) %>%
  dplyr::arrange(
    Sex,
    Exposure
  ) %>%
  dplyr::select(
    Sex,
    Exposure,
    Q2,
    Q3,
    Q4,
    Q5
  )

print(hip_intensity_HR)

save_table_word(
  df = hip_intensity_HR,
  table_number = "S7",
  folder_path = TTE_TABLES_DIR,
  title = paste(
    "Hazard ratios for incident hip fracture according to",
    "physical activity component quintiles"
  )
)


# =========================================================
# Supplementary Figure S4:
# Additional adjustment for musculoskeletal markers
# Total LTPA and hip fracture by sex
# =========================================================
#
# The figure appearance is retained, but it now reads from the single
# combined health-marker sensitivity result object.
# =========================================================

res_health_markers <- readRDS(
  file.path(
    DATA_DERIVED,
    "res_hip_health_marker_sensitivity.rds"
  )
)

health_hr <- res_health_markers$hr_results %>%
  dplyr::filter(
    exposure == "Total PA"
  )

common_y_scale <- ggplot2::scale_y_continuous(
  limits = c(0.5, 1.4),
  breaks = seq(0.5, 1.4, by = 0.1),
  minor_breaks = NULL
)

common_plot_theme <- ggplot2::theme(
  legend.position = "bottom"
)


# ---------------------------------------------------------
# Helper used only to reproduce the existing S4 layout
# ---------------------------------------------------------

make_marker_row <- function(
    marker_name,
    panel_label,
    adjusted_model_name
) {
  
  plot_dat <- health_hr %>%
    dplyr::filter(
      marker == marker_name
    )
  
  model_names <- c(
    "Fully adjusted",
    adjusted_model_name
  )
  
  model_cols <- c(
    "Fully adjusted" = "#0072B2",
    adjusted_model_name = "#D55E00"
  )
  
  names(model_cols)[2] <- adjusted_model_name
  
  p_female <- make_sensitivity_plot(
    tidy_df = plot_dat %>%
      dplyr::filter(
        sex == "Female"
      ),
    plot_title = paste0(
      panel_label,
      "\nFemale"
    ),
    model_names = model_names,
    model_cols = model_cols
  ) +
    common_y_scale
  
  p_male <- make_sensitivity_plot(
    tidy_df = plot_dat %>%
      dplyr::filter(
        sex == "Male"
      ),
    plot_title = "Male",
    model_names = model_names,
    model_cols = model_cols
  ) +
    common_y_scale
  
  (
    p_female +
      p_male +
      patchwork::plot_layout(
        ncol = 2,
        guides = "collect"
      )
  ) &
    common_plot_theme
}


row_falls <- make_marker_row(
  marker_name = "Falls",
  panel_label = "A. Falls",
  adjusted_model_name =
    "Fully adjusted + Falls"
)

row_grip <- make_marker_row(
  marker_name = "Grip strength",
  panel_label = "B. Grip strength",
  adjusted_model_name =
    "Fully adjusted + Grip strength"
)

row_heel_BMD <- make_marker_row(
  marker_name = "Heel BMD",
  panel_label = "C. Heel BMD",
  adjusted_model_name =
    "Fully adjusted + Heel BMD"
)

row_walking_pace <- make_marker_row(
  marker_name = "Walking pace",
  panel_label = "D. Walking pace",
  adjusted_model_name =
    "Fully adjusted + Walking pace"
)


fig_msk_marker_sensitivity <-
  row_falls /
  row_grip /
  row_heel_BMD /
  row_walking_pace +
  patchwork::plot_layout(
    heights = c(1, 1, 1, 1)
  )

fig_msk_marker_sensitivity

ggplot2::ggsave(
  filename = file.path(
    TTE_FIGURES_DIR,
    "Fig_S4_msk_marker_adjusted_sensitivity.png"
  ),
  plot = fig_msk_marker_sensitivity,
  width = 14,
  height = 22,
  dpi = 300
)


# =========================================================
# Supplementary Table S8:
# Health-marker sensitivity analyses across PA exposures
# =========================================================
#
# Exposures:
#   - Total PA
#   - Walking PA
#   - Moderate PA
#   - Vigorous PA
#
# Additional adjustment:
#   - Falls
#   - Grip strength
#   - Walking pace
#   - Heel BMD
#
# Each pair of models uses the same marker-complete sample.
# No p for trend is reported.
# =========================================================


tab_S8_long <- res_health_markers$hr_results %>%
  dplyr::mutate(
    Analysis = dplyr::case_when(
      model == "Fully adjusted" ~
        "Fully adjusted",
      grepl("^Fully adjusted \\+", model) ~
        "Additionally adjusted",
      TRUE ~ as.character(model)
    ),
    
    HR_CI = sprintf(
      "%.2f (%.2f–%.2f)",
      estimate,
      conf.low,
      conf.high
    )
  ) %>%
  dplyr::select(
    exposure,
    marker,
    sex,
    Analysis,
    quintile,
    HR_CI
  )


# ---------------------------------------------------------
# Put Q2-Q5 into columns
# ---------------------------------------------------------

tab_S8_wide <- tab_S8_long %>%
  tidyr::pivot_wider(
    names_from = quintile,
    values_from = HR_CI
  )


# ---------------------------------------------------------
# Put Female and Male alongside each other
# ---------------------------------------------------------

tab_S8 <- tab_S8_wide %>%
  tidyr::pivot_wider(
    names_from = sex,
    values_from = c(
      Q2,
      Q3,
      Q4,
      Q5
    ),
    names_glue = "{sex}_{.value}"
  ) %>%
  dplyr::mutate(
    
    Exposure = factor(
      exposure,
      levels = c(
        "Total PA",
        "Walking PA",
        "Moderate PA",
        "Vigorous PA"
      )
    ),
    
    Marker = factor(
      marker,
      levels = c(
        "Falls",
        "Grip strength",
        "Walking pace",
        "Heel BMD"
      )
    ),
    
    analysis_order = dplyr::if_else(
      Analysis == "Fully adjusted",
      1L,
      2L
    )
  ) %>%
  dplyr::arrange(
    Exposure,
    Marker,
    analysis_order
  ) %>%
  dplyr::select(
    Exposure,
    Marker,
    Analysis,
    
  
    `Female Q2` = Female_Q2,
    `Female Q3` = Female_Q3,
    `Female Q4` = Female_Q4,
    `Female Q5` = Female_Q5,
    

    `Male Q2` = Male_Q2,
    `Male Q3` = Male_Q3,
    `Male Q4` = Male_Q4,
    `Male Q5` = Male_Q5
  )

print(tab_S8)


save_table_word(
  df = tab_S8,
  table_number = "S8",
  folder_path = TTE_TABLES_DIR,
  title = paste(
    "Sensitivity analyses of physical activity and incident hip fracture",
    "with additional adjustment for falls, grip strength, walking pace",
    "and heel bone mineral density"
  )
)


# =========================================================
# Supplementary Table S9:
# Sample sizes for health-marker sensitivity analyses
# =========================================================
#
# This compact table makes the differing complete-case samples explicit.
# N and event counts are identical for the paired base and
# marker-adjusted models, so they are shown once per exposure/marker/sex.
# =========================================================

tab_S9 <- res_health_markers$sample_summary %>%
  dplyr::mutate(
    Exposure = factor(
      exposure,
      levels = c(
        "Total PA",
        "Walking PA",
        "Moderate PA",
        "Vigorous PA"
      )
    ),
    
    Marker = factor(
      marker,
      levels = c(
        "Falls",
        "Grip strength",
        "Walking pace",
        "Heel BMD"
      )
    ),
    
    Sex = factor(
      sex,
      levels = c(
        "Female",
        "Male"
      )
    ),
    
    N = format(
      n,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    
    Events = format(
      events,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )
  ) %>%
  dplyr::arrange(
    Exposure,
    Marker,
    Sex
  ) %>%
  dplyr::transmute(
    Exposure = as.character(Exposure),
    Marker = as.character(Marker),
    Sex = as.character(Sex),
    N,
    Events
  )

print(tab_S9)

save_table_word(
  df = tab_S9,
  table_number = "S9",
  folder_path = TTE_TABLES_DIR,
  title = paste(
    "Sample sizes and hip fracture events in physical activity",
    "health-marker sensitivity analyses"
  )
)

# =========================================================
# Figure: Walking-pace sensitivity across PA exposures
# Hip fracture by sex
# =========================================================

res_health_markers <- readRDS(
  file.path(
    DATA_DERIVED,
    "res_hip_health_marker_sensitivity.rds"
  )
)

walking_pace_dat <- res_health_markers$hr_results %>%
  dplyr::filter(
    marker == "Walking pace"
  )


# ---------------------------------------------------------
# Common settings
# ---------------------------------------------------------

walking_pace_model_names <- c(
  "Fully adjusted",
  "Fully adjusted + Walking pace"
)

walking_pace_model_cols <- c(
  "Fully adjusted" = "#0072B2",
  "Fully adjusted + Walking pace" = "#D55E00"
)

common_y_scale <- ggplot2::scale_y_continuous(
  limits = c(0.5, 1.4),
  breaks = seq(0.5, 1.4, by = 0.1),
  minor_breaks = NULL
)

common_plot_theme <- ggplot2::theme(
  legend.position = "bottom"
)


# ---------------------------------------------------------
# Helper to create one exposure row
# ---------------------------------------------------------

make_walking_pace_exposure_row <- function(
    exposure_name,
    panel_label,
    x_label
) {
  
  plot_dat <- walking_pace_dat %>%
    dplyr::filter(
      exposure == exposure_name
    )
  
  p_female <- make_sensitivity_plot(
    tidy_df = plot_dat %>%
      dplyr::filter(
        sex == "Female"
      ),
    plot_title = paste0(
      panel_label,
      "\nFemale"
    ),
    model_names = walking_pace_model_names,
    model_cols = walking_pace_model_cols
  ) +
    common_y_scale +
    ggplot2::labs(
      x = x_label
    )
  
  p_male <- make_sensitivity_plot(
    tidy_df = plot_dat %>%
      dplyr::filter(
        sex == "Male"
      ),
    plot_title = "Male",
    model_names = walking_pace_model_names,
    model_cols = walking_pace_model_cols
  ) +
    common_y_scale +
    ggplot2::labs(
      x = x_label
    )
  
  (
    p_female +
      p_male +
      patchwork::plot_layout(
        ncol = 2,
        guides = "collect"
      )
  ) &
    common_plot_theme
}


# ---------------------------------------------------------
# A. Total PA
# ---------------------------------------------------------

row_total <- make_walking_pace_exposure_row(
  exposure_name = "Total PA",
  panel_label = "A. Total PA",
  x_label = "Total PA quintile"
)


# ---------------------------------------------------------
# B. Walking PA
# ---------------------------------------------------------

row_walking <- make_walking_pace_exposure_row(
  exposure_name = "Walking PA",
  panel_label = "B. Walking PA",
  x_label = "Walking PA quintile"
)


# ---------------------------------------------------------
# C. Moderate PA
# ---------------------------------------------------------

row_moderate <- make_walking_pace_exposure_row(
  exposure_name = "Moderate PA",
  panel_label = "C. Moderate PA",
  x_label = "Moderate PA quintile"
)


# ---------------------------------------------------------
# D. Vigorous PA
# ---------------------------------------------------------

row_vigorous <- make_walking_pace_exposure_row(
  exposure_name = "Vigorous PA",
  panel_label = "D. Vigorous PA",
  x_label = "Vigorous PA quintile"
)


# ---------------------------------------------------------
# Combine rows
# ---------------------------------------------------------

fig_walking_pace_across_exposures <-
  row_total /
  row_walking /
  row_moderate /
  row_vigorous +
  patchwork::plot_layout(
    heights = c(1, 1, 1, 1)
  )

fig_walking_pace_across_exposures


# ---------------------------------------------------------
# Save figure
# ---------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    TTE_FIGURES_DIR,
    "Fig_walking_pace_sensitivity_across_PA_exposures.png"
  ),
  plot = fig_walking_pace_across_exposures,
  width = 14,
  height = 22,
  dpi = 300
)


