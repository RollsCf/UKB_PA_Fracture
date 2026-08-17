# To create the figures scripts 06-09 need to have been run first

source(here::here("scripts/00_setup.R"))
source(here("scripts/03_helpers_general.R"))
source(here("scripts/04_helpers_CXA.R"))
source(here("scripts/05_helpers_survival.R"))
source(here("scripts/06_helpers_figures_CXA.R"))


# =========================================================
# Figure 1: Exclusion flow table
# =========================================================

tab_exc <- readRDS(file.path(DATA_DERIVED, "tab_exclusions.rds"))

png(
  filename = file.path(FIGURES_DIR, "Figure_1_Exclusion_flow.png"),
  width = 1200,
  height = 800,
  res = 150
)

gridExtra::grid.table(tab_exc)

dev.off()

### As per the manuscript MET_Total will be known as LTPA volume or LTPAv
# Relabeled below.

# =========================================================
# Table 1: Baseline characteristics by sex
# Uses complete-case analytic cohort
# =========================================================

complete_case_dat <- readRDS(
  file.path(DATA_DERIVED, "complete_case_dat.rds")
)

tbl <- complete_case_dat %>% 
  dplyr::select(
    sex,
    age_A0,
    ethnicity_derived,
    tdi,
    education_level,
    weight,
    height,
    BMI,
    SRF,
    Wrist,
    Hip,
    cc_MET_mod_trunc,
    cc_MET_vig_trunc,
    cc_MET_walk_trunc,
    cc_MET_total_trunc
  ) %>% 
  gtsummary::tbl_summary(
    by = sex,
    statistic = list(
      c(age_A0, tdi, weight, height, BMI) ~ "{mean} ({sd})",
      c(cc_MET_mod_trunc, cc_MET_vig_trunc, cc_MET_walk_trunc, cc_MET_total_trunc) ~ "{median} ({p25}, {p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      c(age_A0, tdi, weight, height, BMI) ~ 1,
      c(cc_MET_mod_trunc, cc_MET_vig_trunc, cc_MET_walk_trunc, cc_MET_total_trunc) ~ 0
    ),
    type = gtsummary::all_categorical() ~ "categorical",
    label = list(
      age_A0 ~ "Age (years), mean (SD)",
      ethnicity_derived ~ "Ethnicity, n (%)",
      tdi ~ "Townsend deprivation index, mean (SD)",
      education_level ~ "Education level, n (%)",
      weight ~ "Weight (kg), mean (SD)",
      height ~ "Height (cm), mean (SD)",
      BMI ~ "BMI (kg/m²), mean (SD)",
      SRF ~ "Any fracture, n (%)",
      Wrist ~ "Wrist fracture, n (%)",
      Hip ~ "Hip fracture, n (%)",
      cc_MET_mod_trunc ~ "Moderate activity volume (MET-min/week), median (IQR)",
      cc_MET_vig_trunc ~ "Vigorous activity volume (MET-min/week), median (IQR)",
      cc_MET_walk_trunc ~ "Walking volume (MET-min/week), median (IQR)",
      cc_MET_total_trunc ~ "Leisure Time Physical Activity volume (LTPAv) (MET-min/week), median (IQR)"
    ),
    missing_text = "Missing"
  ) %>%
  gtsummary::add_overall(last = FALSE) %>%
  gtsummary::modify_header(label ~ "**Variable**") %>%
  gtsummary::modify_spanning_header(
    c("stat_0", "stat_1", "stat_2") ~ "**Sex**"
  ) %>%
  gtsummary::modify_caption(
    "**Table 1. Baseline characteristics of the analytic cohort overall and by sex**"
  ) %>%
  gtsummary::modify_footnote(
    gtsummary::all_stat_cols() ~ "Continuous demographic and anthropometric variables are presented as mean (SD). Physical activity variables are presented as median (IQR). Categorical variables are presented as n/N (%)."
  )

tbl

tbl %>%
  gtsummary::as_gt() %>%
  gt::gtsave(
    filename = file.path(TABLES_DIR, "Table1.docx")
  )

# =========================================================
# Supplementary Table 1 and 2: Response rate and degree of missing data
# =========================================================
# =========================================================
# Missingness and activity-data supplementary tables
# Uses FULL processed dataset, not complete-case dataset
# =========================================================

# Load full processed data
full_analysis_dat <- readRDS(file.path(DATA_DERIVED, "analysis_dat.Rds"))

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



vars <- c(
  "ethnicity_derived",
  "SRF",
  "num_days_mod",
  "num_day_walk",
  "num_days_vig",
  "cc_MET_mod_trunc",
  "cc_MET_vig_trunc",
  "cc_MET_walk_trunc",
  "cc_MET_total_trunc",
  "education_level",
  "tdi",
  "weight",
  "height",
  "BMI"
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
    dplyr::mutate(Variable = var) %>%
    dplyr::select(Variable, Answered, Answered_pct, Missing, Missing_pct)
}

Response_rates <- lapply(vars, function(v) response_rate(middle_aged_full_dat, v)) %>%
  dplyr::bind_rows()

Response_rates

save_table_word(
  df = Response_rates,
  table_number = "S1",
  folder_path = TABLES_DIR,
  title = "Response rate and degree of missing data"
)

##################################################################################################
# =========================================================
# Supplementary Table 2: With vs without activity data
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
  "height",
  "BMI",
  "SRF"
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
    height,
    BMI,
    SRF
  ) %>%
  gtsummary::tbl_summary(
    by = activity_data,
    statistic = list(
      c(age_A0, tdi, weight, height, BMI) ~ "{mean} ({sd})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      c(age_A0, tdi, weight, height, BMI) ~ 1
    ),
    type = gtsummary::all_categorical() ~ "categorical",
    label = list(
      sex ~ "Sex, n (%)",
      age_A0 ~ "Age (years), mean (SD)",
      ethnicity_derived ~ "Ethnicity, n (%)",
      tdi ~ "Townsend deprivation index, mean (SD)",
      education_level ~ "Education level, n (%)",
      weight ~ "Weight (kg), mean (SD)",
      height ~ "Height (cm), mean (SD)",
      BMI ~ "BMI (kg/m²), mean (SD)",
      SRF ~ "Any fracture, n (%)"
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
  table_number = "S2",
  folder_path = TABLES_DIR,
  title = "Missing data comparison"
)

# =========================================================
# Supplementary Figure 2: Histogram of PA and log transformed PA
# =========================================================


met_raw_vars <- c(
  "MET_mod",
  "MET_vig",
  "MET_walk",
  "MET_total"
)

met_log_vars <- c(
  "MET_mod_log",
  "MET_vig_log",
  "MET_walk_log",
  "MET_total_log"
)

met_labels <- c(
  MET_mod = "Moderate",
  MET_vig = "Vigorous",
  MET_walk = "Walking",
  MET_total = "LTPAv"
)

hist_raw_df <- IPAQ_dat2 %>%
  dplyr::select(dplyr::all_of(met_raw_vars)) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    variable = factor(variable, levels = met_raw_vars, labels = met_labels[met_raw_vars])
  )

hist_log_df <- IPAQ_dat2 %>%
  dplyr::select(dplyr::all_of(met_log_vars)) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    variable = sub("_log$", "", variable),
    variable = factor(variable, levels = met_raw_vars, labels = met_labels[met_raw_vars])
  )

p_raw <- ggplot2::ggplot(hist_raw_df, ggplot2::aes(x = value)) +
  ggplot2::geom_histogram(bins = 50, colour = "white") +
  ggplot2::facet_wrap(~ variable, scales = "free", ncol = 2) +
  ggplot2::labs(
    title = "Before log + 1 transformation",
    x = "MET-min/week",
    y = "Count"
  ) +
  ggplot2::theme_bw()

p_log <- ggplot2::ggplot(hist_log_df, ggplot2::aes(x = value)) +
  ggplot2::geom_histogram(bins = 50, colour = "white") +
  ggplot2::facet_wrap(~ variable, scales = "free", ncol = 2) +
  ggplot2::labs(
    title = "After log + 1 transformation",
    x = "log(MET-min/week + 1)",
    y = "Count"
  ) +
  ggplot2::theme_bw()

p_supp_s2 <- p_raw / p_log +
  patchwork::plot_annotation(
    
  )

p_supp_s2

ggplot2::ggsave(
  filename = file.path(FIGURES_DIR, "Fig_S2_MET_histograms_raw_log.png"),
  plot = p_supp_s2,
  width = 12,
  height = 10,
  dpi = 300
)

# =========================================================
# Supplementary Figure 3 Spline
# =========================================================

p_spline_LTPAv <- plot_spline_or(
  data = IPAQ_dat2,
  outcome = "SRF",
  exposure_log = "MET_total_log",
  confounders = full_adjustment,
  spline_knots = 4
) +
  ggplot2::labs(
    title = NULL,
    x = "Log-transformed LTPA volume [log(MET-min/week + 1)]"
  )

ggplot2::ggsave(
  filename = file.path(
    FIGURES_DIR,
    "Supplementary_Figure_S3_LTPAv_Spline.png"
  ),
  plot = p_spline_LTPAv,
  width = 10,
  height = 6,
  dpi = 600,
  bg = "white"
)


# # =======================================================================
# Supplementary Table S3: Fracture site distribution by age group and sex
# =======================================================================


IPAQ_dat2 <- readRDS(file.path(DATA_DERIVED, "IPAQ_dat2.rds"))

fracture_vars <- c("Hip", "Wrist", "Leg", "Arm", "Ankle", "Spine", "Other")

fracture_long <- IPAQ_dat2 %>%
  dplyr::select(
    sex,
    agegp_A0,
    MET_total_q,
    dplyr::all_of(fracture_vars)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(fracture_vars),
    names_to = "fracture_type",
    values_to = "fracture"
  ) %>%
  dplyr::filter(fracture == "Yes")



fracture_summary_tbl <- fracture_long %>%
  dplyr::count(sex, agegp_A0, fracture_type) %>%
  dplyr::group_by(sex, agegp_A0) %>%
  dplyr::mutate(
    value = paste0(
      n,
      " (",
      round(100 * n / sum(n), 1),
      "%)"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    sex,
    agegp_A0,
    fracture_type,
    value
  ) %>%
  tidyr::pivot_wider(
    names_from = fracture_type,
    values_from = value
  ) %>%
  dplyr::arrange(sex, agegp_A0)

save_table_word(
  df = fracture_summary_tbl,
  table_number = "S3",
  folder_path = TABLES_DIR,
  title = "S3_Distribution of fracture sites by age group and sex"
)

# ===============================================================
# Supplementary Table S4: Fracture site by activity quintile, 
# ===============================================================

total_fracture_by_quintile <- IPAQ_dat2 %>%
  dplyr::filter(
    !is.na(sex),
    !is.na(MET_total_q),
    !is.na(SRF)
  ) %>%
  dplyr::group_by(sex, MET_total_q) %>%
  dplyr::summarise(
    N = dplyr::n(),
    total_fractures = sum(SRF == "Yes", na.rm = TRUE),
    fractures_per_100 = 100 * total_fractures / N,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    `Total fracture, n (%)` = paste0(
      total_fractures,
      " (",
      sprintf("%.1f", fractures_per_100),
      "%)"
    )
  ) %>%
  dplyr::select(
    sex,
    quintile = MET_total_q,
    N,
    `Total fracture, n (%)`
  ) %>%
  dplyr::arrange(sex, quintile)

save_table_word(
  df = total_fracture_by_quintile,
  table_number = "S4",
  folder_path = TABLES_DIR,
  title = "Total fractures per 100 participants by LTPA quintile"
)


# =========================================================
#  Supplementary Table 5 Activity by quintile
# =========================================================


IPAQ_dat2_prop <- readRDS(file.path(DATA_DERIVED, "IPAQ_dat2_prop.rds"))

tbl_activity_profiles_all_quintiles <- IPAQ_dat2_prop %>%
  dplyr::group_by(agegp_A0, sex, MET_total_q) %>%
  dplyr::summarise(
    n = dplyr::n(),
    
    MET_total = paste0(
      round(stats::median(MET_total, na.rm = TRUE), 0),
      " (",
      round(stats::quantile(MET_total, 0.25, na.rm = TRUE), 0),
      "–",
      round(stats::quantile(MET_total, 0.75, na.rm = TRUE), 0),
      ")"
    ),
    
    mins_wk_total = paste0(
      round(stats::median(mins_wk_Total, na.rm = TRUE), 0),
      " (",
      round(stats::quantile(mins_wk_Total, 0.25, na.rm = TRUE), 0),
      "–",
      round(stats::quantile(mins_wk_Total, 0.75, na.rm = TRUE), 0),
      ")"
    ),
    
    walking_min_week = paste0(
      round(stats::median(mins_wk_walk, na.rm = TRUE), 0),
      " (",
      round(stats::quantile(mins_wk_walk, 0.25, na.rm = TRUE), 0),
      "–",
      round(stats::quantile(mins_wk_walk, 0.75, na.rm = TRUE), 0),
      ")"
    ),
    
    moderate_min_week = paste0(
      round(stats::median(mins_wk_mod, na.rm = TRUE), 0),
      " (",
      round(stats::quantile(mins_wk_mod, 0.25, na.rm = TRUE), 0),
      "–",
      round(stats::quantile(mins_wk_mod, 0.75, na.rm = TRUE), 0),
      ")"
    ),
    
    vigorous_min_week = paste0(
      round(stats::median(mins_wk_vig, na.rm = TRUE), 0),
      " (",
      round(stats::quantile(mins_wk_vig, 0.25, na.rm = TRUE), 0),
      "–",
      round(stats::quantile(mins_wk_vig, 0.75, na.rm = TRUE), 0),
      ")"
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::rename(
    age_group = agegp_A0,
    quintile = MET_total_q
  ) %>%
  dplyr::arrange(age_group, sex, quintile)

tbl_activity_profiles_gt <- tbl_activity_profiles_all_quintiles %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Activity profiles across physical activity quintiles by sex and age group",
    subtitle = "Values are median (IQR)"
  ) %>%
  gt::cols_label(
    age_group = "Age group",
    sex = "Sex",
    quintile = "PA quintile",
    n = "N",
    MET_total = "LTPAv MET-min/week",
    mins_wk_total = "LTPA min/week",
    walking_min_week = "Walking min/week",
    moderate_min_week = "Moderate min/week",
    vigorous_min_week = "Vigorous min/week"
  )


tbl_activity_profiles_save <- tbl_activity_profiles_all_quintiles %>%
  dplyr::rename(
    `Age group` = age_group,
    Sex = sex,
    `PA quintile` = quintile,
    N = n,
    `LTPAv MET-min/week` = MET_total,
    `LTPA min/week` = mins_wk_total,
    `Walking min/week` = walking_min_week,
    `Moderate min/week` = moderate_min_week,
    `Vigorous min/week` = vigorous_min_week
  )

save_table_word(
  df = tbl_activity_profiles_save,
  table_number = "S5",
  folder_path = TABLES_DIR,
  title = "Supplementary_Table_activity_profiles_all_quintiles_by_sex.docx"
)



# =====================================================================================
# Fracture site by activity quintile, age group and sex (To use with figure later)
# =====================================================================================

fracture_summary_q_age <- fracture_long %>%
  dplyr::count(sex, agegp_A0, MET_total_q, fracture_type) %>%
  dplyr::group_by(sex, agegp_A0, MET_total_q) %>%
  dplyr::mutate(
    prop = n / sum(n)
  ) %>%
  dplyr::ungroup()

fracture_order <- fracture_summary_q_age %>%
  dplyr::group_by(fracture_type) %>%
  dplyr::summarise(total_n = sum(n), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total_n)) %>%
  dplyr::pull(fracture_type)

fracture_summary_q_age <- fracture_summary_q_age %>%
  dplyr::mutate(
    agegp_A0 = factor(agegp_A0, levels = c("40-49", "50-59", "60-69")),
    MET_total_q = factor(MET_total_q, levels = paste0("Q", 1:5)),
    fracture_type = factor(fracture_type, levels = rev(fracture_order))
  )

plot_fracture_site_q_age <- function(data, plot_title) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = MET_total_q,
      y = prop,
      fill = fracture_type
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ agegp_A0, nrow = 1) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
      x = "Total activity quintile",
      y = "Percentage of reported fracture sites",
      fill = "Fracture site",
      title = plot_title
    ) +
    ggplot2::theme_bw()
}

p_fracture_female <- plot_fracture_site_q_age(
  data = dplyr::filter(fracture_summary_q_age, sex == "Female"),
  plot_title = "Females: fracture site distribution by age group and activity quintile"
)

p_fracture_male <- plot_fracture_site_q_age(
  data = dplyr::filter(fracture_summary_q_age, sex == "Male"),
  plot_title = "Males: fracture site distribution by age group and activity quintile"
)


ggplot2::ggsave(
  filename = file.path(FIGURES_DIR, "Fig_fracture_site_by_age_activity_quintile_female.png"),
  plot = p_fracture_female,
  width = 11,
  height = 5,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(FIGURES_DIR, "Fig_fracture_site_by_age_activity_quintile_male.png"),
  plot = p_fracture_male,
  width = 11,
  height = 5,
  dpi = 300
)

# =========================================================
# Supplementary table S6: Assessment of linearity
# =========================================================

linearity_results <- readRDS(
  file.path(DATA_DERIVED, "linearity_results.rds")
)

save_table_word(
  df = linearity_results,
  table_number = "S6",
  folder_path = TABLES_DIR,
  title = "Assessment of linearity for continuous physical activity measures (likelihood ratio test comparing linear and spline models)"
)

# =========================================================
# Supplementary Table SX: Alternative Exposure with less missingness
# =========================================================

# sensitivity_quintiles_age_sex <- readRDS(
#   file.path(DATA_DERIVED, "sensitivity_quintiles_age_sex.rds")
# )
# 
# save_table_word(
#   df = sensitivity_quintiles_age_sex,
#   table_number = "SX",
#   folder_path = TABLES_DIR,
#   title = "Sensitivity odds ratios for raw PA quintiles across three models for SRF, stratified by sex and age group"
# )

# ==================================================================
# Supplementary Figure S4: Alternative Exposure with less missingness
# ==================================================================

# Create figure from saved results
model_levels <- c(
  "Unadjusted",
  "Minimal adjustment",
  "Fully adjusted"
)

plot_df_age_sex_total <- prep_quintile_plot_df(
  df = sensitivity_quintiles_age_sex,
  facet_vars = c("sex", "agegp_A0"),
  exposure_filter = "num_days_total",
  outcome_filter = "SRF"
) %>%
  dplyr::filter(agegp_A0 %in% c("40-49", "50-59", "60-69")) %>%
  dplyr::mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    agegp_A0 = factor(agegp_A0, levels = c("40-49", "50-59", "60-69")),
    model = factor(model, levels = model_levels)
  )

p_age_sex_total <- plot_quintile_or(
  data = plot_df_age_sex_total,
  facet_formula = sex ~ agegp_A0,
  title_text = "Odds ratios for Num days total quintiles and SRF by sex, age and model"
)

save_quintile_plot(
  plot_obj = p_age_sex_total,
  filename = "Fig_S4_Sensitivity_numdays_total_age_sex.png",
  figures_dir = FIGURES_DIR
)


# =========================================================
# Supplementary Table S7: Interaction testing for sex and age
# =========================================================

interaction_results <- readRDS(
  file.path(DATA_DERIVED, "interaction_results.rds")
)

save_table_word(
  df = interaction_results,
  table_number = "S7",
  folder_path = TABLES_DIR,
  title = "Assessment of sex and age-group interaction for spline physical activity models"
)

# =======================================================================================
# Figure 2: Fracture-site distribution and activity profiles
# Panel A: Fracture-site heatmap
# Panel B: Absolute activity profiles
# =======================================================================================

# =========================================================
# Load data
# =========================================================

IPAQ_dat2 <- readRDS(
  file.path(DATA_DERIVED, "IPAQ_dat2.rds")
)

# If IPAQ_dat2_prop is saved separately, load it here.
# Otherwise, keep your existing code that creates IPAQ_dat2_prop.

# =========================================================
# Standardise age-group labels
# Important: recode before converting to factor
# =========================================================

age_levels <- c(
  "40-49",
  "50-59",
  "60-65"
)

sex_levels <- c(
  "Female",
  "Male"
)

IPAQ_dat2 <- IPAQ_dat2 %>%
  dplyr::mutate(
    agegp_A0 = dplyr::recode(
      as.character(agegp_A0),
      "60-69" = "60-65"
    ),
    agegp_A0 = factor(
      agegp_A0,
      levels = age_levels
    ),
    sex = factor(
      as.character(sex),
      levels = sex_levels
    )
  )

IPAQ_dat2_prop <- IPAQ_dat2_prop %>%
  dplyr::mutate(
    agegp_A0 = dplyr::recode(
      as.character(agegp_A0),
      "60-69" = "60-65"
    ),
    agegp_A0 = factor(
      agegp_A0,
      levels = age_levels
    ),
    sex = factor(
      as.character(sex),
      levels = sex_levels
    )
  )

# Optional checks
table(IPAQ_dat2$agegp_A0, useNA = "ifany")
table(IPAQ_dat2_prop$agegp_A0, useNA = "ifany")



# =======================================================================================
# Panel A: Fracture-site heatmap
# =======================================================================================

fracture_vars <- c(
  "Hip",
  "Spine",
  "Wrist",
  "Arm",
  "Leg",
  "Ankle",
  "Other"
)

sex_age_levels <- c(
  "Female 40-49",
  "Female 50-59",
  "Female 60-65",
  "Male 40-49",
  "Male 50-59",
  "Male 60-65"
)

fracture_type_levels <- c(
  "Other",
  "Wrist",
  "Ankle",
  "Arm",
  "Leg",
  "Hip",
  "Spine"
)

# =========================================================
# Full cohort N for each sex-age group
# =========================================================

cohort_n_age_sex <- IPAQ_dat2 %>%
  dplyr::filter(
    !is.na(sex),
    !is.na(agegp_A0)
  ) %>%
  dplyr::mutate(
    sex_age = paste(
      sex,
      agegp_A0
    ),
    sex_age = factor(
      sex_age,
      levels = sex_age_levels
    )
  ) %>%
  dplyr::count(
    sex_age,
    name = "cohort_N"
  ) %>%
  tidyr::complete(
    sex_age = factor(
      sex_age_levels,
      levels = sex_age_levels
    ),
    fill = list(
      cohort_N = 0
    )
  )

x_labels <- stats::setNames(
  paste0(
    as.character(cohort_n_age_sex$sex_age),
    "\n(N = ",
    format(
      cohort_n_age_sex$cohort_N,
      big.mark = ",",
      scientific = FALSE
    ),
    ")"
  ),
  as.character(cohort_n_age_sex$sex_age)
)

# =========================================================
# Reshape fracture-site data
# =========================================================

fracture_long_age_sex <- IPAQ_dat2 %>%
  dplyr::filter(
    !is.na(sex),
    !is.na(agegp_A0)
  ) %>%
  dplyr::select(
    sex,
    agegp_A0,
    dplyr::all_of(fracture_vars)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(fracture_vars),
    names_to = "fracture_type",
    values_to = "fracture"
  ) %>%
  dplyr::filter(
    fracture == "Yes"
  ) %>%
  dplyr::mutate(
    sex_age = paste(
      sex,
      agegp_A0
    ),
    sex_age = factor(
      sex_age,
      levels = sex_age_levels
    ),
    fracture_type = factor(
      fracture_type,
      levels = fracture_type_levels
    )
  )

# =========================================================
# Calculate percentage distribution of fracture sites
# =========================================================

fracture_heatmap_age_sex_df <- fracture_long_age_sex %>%
  dplyr::count(
    sex_age,
    fracture_type,
    name = "n"
  ) %>%
  tidyr::complete(
    sex_age = factor(
      sex_age_levels,
      levels = sex_age_levels
    ),
    fracture_type = factor(
      fracture_type_levels,
      levels = fracture_type_levels
    ),
    fill = list(
      n = 0
    )
  ) %>%
  dplyr::group_by(
    sex_age
  ) %>%
  dplyr::mutate(
    total_reported_fracture_sites = sum(n),
    prop = dplyr::if_else(
      total_reported_fracture_sites > 0,
      n / total_reported_fracture_sites,
      NA_real_
    ),
    prop_pct = 100 * prop,
    label = dplyr::if_else(
      is.na(prop_pct),
      "",
      paste0(
        round(prop_pct, 1),
        "%"
      )
    ),
    text_colour = dplyr::if_else(
      !is.na(prop_pct) & prop_pct >= 30,
      "white",
      "black"
    )
  ) %>%
  dplyr::ungroup()

# =========================================================
# Create heatmap
# =========================================================

p_heatmap_panel <- ggplot2::ggplot(
  fracture_heatmap_age_sex_df,
  ggplot2::aes(
    x = sex_age,
    y = fracture_type,
    fill = prop_pct
  )
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = label,
      colour = text_colour
    ),
    size = 5.5,
    fontface = "bold"
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_x_discrete(
    labels = x_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_discrete(
    limits = rev,
    drop = FALSE
  ) +
  ggplot2::scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = "grey95",
    labels = function(x) paste0(x, "%")
  ) +
  ggplot2::labs(
    title = "Fracture site distribution by sex and age group",
    subtitle = paste(
      "Cells show the percentage of reported fracture sites",
      "within each sex and age group"
    ),
    x = NULL,
    y = "Fracture site",
    fill = "% of fracture sites",
    caption = paste(
      "N denotes the number of participants in each sex-age group.",
      "Percentages are calculated among reported fracture sites."
    )
  ) +
  ggplot2::theme_bw(
    base_size = 18
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      size = 20
    ),
    plot.subtitle = ggplot2::element_text(
      size = 16,
      hjust = 0.5
    ),
    plot.caption = ggplot2::element_text(
      size = 12,
      hjust = 0
    ),
    axis.title.y = ggplot2::element_text(
      size = 18
    ),
    axis.text.x = ggplot2::element_text(
      size = 14,
      angle = 35,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_text(
      size = 16
    ),
    legend.title = ggplot2::element_text(
      size = 16
    ),
    legend.text = ggplot2::element_text(
      size = 14
    ),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  )

# =======================================================================================
# Panel B: Activity profiles
# =======================================================================================

plot_absolute_activity <- prep_absolute_activity_profile_df(
  IPAQ_dat2_prop
) %>%
  dplyr::mutate(
    agegp_A0 = dplyr::recode(
      as.character(agegp_A0),
      "60-69" = "60-65"
    ),
    agegp_A0 = factor(
      agegp_A0,
      levels = age_levels
    ),
    sex = factor(
      as.character(sex),
      levels = sex_levels
    )
  ) %>%
  dplyr::filter(
    !is.na(agegp_A0),
    !is.na(sex)
  )

# Check there are no NA age facets
table(plot_absolute_activity$agegp_A0, useNA = "ifany")

p_activity_profiles <- ggplot2::ggplot(
  plot_absolute_activity,
  ggplot2::aes(
    y = forcats::fct_rev(MET_total_q),
    x = minutes,
    fill = intensity
  )
) +
  ggplot2::geom_col(
    width = 0.7,
    position = ggplot2::position_stack(
      reverse = TRUE
    )
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(agegp_A0),
    cols = ggplot2::vars(sex),
    drop = FALSE,
    labeller = ggplot2::labeller(
      agegp_A0 = function(x) paste0("Age ", x)
    )
  ) +
  ggplot2::scale_fill_manual(
    values = activity_cols,
    name = "Activity type"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1800),
    breaks = seq(
      0,
      1800,
      by = 300
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.02)
    )
  ) +
  ggplot2::labs(
    title = "Median amount of activity by intensity, age and sex",
    x = "Median minutes/week",
    y = "Activity quintile"
  ) +
  ggplot2::theme_minimal(
    base_size = 18
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 20,
      hjust = 0.5
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 16
    ),
    axis.title = ggplot2::element_text(
      size = 17
    ),
    axis.text = ggplot2::element_text(
      size = 15
    ),
    legend.title = ggplot2::element_text(
      size = 16
    ),
    legend.text = ggplot2::element_text(
      size = 14
    )
  )

# =======================================================================================
# Combine panels into Figure 2
# =======================================================================================

p_figure2 <-
  p_heatmap_panel /
  p_activity_profiles +
  patchwork::plot_layout(
    heights = c(1, 1.4)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 22
      )
    )
  )

p_figure2

# =======================================================================================
# Save Figure 2
# =======================================================================================

ggplot2::ggsave(
  filename = file.path(
    FIGURES_DIR,
    "Figure_2_FractureHeatmap_and_ActivityProfiles.png"
  ),
  plot = p_figure2,
  width = 14,
  height = 16,
  dpi = 600,
  bg = "white"
)

# =========================================================
# Figure 3: LTPAv and fracture risk by sex and age group
# =========================================================

# Calculate age/sex N
age_sex_n <- IPAQ_dat2 %>%
  dplyr::mutate(
    agegp_A0 = dplyr::recode(
      as.character(agegp_A0),
      "60-69" = "60-65"
    ),
    agegp_A0 = factor(
      agegp_A0,
      levels = c("40-49", "50-59", "60-65")
    ),
    sex = factor(
      as.character(sex),
      levels = c("Female", "Male")
    )
  ) %>%
  dplyr::filter(
    !is.na(sex),
    !is.na(agegp_A0)
  ) %>%
  dplyr::count(
    sex,
    agegp_A0,
    name = "N"
  )

plot_df_age_sex_total_no_ref <- plot_df_age_sex_total %>%
  dplyr::filter(Quintile != "Q1") %>%
  dplyr::mutate(
    Quintile = factor(
      as.character(Quintile),
      levels = c("Q2", "Q3", "Q4", "Q5")
    )
  )

p_or_female <- plot_total_or_age(
  data = plot_df_age_sex_total_no_ref,
  sex_value = "Female",
  title_text = "Female",
  age_labels = female_age_labels,
  y_limits = c(0.75, 2.1)
) +
  ggplot2::labs(
    tag = "A"
  )

p_or_male <- plot_total_or_age(
  data = plot_df_age_sex_total_no_ref,
  sex_value = "Male",
  title_text = "Male",
  age_labels = male_age_labels,
  y_limits = c(0.75, 2.1)
) +
  ggplot2::labs(
    tag = "B"
  )

p_figure3_forest <- patchwork::wrap_plots(
  p_or_female,
  p_or_male,
  ncol = 2,
  widths = c(1, 1),
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = paste(
      "Association between leisure time physical activity volume",
      "(MET-min/week) and fracture risk by sex and age group"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 22,
        hjust = 0.5
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 22
      )
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

p_figure3_forest

ggplot2::ggsave(
  filename = file.path(
    FIGURES_DIR,
    "Figure_3_LTPAv_FractureRisk_by_Sex_Age.png"
  ),
  plot = p_figure3_forest,
  width = 22,
  height = 8.5,
  dpi = 600,
  bg = "white"
)


# =========================================================
# Supplementary Table S8 OR for age and sex stratified analysis
# =========================================================


results_quintiles_age_sex <- readRDS(
  file.path(DATA_DERIVED, "results_quintiles_age_sex.rds")
)

results_quintiles_age_sex_table <- results_quintiles_age_sex %>%
  dplyr::mutate(
    Quintile = stringr::str_extract(term, "Q[2-5]"),
    Quintile = stringr::str_remove(Quintile, "Q"),
    Odds_Ratio = round(Odds_Ratio, 2),
    CI_lower = round(CI_lower, 2),
    CI_upper = round(CI_upper, 2),
    p_value = dplyr::case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    Global_p_value = dplyr::case_when(
      is.na(Global_p_value) ~ NA_character_,
      Global_p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", Global_p_value)
    )
  ) %>%
  dplyr::select(
    Sex = sex,
    `Age Group` = agegp_A0,
    Model = model,
    Quintile,
    `Odds ratio` = Odds_Ratio,
    `CI lower` = CI_lower,
    `CI upper` = CI_upper,
    `P value` = p_value,
    `Global P value` = Global_p_value
  )

save_table_word(
  df = results_quintiles_age_sex_table,
  table_number = "S8",
  folder_path = TABLES_DIR,
  title = "Odds ratios for LTPAv quintiles and self-reported fracture, stratified by sex and age group"
)

# =========================================================
# Supplementary figure S5 Weight/height sensitivity
# =========================================================


results_weight_sensitivity <- readRDS(
  file.path(DATA_DERIVED, "results_weight_sensitivity.rds")
)

save_table_word(
  df = results_weight_sensitivity,
  table_number = "weight_height_sensitivity",
  folder_path = TABLES_DIR,
  title = "Sensitivity analysis using alternative body size and adiposity adjustments, fully adjusted models"
)

# =========================================================
# Plot 
# =========================================================

plot_df_weight_sens <- prep_quintile_plot_df(
  df = results_weight_sensitivity,
  facet_vars = c("adiposity_adjustment"),
  exposure_filter = "MET_total",
  outcome_filter = "SRF"
) %>%
  dplyr::mutate(
    adiposity_adjustment = factor(
      adiposity_adjustment,
      levels = c(
        "Weight + height",
        "BMI",
        "Waist circumference + height"
      )
    )
  )

p_weight_sens <- ggplot2::ggplot(
  plot_df_weight_sens,
  ggplot2::aes(
    x = Quintile,
    y = Odds_Ratio,
    group = adiposity_adjustment
  )
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "black"
  ) +
  ggplot2::facet_grid(exposure ~ adiposity_adjustment) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = "Activity quintile",
    y = "Odds ratio (log scale, 95% CI)",
    title = "Fully adjusted sensitivity analysis using alternative adiposity adjustments"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "lightgrey", colour = "black"),
    strip.text = ggplot2::element_text(face = "bold")
  )

print(p_weight_sens)

save_quintile_plot(
  plot_obj = p_weight_sens,
  filename = "Fig_S5_adiposity_sensitivity_fully_adjusted.png",
  figures_dir = FIGURES_DIR
)

# =======================================================================================
# Figure 4: Mutually adjusted activity intensities and fracture risk
# =======================================================================================

results_mutual_age_sex <- readRDS(
  file.path(
    DATA_DERIVED,
    "results_mutual_quintiles_age_sex.rds"
  )
)

results_mutual_age_sex_plot_df <- prep_mutual_age_sex_plot_df(
  results_mutual_age_sex
) %>%
  dplyr::filter(Quintile != "Q1") %>%
  dplyr::mutate(
    Quintile = factor(
      Quintile,
      levels = c("Q2", "Q3", "Q4", "Q5")
    )
  )

p_mutual_female <- plot_mutual_or_age(
  data = results_mutual_age_sex_plot_df,
  sex_value = "Female",
  title_text = "Female",
  age_labels = female_age_labels,
  y_limits = c(0.8, 2.0)
) +
  ggplot2::labs(
    tag = "A"
  )

p_mutual_male <- plot_mutual_or_age(
  data = results_mutual_age_sex_plot_df,
  sex_value = "Male",
  title_text = "Male",
  age_labels = male_age_labels,
  y_limits = c(0.8, 2.0)
) +
  ggplot2::labs(
    tag = "B"
  )

p_figure4_mutual <- patchwork::wrap_plots(
  p_mutual_female,
  p_mutual_male,
  ncol = 2,
  widths = c(1, 1),
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = paste(
      "Mutually adjusted associations between activity intensity",
      "quintiles and fracture risk"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 22,
        hjust = 0.5
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 22
      )
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

p_figure4_mutual

ggplot2::ggsave(
  filename = file.path(
    FIGURES_DIR,
    "Figure_4_MutuallyAdjusted_Activity_FractureRisk.png"
  ),
  plot = p_figure4_mutual,
  width = 22,
  height = 8.5,
  dpi = 600,
  bg = "white"
)



# =========================================================
# Supp Table S9
# =========================================================

table_mutual <- results_mutual_age_sex_plot_df %>%
  dplyr::filter(Quintile != "Q1") %>%
  dplyr::mutate(
    Estimate = paste0(
      round(OR, 2),
      " (",
      round(CI_lower, 2),
      "–",
      round(CI_upper, 2),
      ")"
    )
  ) %>%
  dplyr::select(
    Sex = sex,
    Age_group = agegp_A0,
    Exposure = exposure,
    Quintile,
    Estimate
  )

table_mutual_wide <- table_mutual %>%
  tidyr::pivot_wider(
    names_from = Quintile,
    values_from = Estimate
  )

save_table_word(
  df = table_mutual_wide,
  table_number = "S9",
  folder_path = TABLES_DIR,
  title = "Mutually_Adjusted_Activity_ORs"
)


# =========================================================
# Supp Table S8 Proportion models within MET_total quintiles
# =========================================================

results_prop_within_q_df <- readRDS(
  file.path(DATA_DERIVED, "results_prop_quintile_activity.rds")
)

results_prop_within_q_table <- results_prop_within_q_df %>%
  dplyr::select(-term) %>%
  dplyr::relocate(model, sex) %>%
  dplyr::mutate(
    OR = round(OR, 2),
    CI_lower = round(CI_lower, 2),
    CI_upper = round(CI_upper, 2),
    p_value = dplyr::case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.2f", p_value)
    )
  )

save_table_word(
  df = results_prop_within_q_table,
  table_number = "S10",
  folder_path = TABLES_DIR,
  title = "Results_Prop_Quintile_Activity"
)


# =====================================================================
# Figure 4: Activity composition and fracture
# =====================================================================

results_prop_within_q_df <- readRDS(
  file.path(DATA_DERIVED, "results_prop_quintile_activity.rds")
)

# -----------------------------
# Panel A: OR per 10% higher activity proportion
# -----------------------------

prop_q_plot_df <- results_prop_within_q_df %>%
  dplyr::mutate(
    MET_total_q = stringr::str_extract(model, "Q[1-5]"),
    activity_prop = dplyr::case_when(
      term == "prop_walk_10" ~ "Walking",
      term == "prop_mod_10"  ~ "Moderate",
      term == "prop_vig_10"  ~ "Vigorous",
      TRUE ~ term
    ),
    MET_total_q = factor(MET_total_q, levels = paste0("Q", 1:5)),
    activity_prop = factor(activity_prop, levels = c("Walking", "Moderate", "Vigorous")),
    sex = factor(sex, levels = c("Female", "Male"))
  )

p_prop_forest <- plot_prop_within_q(
  data = prop_q_plot_df,
  title_text = "Associations per 10% higher activity proportion",
  y_limits = c(0.85, 1.15)
) +
  ggplot2::facet_grid(. ~ sex) +
  ggplot2::theme(
    legend.position = "none"
  )

# -----------------------------
# Panel B: Percentage composition of activity by total activity quintile
# -----------------------------

activity_prop_df <- IPAQ_dat2_prop %>%
  dplyr::mutate(
    total_mins = mins_wk_walk + mins_wk_mod + mins_wk_vig
  ) %>%
  dplyr::filter(
    !is.na(total_mins),
    total_mins > 0,
    !is.na(mins_wk_walk),
    !is.na(mins_wk_mod),
    !is.na(mins_wk_vig)
  ) %>%
  dplyr::mutate(
    Walking  = mins_wk_walk / total_mins,
    Moderate = mins_wk_mod  / total_mins,
    Vigorous = mins_wk_vig  / total_mins
  ) %>%
  dplyr::group_by(sex, MET_total_q) %>%
  dplyr::summarise(
    Walking  = mean(Walking, na.rm = TRUE),
    Moderate = mean(Moderate, na.rm = TRUE),
    Vigorous = mean(Vigorous, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    total_prop = Walking + Moderate + Vigorous,
    Walking  = Walking  / total_prop * 100,
    Moderate = Moderate / total_prop * 100,
    Vigorous = Vigorous / total_prop * 100
  ) %>%
  tidyr::pivot_longer(
    cols = c(Walking, Moderate, Vigorous),
    names_to = "activity_prop",
    values_to = "prop_pct"
  ) %>%
  dplyr::mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    MET_total_q = factor(MET_total_q, levels = paste0("Q", 1:5)),
    activity_prop = factor(activity_prop, levels = c("Walking", "Moderate", "Vigorous"))
  )

p_activity_prop_stack <- ggplot2::ggplot(
  activity_prop_df,
  ggplot2::aes(
    x = MET_total_q,
    y = prop_pct,
    fill = activity_prop
  )
) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::facet_grid(. ~ sex) +
  ggplot2::scale_fill_manual(
    values = activity_cols,
    name = "Activity type",
    labels = c(
      "Walking" = "Walking",
      "Moderate" = "Moderate",
      "Vigorous" = "Vigorous"
    )
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 100)
  ) +
  ggplot2::labs(
    title = "Mean activity composition within total activity quintiles",
    x = "Total activity quintile",
    y = "Mean contribution to total activity (%)"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "plain"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 12)
  ) +
  theme_figure_titles()

# -----------------------------
# Combined figure
# -----------------------------

p_fig4 <-
  p_prop_forest /
  patchwork::plot_spacer() /
  p_activity_prop_stack +
  patchwork::plot_layout(
    heights = c(1.25, 0.08, 0.9)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

p_fig4

ggplot2::ggsave(
  filename = file.path(FIGURES_DIR, "Fig_4_Activity_Composition_and_Fracture.png"),
  plot = p_fig4,
  width = 12,
  height = 8.5,
  dpi = 300
)


# =========================================================
# Figure 5 Mapping to Public Health - WHO category and Supp Table 10
# =========================================================
results_WHO_sex_df <- readRDS(file.path(DATA_DERIVED, "results_WHO_sex.rds"))
results_WHO_age_sex_df <- readRDS(file.path(DATA_DERIVED, "results_WHO_age_sex.rds"))

results_WHO_sex_table <- results_WHO_sex_df %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      grepl("Meets_guidelines", term) ~ "Meets guidelines",
      grepl("Exceeds_guidelines", term) ~ "Exceeds guidelines",
      TRUE ~ term
    ),
    OR = round(OR, 2),
    CI_lower = round(CI_lower, 2),
    CI_upper = round(CI_upper, 2),
    p_value = dplyr::case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    Global_p_value = dplyr::case_when(
      is.na(Global_p_value) ~ NA_character_,
      Global_p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", Global_p_value)
    )
  ) %>%
  dplyr::select(
    Sex = sex,
    Model = model,
    Category,
    `Odds ratio` = OR,
    `CI lower` = CI_lower,
    `CI upper` = CI_upper,
    `P value` = p_value,
    `Global P value` = Global_p_value
  )

save_table_word(
  df = results_WHO_sex_table,
  table_number = "S9",
  folder_path = TABLES_DIR,
  title = "Adjusted odds ratios for WHO physical activity categories and self-reported fracture, stratified by sex"
)




# =========================================================
# Figure 5: WHO categories by sex and age group
# Fully adjusted model only
# =========================================================

results_WHO_age_sex_df <- readRDS(
  file.path(DATA_DERIVED, "results_WHO_age_sex.rds")
)

results_WHO_age_sex_plot_df <- results_WHO_age_sex_df %>%
  dplyr::filter(model == "Fully adjusted") %>%
  dplyr::mutate(
    WHO_cat_label = sub("^WHO_MVPA_equiv_cat", "", term),
    WHO_cat_label = factor(
      WHO_cat_label,
      levels = c(
        "Meets_guidelines",
        "Exceeds_guidelines"
      ),
      labels = c(
        "Meets guidelines",
        "Exceeds guidelines"
      )
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    agegp_A0 = factor(agegp_A0, levels = c("40-49", "50-59", "60-65"))
  )

p_WHO_age_sex_fully <- ggplot2::ggplot(
  results_WHO_age_sex_plot_df,
  ggplot2::aes(
    x = WHO_cat_label,
    y = OR
  )
) +
  ggplot2::geom_point(
    size = 3,
    colour = model_cols["Fully adjusted"]
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = CI_lower,
      ymax = CI_upper
    ),
    width = 0.2,
    colour = model_cols["Fully adjusted"]
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  ggplot2::scale_y_log10() +
  ggplot2::facet_grid(
    rows = ggplot2::vars(agegp_A0),
    cols = ggplot2::vars(sex),
    labeller = ggplot2::labeller(
      agegp_A0 = function(x) paste0("Age ", x)
    )
  ) +
  ggplot2::labs(
    x = "WHO physical activity category",
    y = "Odds ratio for fracture",
    title = "Fully adjusted associations between WHO physical activity categories and fracture"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      lineheight = 0.9,
      hjust = 0.5
    ),
    strip.background = ggplot2::element_blank()
  ) +
  theme_figure_titles()

p_WHO_age_sex_fully

ggplot2::ggsave(
  filename = file.path(FIGURES_DIR, "Fig_5_alternative_WHO_PA_categories_fully_adjusted_age_sex.png"),
  plot = p_WHO_age_sex_fully,
  width = 10,
  height = 8,
  dpi = 300
)

# =========================================================
# Figure 5 Mapping to Public Health - WHO category by age group and sex
# Supplementary Table S10
# =========================================================

results_WHO_age_sex_df <- readRDS(file.path(DATA_DERIVED, "results_WHO_age_sex.rds"))

results_WHO_age_sex_table <- results_WHO_age_sex_df %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      grepl("Meets_guidelines", term) ~ "Meets guidelines",
      grepl("Exceeds_guidelines", term) ~ "Exceeds guidelines",
      TRUE ~ term
    ),
    OR = round(OR, 2),
    CI_lower = round(CI_lower, 2),
    CI_upper = round(CI_upper, 2),
    p_value = dplyr::case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    Global_p_value = dplyr::case_when(
      is.na(Global_p_value) ~ NA_character_,
      Global_p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", Global_p_value)
    )
  ) %>%
  dplyr::select(
    `Age group` = agegp_A0,
    Sex = sex,
    Model = model,
    Category,
    `Odds ratio` = OR,
    `CI lower` = CI_lower,
    `CI upper` = CI_upper,
    `P value` = p_value,
    `Global P value` = Global_p_value
  )

save_table_word(
  df = results_WHO_age_sex_table,
  table_number = "S11",
  folder_path = TABLES_DIR,
  title = "Adjusted odds ratios for WHO physical activity categories and self-reported fracture, stratified by age group and sex"
)

results_plot_age_sex_df <- results_WHO_age_sex_df %>%
  dplyr::mutate(
    WHO_cat_label = sub("^WHO_MVPA_equiv_cat", "", term),
    WHO_cat_label = factor(
      WHO_cat_label,
      levels = c(
        "Meets_guidelines",
        "Exceeds_guidelines"
      ),
      labels = c(
        "Meets guidelines",
        "Exceeds guidelines"
      )
    ),
    model = factor(
      model,
      levels = c("Unadjusted", "Minimal adjustment", "Fully adjusted")
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    agegp_A0 = factor(agegp_A0)
  )

p_WHO_age_sex <- ggplot2::ggplot(
  results_plot_age_sex_df,
  ggplot2::aes(
    x = WHO_cat_label,
    y = OR,
    colour = model
  )
) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.6),
    size = 3
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    position = ggplot2::position_dodge(width = 0.6)
  ) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(
    values = model_cols
  ) +
  ggplot2::facet_grid(agegp_A0 ~ sex) +
  ggplot2::labs(
    x = "WHO physical activity category",
    y = "Odds ratio for fracture",
    colour = "Model",
    title = "Odds ratios for WHO physical activity categories, stratified by age group and sex"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(lineheight = 0.9, hjust = 0.5),
    legend.position = "bottom"
  ) +
  theme_figure_titles()

print(p_WHO_age_sex)

ggplot2::ggsave(
  filename = file.path(FIGURES_DIR, "Fig_5_WHO_PA_categories_SRF_age_sex.png"),
  plot = p_WHO_age_sex,
  width = 10,
  height = 8,
  dpi = 300
)

# =======================================================================================
# Figure 5 and Supplementary Tables S9 and S11
# WHO physical activity categories and self-reported fracture
# =======================================================================================

# =========================================================
# Load data once
# =========================================================

IPAQ_dat2 <- readRDS(
  file.path(DATA_DERIVED, "IPAQ_dat2.rds")
)

results_WHO_sex_df <- readRDS(
  file.path(DATA_DERIVED, "results_WHO_sex.rds")
)

results_WHO_age_sex_df <- readRDS(
  file.path(DATA_DERIVED, "results_WHO_age_sex.rds")
)

# =========================================================
# Standardise age and sex labels
# =========================================================

age_levels <- c(
  "40-49",
  "50-59",
  "60-65"
)

sex_levels <- c(
  "Female",
  "Male"
)

model_levels <- c(
  "Unadjusted",
  "Minimal adjustment",
  "Fully adjusted"
)

IPAQ_dat2 <- IPAQ_dat2 %>%
  dplyr::mutate(
    agegp_A0 = dplyr::recode(
      as.character(agegp_A0),
      "60-69" = "60-65"
    ),
    agegp_A0 = factor(
      agegp_A0,
      levels = age_levels
    ),
    sex = factor(
      as.character(sex),
      levels = sex_levels
    )
  )

results_WHO_age_sex_df <- results_WHO_age_sex_df %>%
  dplyr::mutate(
    agegp_A0 = dplyr::recode(
      as.character(agegp_A0),
      "60-69" = "60-65"
    ),
    agegp_A0 = factor(
      agegp_A0,
      levels = age_levels
    ),
    sex = factor(
      as.character(sex),
      levels = sex_levels
    ),
    model = factor(
      model,
      levels = model_levels
    )
  )

results_WHO_sex_df <- results_WHO_sex_df %>%
  dplyr::mutate(
    sex = factor(
      as.character(sex),
      levels = sex_levels
    ),
    model = factor(
      model,
      levels = model_levels
    )
  )

# =======================================================================================
# Supplementary Table S9: WHO categories stratified by sex
# =======================================================================================

results_WHO_sex_table <- results_WHO_sex_df %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      grepl("Meets_guidelines", term) ~ "Meets guidelines",
      grepl("Exceeds_guidelines", term) ~ "Exceeds guidelines",
      TRUE ~ term
    ),
    OR = round(OR, 2),
    CI_lower = round(CI_lower, 2),
    CI_upper = round(CI_upper, 2),
    p_value = dplyr::case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    Global_p_value = dplyr::case_when(
      is.na(Global_p_value) ~ NA_character_,
      Global_p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", Global_p_value)
    )
  ) %>%
  dplyr::select(
    Sex = sex,
    Model = model,
    Category,
    `Odds ratio` = OR,
    `CI lower` = CI_lower,
    `CI upper` = CI_upper,
    `P value` = p_value,
    `Global P value` = Global_p_value
  )

save_table_word(
  df = results_WHO_sex_table,
  table_number = "S9",
  folder_path = TABLES_DIR,
  title = paste(
    "Adjusted odds ratios for WHO physical activity categories",
    "and self-reported fracture, stratified by sex"
  )
)

# =======================================================================================
# Supplementary Table S11: WHO categories stratified by age and sex
# =======================================================================================

results_WHO_age_sex_table <- results_WHO_age_sex_df %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      grepl("Meets_guidelines", term) ~ "Meets guidelines",
      grepl("Exceeds_guidelines", term) ~ "Exceeds guidelines",
      TRUE ~ term
    ),
    OR = round(OR, 2),
    CI_lower = round(CI_lower, 2),
    CI_upper = round(CI_upper, 2),
    p_value = dplyr::case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    Global_p_value = dplyr::case_when(
      is.na(Global_p_value) ~ NA_character_,
      Global_p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", Global_p_value)
    )
  ) %>%
  dplyr::select(
    `Age group` = agegp_A0,
    Sex = sex,
    Model = model,
    Category,
    `Odds ratio` = OR,
    `CI lower` = CI_lower,
    `CI upper` = CI_upper,
    `P value` = p_value,
    `Global P value` = Global_p_value
  )

save_table_word(
  df = results_WHO_age_sex_table,
  table_number = "S11",
  folder_path = TABLES_DIR,
  title = paste(
    "Adjusted odds ratios for WHO physical activity categories",
    "and self-reported fracture, stratified by age group and sex"
  )
)

# =======================================================================================
# Calculate N within each age-sex stratum
# =======================================================================================

# These are participants with non-missing exposure and outcome data.
# Add the model covariates to this filter if you need exact complete-case
# sample sizes for the fully adjusted models.

# =========================================================
# Calculate analysis N within each age-sex stratum
# =========================================================

WHO_strata_n <- IPAQ_dat2 %>%
  dplyr::filter(
    !is.na(sex),
    !is.na(agegp_A0),
    !is.na(WHO_MVPA_equiv_cat),
    !is.na(SRF)
  ) %>%
  dplyr::count(
    sex,
    agegp_A0,
    name = "N"
  )

WHO_strata_n

# =========================================================
# Calculate analysis N within each age-sex stratum
# =========================================================

WHO_strata_n <- IPAQ_dat2 %>%
  dplyr::filter(
    !is.na(sex),
    !is.na(agegp_A0),
    !is.na(WHO_MVPA_equiv_cat),
    !is.na(SRF)
  ) %>%
  dplyr::count(
    sex,
    agegp_A0,
    name = "N"
  )

WHO_strata_n

p_WHO_female <- plot_WHO_or_age(
  data = results_WHO_age_sex_plot_df,
  sex_value = "Female",
  title_text = "Female",
  age_labels = female_WHO_age_labels,
  y_limits = c(0.75, 1.7)
) +
  ggplot2::labs(
    tag = "A"
  )

p_WHO_male <- plot_WHO_or_age(
  data = results_WHO_age_sex_plot_df,
  sex_value = "Male",
  title_text = "Male",
  age_labels = male_WHO_age_labels,
  y_limits = c(0.75, 1.7)
) +
  ggplot2::labs(
    tag = "B"
  )

p_figure5_WHO <- patchwork::wrap_plots(
  p_WHO_female,
  p_WHO_male,
  ncol = 2,
  widths = c(1, 1),
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = paste(
      "Association between WHO physical activity categories",
      "and fracture risk by sex and age group"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 22,
        hjust = 0.5
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 22
      )
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

p_figure5_WHO

ggplot2::ggsave(
  filename = file.path(
    FIGURES_DIR,
    "Figure_5_WHO_PA_categories_SRF_by_age_sex.png"
  ),
  plot = p_figure5_WHO,
  width = 22,
  height = 8.5,
  dpi = 600,
  bg = "white"
)



