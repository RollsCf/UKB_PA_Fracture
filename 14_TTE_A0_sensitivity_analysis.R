# =========================================================
# TIME-TO-EVENT SENSITIVITY DATASET:
# Alternative exposure = total number of activity days
# Larger cohort because MET_total is not required
# =========================================================

source(here::here("scripts/00_setup.R"))
source(here("scripts/03_helpers_general.R"))
source(here("scripts/04_helpers_CXA.R"))
source(here("scripts/05_helpers_survival.R"))

TTE_sensitivity_dat <- readRDS(file.path(DATA_DERIVED, "TTE_analysis_dat.Rds"))

# ---------------------------------------------------------
# Keep variable names consistent with TTE models
# Do NOT strip _raw or _clean if your Cox functions use them
# ---------------------------------------------------------

# Create middle aged cohort A0 for all subsequent analysis of self-reported activity

TTE_sensitivity_dat <- TTE_sensitivity_dat %>%
  dplyr::filter(dplyr::between(age_A0_raw, 40, 65))

# Ensure sex is a factor
TTE_sensitivity_dat$sex <- factor(TTE_sensitivity_dat$sex)

TTE_sensitivity_dat <- TTE_sensitivity_dat %>%
  mutate(
    fragility = factor(event_fragility_PA,   levels = c(0, 1), labels = c("No", "Yes")),
  )

check_eid(TTE_sensitivity_dat)

# ---------------------------------------------------------
# Create alternative activity exposure
# ---------------------------------------------------------

TTE_sensitivity_dat <- TTE_sensitivity_dat %>%
  dplyr::mutate(
    total_num_days = dplyr::if_else(
      !is.na(num_day_walk_clean) &
        !is.na(num_days_mod_clean) &
        !is.na(num_days_vig_clean),
      num_day_walk_clean + num_days_mod_clean + num_days_vig_clean,
      NA_real_
    ),
    log_total_num_days = log(total_num_days + 1),
    total_num_days_quintile = dplyr::ntile(total_num_days, 5),
    total_num_days_quintile_f = factor(
      total_num_days_quintile,
      levels = 1:5,
      labels = c("Q1", "Q2", "Q3", "Q4", "Q5")
    ),
    sex_raw = factor(sex_raw),
    ethnicity_derived = factor(ethnicity_derived),
    education_level = factor(education_level)
  )

# =========================================================
# EXCLUSION TABLE
# =========================================================

tab_exc <- data.frame(
  Exclusion = "Starting cohort age 40-65",
  Number_excluded = NA,
  Number_remaining = nrow(TTE_sensitivity_dat)
)

add_exclusion <- function(data, condition, label, tab) {
  nb <- nrow(data)
  data <- data[condition, ]
  
  tab <- rbind(
    tab,
    data.frame(
      Exclusion = label,
      Number_excluded = nb - nrow(data),
      Number_remaining = nrow(data)
    )
  )
  
  list(data = data, tab = tab)
}

# ---------------------------------------------------------
# HES outcome variables required for time-to-event analysis
# ---------------------------------------------------------

tmp <- add_exclusion(
  TTE_sensitivity_dat,
  !is.na(TTE_sensitivity_dat$age_entry_PA) &
    !is.na(TTE_sensitivity_dat$age_exit_fragility_PA) &
    !is.na(TTE_sensitivity_dat$event_fragility_PA),
  "Missing HES time-to-event outcome data",
  tab_exc
)

TTE_sensitivity_dat <- tmp$data
tab_exc <- tmp$tab

# ---------------------------------------------------------
# Covariate exclusions
# ---------------------------------------------------------

tmp <- add_exclusion(
  TTE_sensitivity_dat,
  !is.na(TTE_sensitivity_dat$ethnicity_derived),
  "Missing ethnicity data",
  tab_exc
)

TTE_sensitivity_dat <- tmp$data
tab_exc <- tmp$tab

tmp <- add_exclusion(
  TTE_sensitivity_dat,
  !is.na(TTE_sensitivity_dat$tdi_raw),
  "Missing deprivation data",
  tab_exc
)

TTE_sensitivity_dat <- tmp$data
tab_exc <- tmp$tab

tmp <- add_exclusion(
  TTE_sensitivity_dat,
  !is.na(TTE_sensitivity_dat$education_level),
  "Missing education data",
  tab_exc
)

TTE_sensitivity_dat <- tmp$data
tab_exc <- tmp$tab

tmp <- add_exclusion(
  TTE_sensitivity_dat,
  !is.na(TTE_sensitivity_dat$height_clean) &
    !is.na(TTE_sensitivity_dat$weight_clean),
  "Missing height or weight data",
  tab_exc
)

TTE_sensitivity_dat <- tmp$data
tab_exc <- tmp$tab

# ---------------------------------------------------------
# Alternative exposure exclusion
# ---------------------------------------------------------

tmp <- add_exclusion(
  TTE_sensitivity_dat,
  !is.na(TTE_sensitivity_dat$total_num_days),
  "Missing total activity days data",
  tab_exc
)

TTE_sensitivity_dat <- tmp$data
tab_exc <- tmp$tab

tab_exc

# =========================================================
# SENSITIVITY ANALYSIS:
# total_num_days as alternative exposure
# Sex-stratified quintile models to mirror main analysis
# =========================================================

# ---------------------------------------------------------
# 1. Check exposure distribution
# ---------------------------------------------------------

summary(TTE_sensitivity_dat$total_num_days)

quantile(
  TTE_sensitivity_dat$total_num_days,
  probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
  na.rm = TRUE
)

ggplot2::ggplot(TTE_sensitivity_dat, ggplot2::aes(x = total_num_days)) +
  ggplot2::geom_histogram(binwidth = 1, fill = "steelblue", colour = "white") +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Distribution of total_num_days",
    x = "Total activity days",
    y = "Count"
  )

# =========================================================
# Sensitivity: total_num_days quintiles
# Sex-stratified, fully adjusted only
# =========================================================


# Create quintiles
TTE_sensitivity_dat <- TTE_sensitivity_dat %>%
  dplyr::mutate(
    total_num_days_quintile = dplyr::ntile(total_num_days, 5),
    total_num_days_quintile_f = factor(
      total_num_days_quintile,
      levels = 1:5,
      labels = paste0("Q", 1:5)
    ),
    total_num_days_quintile_f = stats::relevel(total_num_days_quintile_f, ref = "Q1")
  )

# Sex-specific datasets
TTE_sensitivity_female <- TTE_sensitivity_dat %>%
  dplyr::filter(sex_raw == "Female")

TTE_sensitivity_male <- TTE_sensitivity_dat %>%
  dplyr::filter(sex_raw == "Male")

# Fully adjusted quintile formula
f_total_days_quintile <- survival::Surv(
  age_entry_PA,
  age_exit_fragility_PA,
  event_fragility_PA
) ~ total_num_days_quintile_f +
  ethnicity_derived +
  height_clean +
  weight_clean +
  tdi_raw +
  education_level

# Fit models
m_total_days_quintile_female <- survival::coxph(
  f_total_days_quintile,
  data = TTE_sensitivity_female
)

m_total_days_quintile_male <- survival::coxph(
  f_total_days_quintile,
  data = TTE_sensitivity_male
)

# Extract HRs
extract_total_days_quintiles <- function(model, sex_label) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(grepl("^total_num_days_quintile_f", term)) %>%
    dplyr::mutate(
      Sex = sex_label,
      Exposure = "Total activity days",
      Quintile = sub("^total_num_days_quintile_f", "", term)
    ) %>%
    dplyr::select(Sex, Exposure, Quintile, estimate, conf.low, conf.high, p.value)
}

total_days_quintile_results <- dplyr::bind_rows(
  extract_total_days_quintiles(m_total_days_quintile_female, "Female"),
  extract_total_days_quintiles(m_total_days_quintile_male, "Male")
)

# HR table
total_days_hr_table <- total_days_quintile_results %>%
  dplyr::mutate(
    HR_95CI = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high),
    p_value = formatC(p.value, format = "f", digits = 3)
  ) %>%
  dplyr::select(Sex, Quintile, HR_95CI, p_value)

total_days_hr_table

write.csv(
  total_days_hr_table,
  file.path(TTE_TABLES_DIR, "Table_sensitivity_total_num_days_quintiles.csv"),
  row.names = FALSE
)

# =========================================================
# Plot total_num_days beside original MET quintile plots
# Requires:
# p_fragility_quint_female
# p_fragility_quint_male
# =========================================================

make_total_days_plot <- function(results_df, sex_label) {
  
  plot_df <- results_df %>%
    dplyr::filter(Sex == sex_label) %>%
    dplyr::rename(
      HR = estimate,
      HR_low = conf.low,
      HR_high = conf.high
    ) %>%
    dplyr::select(Quintile, HR, HR_low, HR_high) %>%
    dplyr::bind_rows(
      tibble::tibble(
        Quintile = "Q1",
        HR = 1,
        HR_low = 1,
        HR_high = 1
      ),
      .
    ) %>%
    dplyr::mutate(
      Quintile = factor(Quintile, levels = paste0("Q", 1:5))
    )
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = Quintile, y = HR)) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = HR_low, ymax = HR_high),
      width = 0.15
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Total activity days quintile",
      y = "Hazard ratio vs Q1",
      title = paste0(sex_label, ": total activity days")
    ) +
    ggplot2::theme_minimal()
}

p_total_days_female <- make_total_days_plot(total_days_quintile_results, "Female")
p_total_days_male <- make_total_days_plot(total_days_quintile_results, "Male")

p_compare_all <- 
  (p_fragility_quint_female + p_total_days_female) /
  (p_fragility_quint_male + p_total_days_male) +
  patchwork::plot_annotation(
    title = "Comparison of MET total and total activity days quintile models"
  )

print(p_compare_all)

ggplot2::ggsave(
  filename = file.path(TTE_FIGURES_DIR, "Fig_MET_vs_total_days_quintiles_side_by_side.png"),
  plot = p_compare_all,
  width = 12,
  height = 8,
  dpi = 300
)
