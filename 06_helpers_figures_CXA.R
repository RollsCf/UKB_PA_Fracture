# =========================================================
# Plotting helpers
# =========================================================
# =========================================================
# Colour schemes
# =========================================================

activity_cols <- c(
  "Walking"  = "#E69F00",
  "Moderate" = "#CC79A7",
  "Vigorous" = "#000000"
)

model_cols <- c(
  "Unadjusted" = "#4D4D4D",
  "Minimal adjustment" = "#56B4E9",
  "Fully adjusted" = "#0072B2"
)

sex_cols <- c(
  "Female" = "#CC79A7",
  "Male" = "#56B4E9"
)

fracture_cols <- c(
  "Other" = "#999999",
  "Wrist" = "#0072B2",
  "Ankle" = "#E69F00",
  "Arm" = "#009E73",
  "Leg" = "#CC79A7",
  "Spine" = "#F0E442",
  "Hip" = "#56B4E9"
)

# =========================================================
# Shared figure styling
# =========================================================

figure_title_size <- 15
figure_title_face <- "bold"
facet_text_size <- 11

theme_figure_titles <- function() {
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = figure_title_face,
      size = figure_title_size,
      hjust = 0.5
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = facet_text_size
    ),
    legend.title = ggplot2::element_text(face = "plain")
  )
}



plot_quintile_or <- function(
    data,
    facet_formula,
    title_text,
    facet_labeller = ggplot2::label_value
) {
  
  dodge_position <- ggplot2::position_dodge(
    width = 0.6
  )
  
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = Quintile,
      y = Odds_Ratio,
      colour = model,
      group = model
    )
  ) +
    ggplot2::geom_point(
      position = dodge_position,
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      width = 0.2,
      position = dodge_position
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "black"
    ) +
    ggplot2::facet_grid(
      facet_formula,
      labeller = facet_labeller
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(
      values = model_cols,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = "Activity quintile",
      y = "Odds ratio (95% CI)",
      colour = "Model",
      title = title_text
    ) +
    ggplot2::theme_minimal(
      base_size = 18
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "lightgrey",
        colour = "black"
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 14,
        lineheight = 1.05
      ),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    theme_figure_titles()
}

save_quintile_plot <- function(plot_obj, filename, figures_dir, width = 2400, height = 1800, res = 150) {
  png(file.path(figures_dir, filename), width = width, height = height, res = res)
  print(plot_obj)
  dev.off()
}

plot_quintile_forest <- function(data, title_text) {
  
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = OR,
      y = Quintile
    )
  ) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = CI_lower,
        xmax = CI_upper
      ),
      height = 0.2
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "black"
    ) +
    ggplot2::facet_wrap(
      ~ exposure,
      ncol = 3
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Odds ratio (log scale, 95% CI)",
      y = "Quintile",
      title = title_text
    ) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "grey90",
        colour = "black"
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_grouped_or <- function(data, x_var, title_text, facet_formula) {
  
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[x_var]],
      y = OR,
      colour = model
    )
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.6),
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.6)
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "black"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(
      values = model_cols
    ) +
    ggplot2::facet_grid(facet_formula) +
    ggplot2::labs(
      x = x_var,
      y = "Odds ratio (log scale, 95% CI)",
      colour = "Model",
      title = title_text
    ) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "lightgrey",
        colour = "black"
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# =========================================================
# Additional helper functions for sensitivity analyses
# =========================================================

plot_mutual_quintile_or <- function(data, title_text, y_limits = c(0.8, 2.0)) {
  
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = Quintile,
      y = OR,
      colour = exposure,
      group = exposure
    )
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.6),
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.6)
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "black"
    ) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = seq(
        from = y_limits[1],
        to = y_limits[2],
        by = 0.2
      )
    ) +
    ggplot2::scale_colour_manual(
      name = "Mutually adjusted activity",
      values = activity_cols,
      labels = c(
        "Walking" = "Walking",
        "Moderate" = "Moderate",
        "Vigorous" = "Vigorous"
      )
    ) +
    ggplot2::labs(
      x = "Activity quintile",
      y = "Odds ratio (95% CI)",
      title = title_text
    ) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "plain")
    )
}

plot_prop_within_q <- function(data, title_text, y_limits = c(0.85, 1.2)) {
  
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = MET_total_q,
      y = OR,
      colour = activity_prop,
      group = activity_prop
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "black"
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.6),
      size = 2.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      linewidth = 0.7
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::scale_colour_manual(
      name = "Activity type",
      values = activity_cols,
      labels = c(
        "Walking" = "Walking",
        "Moderate" = "Moderate",
        "Vigorous" = "Vigorous"
      )
    ) +
    ggplot2::labs(
      x = "Total activity quintile",
      y = "OR per 10% higher activity proportion",
      title = title_text
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    theme_figure_titles()
}

plot_quintile_forest_sex_combined <- function(data, title_text = NULL) {
  
  dodge <- ggplot2::position_dodge(width = 0.55)
  
  data %>%
    dplyr::mutate(
      sex = factor(sex, levels = c("Female", "Male")),
      Quintile = factor(Quintile, levels = paste0("Q", 1:5))
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = OR,
        y = Quintile,
        colour = sex
      )
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "grey50"
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = CI_lower, xmax = CI_upper),
      height = 0.18,
      position = dodge,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      position = dodge,
      size = 2.4
    ) +
    ggplot2::scale_colour_manual(
      values = sex_cols
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_grid(
      outcome ~ exposure,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Odds ratio (log scale, 95% CI)",
      y = "Quintile",
      colour = "Sex"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      strip.text.x = ggplot2::element_text(size = 9),
      strip.text.y = ggplot2::element_text(size = 9)
    )
}
  
  # =========================================================
  # OR plot by sex and age group
  # =========================================================
  
plot_total_or_age <- function(
    data,
    sex_value,
    title_text,
    age_labels,
    y_limits = c(0.75, 2.1)
) {
  
  plot_quintile_or(
    data = data %>%
      dplyr::filter(sex == sex_value),
    facet_formula = . ~ agegp_A0,
    title_text = title_text,
    facet_labeller = ggplot2::labeller(
      agegp_A0 = ggplot2::as_labeller(
        age_labels
      )
    )
  ) +
    ggplot2::coord_cartesian(
      ylim = y_limits
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Odds ratio for fracture"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      strip.text.y = ggplot2::element_blank(),
      strip.background.y = ggplot2::element_blank()
    )
}

  # =========================================================
  # Helper 2: fracture-site stacked bar by sex and age group
  # =========================================================
  
plot_fracture_sites_age_q <- function(data, sex_value, show_legend = TRUE) {
  
  ggplot2::ggplot(
    data %>% dplyr::filter(sex == sex_value),
    ggplot2::aes(
      x = MET_total_q,
      y = prop,
      fill = fracture_type
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(. ~ agegp_A0) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format()
    ) +
    ggplot2::scale_fill_manual(
      values = fracture_cols
    ) +
    ggplot2::labs(
      x = "Activity quintile",
      y = "Reported fracture sites (%)",
      fill = "Fracture site"
    ) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "grey90",
        colour = "black"
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = ifelse(show_legend, "right", "none")
    )
}
  
### Figure 4 helper
prep_mutual_age_sex_plot_df <- function(df) {
  
  df <- df %>%
    dplyr::mutate(
      Quintile = sub(".*(Q[1-5])$", "\\1", term),
      exposure = dplyr::case_when(
        exposure == "MET_vig"  ~ "Vigorous",
        exposure == "MET_mod"  ~ "Moderate",
        exposure == "MET_walk" ~ "Walking",
        TRUE ~ exposure
      ),
      sex = factor(sex, levels = c("Female", "Male")),
      agegp_A0 = factor(agegp_A0, levels = c("40-49", "50-59", "60-65")),
      exposure = factor(exposure, levels = c("Walking", "Moderate", "Vigorous"))
    ) %>%
    dplyr::filter(Quintile %in% paste0("Q", 2:5))
  
  ref_rows <- df %>%
    dplyr::distinct(sex, agegp_A0, outcome, exposure, model) %>%
    dplyr::mutate(
      term = NA_character_,
      Quintile = "Q1",
      OR = 1,
      CI_lower = 1,
      CI_upper = 1,
      p_value = NA_real_
    )
  
  df %>%
    dplyr::bind_rows(ref_rows) %>%
    dplyr::mutate(
      Quintile = factor(Quintile, levels = paste0("Q", 1:5))
    )
}


prep_absolute_activity_profile_df <- function(data) {
  
  data %>%
    dplyr::group_by(agegp_A0, sex, MET_total_q) %>%
    dplyr::summarise(
      Walking  = stats::median(mins_wk_walk, na.rm = TRUE),
      Moderate = stats::median(mins_wk_mod, na.rm = TRUE),
      Vigorous = stats::median(mins_wk_vig, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(Walking, Moderate, Vigorous),
      names_to = "intensity",
      values_to = "minutes"
    ) %>%
    dplyr::mutate(
      intensity = factor(intensity, levels = c("Walking", "Moderate", "Vigorous")),
      sex = factor(sex, levels = c("Female", "Male")),
      agegp_A0 = factor(agegp_A0, levels = c("40-49", "50-59", "60-65")),
      MET_total_q = factor(MET_total_q, levels = paste0("Q", 1:5))
    )
}

#### dxa forests

make_dxa_forest_plot <- function(tidy_no_bmd, tidy_bmd, plot_title) {
  
  plot_dat <- dplyr::bind_rows(
    tidy_no_bmd %>%
      dplyr::filter(model == "Fully adjusted") %>%
      dplyr::mutate(model_set = "DXA cohort"),
    
    tidy_bmd %>%
      dplyr::filter(model == "Fully adjusted + FN BMD") %>%
      dplyr::mutate(model_set = "DXA cohort + FN BMD")
  ) %>%
    dplyr::mutate(
      quintile = factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5")),
      model_set = factor(
        model_set,
        levels = c("DXA cohort", "DXA cohort + FN BMD")
      )
    )
  
  ggplot(
    plot_dat,
    aes(
      x = quintile,
      y = estimate,
      ymin = conf.low,
      ymax = conf.high,
      colour = model_set,
      group = model_set
    )
  ) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_pointrange(
      position = position_dodge(width = 0.6),
      linewidth = 0.7
    ) +
    labs(
      title = plot_title,
      x = "Total PA quintile",
      y = "Hazard ratio (95% CI)",
      colour = "Model"
    ) +
    forest_panel_theme
}


make_age_strip_labels <- function(n_data, sex_value) {
  
  labels_df <- n_data %>%
    dplyr::filter(sex == sex_value) %>%
    dplyr::mutate(
      strip_label = paste0(
        agegp_A0,
        "\nN = ",
        format(
          N,
          big.mark = ",",
          scientific = FALSE
        )
      )
    )
  
  stats::setNames(
    labels_df$strip_label,
    as.character(labels_df$agegp_A0)
  )
}

female_age_labels <- make_age_strip_labels(
  age_sex_n,
  "Female"
)

male_age_labels <- make_age_strip_labels(
  age_sex_n,
  "Male"
)

plot_mutual_or_age <- function(
    data,
    sex_value,
    title_text,
    age_labels,
    y_limits = c(0.8, 2.0)
) {
  
  dodge_position <- ggplot2::position_dodge(
    width = 0.6
  )
  
  ggplot2::ggplot(
    data %>%
      dplyr::filter(sex == sex_value),
    ggplot2::aes(
      x = Quintile,
      y = OR,
      colour = exposure,
      group = exposure
    )
  ) +
    ggplot2::geom_point(
      position = dodge_position,
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      position = dodge_position,
      width = 0.2
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "black"
    ) +
    ggplot2::facet_grid(
      . ~ agegp_A0,
      labeller = ggplot2::labeller(
        agegp_A0 = ggplot2::as_labeller(
          age_labels
        )
      )
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::coord_cartesian(
      ylim = y_limits
    ) +
    ggplot2::scale_colour_manual(
      values = activity_cols,
      breaks = c(
        "Walking",
        "Moderate",
        "Vigorous"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = "Odds ratio for fracture",
      colour = "Mutually adjusted activity"
    ) +
    ggplot2::theme_minimal(
      base_size = 18
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 18,
        hjust = 0.5
      ),
      strip.background = ggplot2::element_rect(
        fill = "lightgrey",
        colour = "black"
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 14,
        lineheight = 1.05
      ),
      axis.title = ggplot2::element_text(
        size = 18
      ),
      axis.text = ggplot2::element_text(
        size = 16
      ),
      legend.title = ggplot2::element_text(
        size = 16
      ),
      legend.text = ggplot2::element_text(
        size = 16
      ),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    theme_figure_titles()
}


make_WHO_age_strip_labels <- function(n_data, sex_value) {
  
  labels_df <- n_data %>%
    dplyr::filter(
      sex == sex_value
    ) %>%
    dplyr::mutate(
      strip_label = paste0(
        agegp_A0,
        "\nN = ",
        format(
          N,
          big.mark = ",",
          scientific = FALSE
        )
      )
    )
  
  stats::setNames(
    labels_df$strip_label,
    as.character(labels_df$agegp_A0)
  )
}

female_WHO_age_labels <- make_WHO_age_strip_labels(
  WHO_strata_n,
  "Female"
)

male_WHO_age_labels <- make_WHO_age_strip_labels(
  WHO_strata_n,
  "Male"
)

plot_WHO_or_age <- function(
    data,
    sex_value,
    title_text,
    age_labels,
    y_limits = c(0.75, 1.7)
) {
  
  dodge_position <- ggplot2::position_dodge(
    width = 0.65
  )
  
  ggplot2::ggplot(
    data %>%
      dplyr::filter(
        sex == sex_value
      ),
    ggplot2::aes(
      x = WHO_cat_label,
      y = OR,
      colour = model,
      group = model
    )
  ) +
    ggplot2::geom_point(
      position = dodge_position,
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      width = 0.12,
      position = dodge_position
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.6
    ) +
    ggplot2::facet_grid(
      . ~ agegp_A0,
      drop = FALSE,
      labeller = ggplot2::labeller(
        agegp_A0 = ggplot2::as_labeller(
          age_labels
        )
      )
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::coord_cartesian(
      ylim = y_limits
    ) +
    ggplot2::scale_colour_manual(
      values = model_cols,
      breaks = model_levels,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = "Odds ratio for fracture",
      colour = "Model"
    ) +
    ggplot2::theme_minimal(
      base_size = 18
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 18,
        hjust = 0.5
      ),
      strip.background = ggplot2::element_rect(
        fill = "lightgrey",
        colour = "black"
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 14,
        lineheight = 1.05
      ),
      axis.title = ggplot2::element_text(
        size = 18
      ),
      axis.text.y = ggplot2::element_text(
        size = 16
      ),
      axis.text.x = ggplot2::element_text(
        size = 14,
        lineheight = 0.9,
        hjust = 0.5
      ),
      legend.title = ggplot2::element_text(
        size = 16
      ),
      legend.text = ggplot2::element_text(
        size = 16
      ),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(
        0.8,
        "cm"
      )
    ) +
    theme_figure_titles()
}

