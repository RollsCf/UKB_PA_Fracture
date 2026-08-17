# This code describes the derivation of variables for use in the final analysis - it uses variables created in 01_data_cleaning.R
# Variables that do not need derivation are renamed for ease in analysis ie BMI_raw to BMI

source(here::here("scripts/00_setup.R"))

# load the files for Assessment centre 0, HES and death, 

dat_clean   <- readRDS(file.path(DATA_DERIVED, "baseline_clean.Rds"))

#### Derivations for cross-sectional analysis ##################


# Print column names
names(dat_clean)

# If no derivation needed then drop the word _raw for the analysis to keep coding tight

dat_clean$ethnicity_derived <- case_when(
  dat_clean$ethnicity_clean %in% c(
    "White", "British", "Irish", "Any other white background"
  ) ~ "White",
  
  dat_clean$ethnicity_clean %in% c(
    "Mixed",
    "White and Black Caribbean",
    "White and Black African",
    "White and Asian",
    "Any other mixed background"
  ) ~ "Mixed",
  
  dat_clean$ethnicity_clean %in% c(
    "Asian or Asian British",
    "Indian",
    "Pakistani",
    "Bangladeshi",
    "Any other Asian background"
  ) ~ "Asian",
  
  dat_clean$ethnicity_clean %in% c(
    "Black or Black British",
    "Caribbean",
    "African",
    "Any other Black background",
    "Chinese",
    "Other ethnic group"
  ) ~ "Black or other",
  
  TRUE ~ NA_character_
)

## Check
table(dat_clean$ethnicity_derived, useNA = "ifany")


# UKB assessment centre
# Field 54	# Assessment Centre location (A0)
# To allow me to look at geographic variation I re-coded Assessment centre by Country it was located in.
# Codes were as follows

dat_clean <- dat_clean %>%
  dplyr::mutate(
    assessment_country = dplyr::case_when(
      ukb_assess_center_raw %in% c(
        "Barts",
        "Birmingham",
        "Bristol",
        "Bury",
        "Cheadle (revisit)",
        "Croydon",
        "Hounslow",
        "Leeds",
        "Liverpool",
        "Manchester",
        "Middlesbrough",
        "Newcastle",
        "Nottingham",
        "Oxford",
        "Reading",
        "Sheffield",
        "Stockport (pilot)",
        "Stoke",
        "Cheadle (imaging)",
        "Reading (imaging)",
        "Newcastle (imaging)",
        "Bristol (imaging)"
      ) ~ "England",
      
      ukb_assess_center_raw %in% c(
        "Edinburgh",
        "Glasgow"
      ) ~ "Scotland",
      
      ukb_assess_center_raw %in% c(
        "Cardiff",
        "Swansea",
        "Wrexham"
      ) ~ "Wales",
      
      TRUE ~ NA_character_
    )
  )


# create agegroup_A0 based on age_A0_raw

middle_aged_A0_dat <- dat_clean %>%
  filter(between(age_A0_raw, 40, 65))

# Create age at A0 groups
dat_clean$agegp_A0 <-
  cut(
    dat_clean$age_A0_raw,
    breaks = c(40, 50, 60, 70, 80),
    right = FALSE,
    labels = c("40-49", "50-59", "60-69", "70-79")
  )

# Create binary Over or under age 50 at baseline
dat_clean$age_over_50 <- factor(
  ifelse(dat_clean$age_A0_raw > 50, "Over50", "50_or_less"),
  levels = c("50_or_less", "Over50")
)

# UK Biobank doesn't contain day of birth as it would be unnecessary identifying information.
# so we roughly impute it as the 15th of the birth month.

# Add date of birth
dat_clean$approx_dob_derived <-
  as.Date(paste(dat_clean$year_birth_raw, dat_clean$month_birth_raw, "15", sep = "-"),
          "%Y-%B-%d") 

# derive new education levels
# define education level ranking
edu_levels <- c(
  "School education or less" = 1,
  "College or vocational training" = 2,
  "University degree or above" = 3
)

# derive new education levels
derive_education_level <- function(qual_vec) {
  sapply(qual_vec, function(x) {
    
    if (is.na(x)) return(NA_character_)
    
    # Split multiple entries by '|', trim whitespace, lowercase
    quals <- trimws(tolower(unlist(strsplit(x, "\\|"))))
    
    # Map each qualification to a category
    mapped <- sapply(quals, function(q) {
      if (q %in% tolower("College or University degree")) {
        "University degree or above"
      } else if (q %in% tolower(c(
        "NVQ or HND or HNC or equivalent",
        "Other professional qualifications eg: nursing, teaching"
      ))) {
        "College or vocational training"
      } else if (q %in% tolower(c(
        "A levels/AS levels or equivalent",
        "O levels/GCSEs or equivalent",
        "CSEs or equivalent",
        "None of the above"
      ))) {
        "School education or less"
      } else {
        NA_character_
      }
    })
    
    # Pick the highest-ranked category
    if (all(is.na(mapped))) {
      NA_character_
    } else {
      names(edu_levels)[which.max(edu_levels[mapped])]
    }
  })
}

# apply function
dat_clean$education_level <- derive_education_level(dat_clean$qualification_raw)

# check
table(dat_clean$education_level, useNA = "ifany")

# convert to ordered factor
dat_clean$education_level <- factor(
  dat_clean$education_level,
  levels = names(edu_levels),
  ordered = TRUE
)


# Anthropometrics

# # Derived variables based on WHO guidance
# WHO recommends WC as a measure of abdominal (central) obesity, because it strongly correlates with visceral fat and cardiometabolic risk.

# Cut-offs for Europeans (typical WHO guidance):
# Sex	  Normal	  Increased risk	    Substantially increased risk
# Men	    <94 cm	  94–102 cm	          >102 cm
# Women	  <80 cm	  80–88 cm	          >88 cm


dat_clean <- dat_clean %>%
  mutate(
    central_obesity_derived = case_when(
      sex_raw == "Male" & waist_circ_clean >= 102 ~ "Substantially increased risk",
      sex_raw == "Male" & waist_circ_clean >= 94  ~ "Increased risk",
      sex_raw == "Male" & waist_circ_clean < 94   ~ "Normal",
      sex_raw == "Female" & waist_circ_clean >= 88 ~ "Substantially increased risk",
      sex_raw == "Female" & waist_circ_clean >= 80 ~ "Increased risk",
      sex_raw == "Female" & waist_circ_clean < 80  ~ "Normal",
      TRUE ~ NA_character_
    ),
    BMI_category_derived = case_when(
      BMI_clean < 25                     ~ "Normal",
      BMI_clean >= 25 & BMI_clean < 30  ~ "Overweight",
      BMI_clean >= 30                    ~ "Obese",
      TRUE ~ NA_character_
    )
  )

#### Physical activity measures ####


#### Helper functions ####

# Truncate duration at 180 min/day (IPAQ rule)
truncate_duration <- function(duration, max_per_day = 180) {
  case_when(
    is.na(duration) ~ NA_real_,
    TRUE ~ pmin(duration, max_per_day)
  )
}

# Calculate weekly minutes
calc_mins_wk <- function(days, duration) {
  case_when(
    is.na(days) ~ NA_real_,
    days == 0 ~ 0,
    is.na(duration) ~ NA_real_,
    TRUE ~ days * duration
  )
}

# Calculate MET-minutes
calc_MET <- function(mins, factor) {
  if_else(!is.na(mins), mins * factor, NA_real_)
}

# Strict complete-case sum
strict_sum <- function(...) {
  vals <- cbind(...)
  if_else(rowSums(is.na(vals)) == 0, rowSums(vals), NA_real_)
}

# ---- Main pipeline ----

dat_clean <- dat_clean %>%
  mutate(
    # =========================
    # 1. Truncate duration (per day)
    # =========================
    duration_mod_trunc  = truncate_duration(duration_mod_clean),
    duration_vig_trunc  = truncate_duration(duration_vig_clean),
    duration_walk_trunc = truncate_duration(duration_walk_clean),
    
    # =========================
    # 2. Weekly minutes (original)
    # =========================
    cc_mins_wk_mod  = calc_mins_wk(num_days_mod_clean, duration_mod_clean),
    cc_mins_wk_vig  = calc_mins_wk(num_days_vig_clean, duration_vig_clean),
    cc_mins_wk_walk = calc_mins_wk(num_day_walk_clean, duration_walk_clean),
    
    # =========================
    # 3. Weekly minutes (truncated)
    # =========================
    cc_mins_wk_mod_trunc  = calc_mins_wk(num_days_mod_clean, duration_mod_trunc),
    cc_mins_wk_vig_trunc  = calc_mins_wk(num_days_vig_clean, duration_vig_trunc),
    cc_mins_wk_walk_trunc = calc_mins_wk(num_day_walk_clean, duration_walk_trunc),
    
    # =========================
    # 4. MET-minutes (original)
    # =========================
    cc_MET_mod  = calc_MET(cc_mins_wk_mod, 4),
    cc_MET_vig  = calc_MET(cc_mins_wk_vig, 8),
    cc_MET_walk = calc_MET(cc_mins_wk_walk, 3.3),
    
    # =========================
    # 5. MET-minutes (truncated)
    # =========================
    cc_MET_mod_trunc  = calc_MET(cc_mins_wk_mod_trunc, 4),
    cc_MET_vig_trunc  = calc_MET(cc_mins_wk_vig_trunc, 8),
    cc_MET_walk_trunc = calc_MET(cc_mins_wk_walk_trunc, 3.3),
    
    # =========================
    # 6. MVPA (moderate + vigorous only)
    # =========================
    cc_mins_wk_MVPA        = strict_sum(cc_mins_wk_mod, cc_mins_wk_vig),
    cc_MET_MVPA            = strict_sum(cc_MET_mod, cc_MET_vig),
    
    cc_mins_wk_MVPA_trunc  = strict_sum(cc_mins_wk_mod_trunc, cc_mins_wk_vig_trunc),
    cc_MET_MVPA_trunc      = strict_sum(cc_MET_mod_trunc, cc_MET_vig_trunc),
    
    # =========================
    # 7. Total activity (includes walking)
    # =========================
    cc_mins_wk_total       = strict_sum(cc_mins_wk_mod, cc_mins_wk_vig, cc_mins_wk_walk),
    cc_MET_total           = strict_sum(cc_MET_mod, cc_MET_vig, cc_MET_walk),
    
    cc_mins_wk_total_trunc = strict_sum(cc_mins_wk_mod_trunc, cc_mins_wk_vig_trunc, cc_mins_wk_walk_trunc),
    cc_MET_total_trunc     = strict_sum(cc_MET_mod_trunc, cc_MET_vig_trunc, cc_MET_walk_trunc),
    
    # =========================
    # 8. Binary indicators (clean + meaningful)
    # =========================
    any_mod        = as.integer(cc_mins_wk_mod > 0),
    any_vig        = as.integer(cc_mins_wk_vig > 0),
    any_walk       = as.integer(cc_mins_wk_walk > 0),
    
    any_MVPA       = as.integer(cc_mins_wk_MVPA > 0),
    any_MVPA_trunc = as.integer(cc_mins_wk_MVPA_trunc > 0)
  )



# =========================================================
# CREATE WHO MODERATE-EQUIVALENT PHYSICAL ACTIVITY VARIABLES
# =========================================================
# creates a continuous mins/week measure (cont), and a catagorical based on meeting thresholds

dat_clean <- dat_clean %>%
  dplyr::mutate(
    WHO_MVPA_equiv_cont =
      as.numeric(as.character(cc_mins_wk_walk_trunc)) +
      as.numeric(as.character(cc_mins_wk_mod_trunc)) +
      2 * as.numeric(as.character(cc_mins_wk_vig_trunc))
  ) %>%
  dplyr::mutate(
    WHO_MVPA_equiv_cat = cut(
      WHO_MVPA_equiv_cont,
      breaks = c(-Inf, 150, 300, Inf),
      labels = c(
        "Below_guidelines",
        "Meets_guidelines",
        "Exceeds_guidelines"
      ),
      include.lowest = TRUE,
      right = TRUE
    ),
    WHO_MVPA_equiv_cat = factor(
      WHO_MVPA_equiv_cat,
      levels = c(
        "Below_guidelines",
        "Meets_guidelines",
        "Exceeds_guidelines"
      )
    )
  )

###### Outcome data variables

# Create binary fracture indicator for a given bone

# Define fracture sites
fracture_sites <- c("Ankle", "Leg", "Hip", "Spine", "Wrist", "Arm", "Other bones")

# Helper function
fracture_indicator <- function(fracture_list, bone_name, self_reported) {
  case_when(
    map_lgl(fracture_list, ~ bone_name %in% .x) ~ 1L,
    self_reported %in% c("Yes", "No") ~ 0L,
    TRUE ~ NA_integer_
  )
}

dat_clean <- dat_clean %>%
  dplyr::mutate(fracture_list = str_split(fracture_location_A0_clean, "\\|")) %>%
  # Dynamically create columns for all bones
  { 
    tmp <- .
    for(bone in fracture_sites) {
      tmp[[bone]] <- fracture_indicator(tmp$fracture_list, bone, tmp$self_reported_fracture_A0_clean)
    }
    tmp
  } %>%
  dplyr::select(-fracture_list) 

dat_clean %>%
  dplyr::summarise(
    dplyr::across(
      c(Ankle, Leg, Hip, Spine, Wrist, Arm, `Other bones`),
      ~ sum(.x, na.rm = TRUE)
    )
  )


#### Additional descriptive variables for cross-sectional analysis exploration

# Create binary variables for use in prevalence/radar plot

# Smoking variable derived as ever or never
dat_clean <- dat_clean %>%
  mutate(
    smoking_binary = case_when(
      smoking_clean == "Never" ~ "Never",
      is.na(smoking_clean) ~ NA_character_,
      TRUE ~ "Ever"
    )
  )


# risk taking

dat_clean <- dat_clean %>%
  mutate(
    risk_binary = case_when(
      risk_taking_clean == "No" ~ "No",
      is.na(risk_taking_clean) ~ NA_character_,
      TRUE ~ "Yes"
    )
  )

# drive fast

dat_clean <- dat_clean %>%
  mutate(
    drive_fast_clean = ifelse(drive_fast_clean == "do not drive on the motorway", NA, drive_fast_clean),
    drive_fast_binary = case_when(
      is.na(drive_fast_clean) ~ NA_character_,
      drive_fast_clean == "Never/rarely" ~ "No",
      TRUE ~ "Yes"
    )
  )


# heel bmd
# Not DXA so can't use T-scores
# Assume bottom 20% of values are low BMD
# Remove 0 values which most likely represent measurement error

dat_clean <- dat_clean %>%
  mutate(
    # Step 1: remove invalid values
    heel_BMD_clean = ifelse(heel_BMD_clean == 0, NA, heel_BMD_clean)
  )

# Step 2: compute 20th percentile
bmd_cut <- quantile(dat_clean$heel_BMD_clean, probs = 0.20, na.rm = TRUE)

# Step 3: create binary variable
dat_clean <- dat_clean %>%
  mutate(
    low_BMD = ifelse(heel_BMD_clean <= bmd_cut, 1, 0)
  )


# Low BMI
dat_clean <- dat_clean %>%
  mutate(
    bmi_reds_binary = case_when(
      is.na(BMI_clean) ~ NA_character_,
      BMI_clean < 18.5 ~ "Yes",
      BMI_clean >= 18.5 ~ "No"
    ),
    bmi_reds_binary = factor(bmi_reds_binary, levels = c("No", "Yes"))
  )

# Falls last year
dat_clean <-dat_clean %>%
  mutate(
    falls_binary = case_when(
      is.na(self_reported_falls_clean) ~ NA_character_,
      self_reported_falls_clean == "No falls" ~ "No",
      TRUE ~"Yes"
    )
  )

# Create WHO catagories based on MET mins MVPA

# -------------------------------
# 1. WHO cut-points
# -------------------------------

met_cuts <- c(-Inf, 600, 1200, Inf)
met_labels <- c("Below_guidelines",
                "Meets_guidelines",
                "Exceeds_guidelines")

# -------------------------------
# 2. Create WHO activity category
# -------------------------------

dat_clean <- dat_clean %>%
  mutate(
    MET_MVPA_WHO = cut(
      cc_MET_MVPA_trunc,
      breaks = met_cuts,
      labels = met_labels,
      include.lowest = TRUE,
      right = FALSE
    )
  )

# Quick checks
sum(duplicated(dat_clean$eid))
dplyr::n_distinct(dat_clean$eid)
nrow(dat_clean)

# Save dataset with derivations for analysis


saveRDS(
  dat_clean,
  file.path(DATA_DERIVED, "analysis_dat.rds")
)


########################### Measures for time-to-event analysis #############################

# Time to event measure
# Time to fracture
# Person time at risk 


# Additional derivations from baseline dataset

dat_clean <- dat_clean %>%
  mutate(date_end_accel_derived = as.Date(date_end_accel_raw))

head(dat_clean$date_end_accel_derived, 10)

dat_clean <- dat_clean %>%
  mutate(
    # Age at accelerometer wear
    age_accel_entry_days_derived = as.numeric(
      difftime(date_end_accel_raw, approx_dob_derived, units = "days")
    ),
    age_accel_entry_years = age_accel_entry_days_derived / 365.25,
    
    # Age group at accelerometer wear
    age_gp_accel = cut(
      age_accel_entry_years,
      breaks = c(40, 50, 60, 70, 80),
      right = FALSE,
      labels = c("40-49", "50-59", "60-69", "70-79")
    )
  )


#### Derived accelerometry measures for creating catagorical WHO variables later

# Minutes/wk = 7 days/wk * 24 hours/day * 60 minutes/hr = 10080 minutes/wk

dat_clean <- dat_clean %>%
  mutate(
    light_accel_minwk = if_else(!is.na(light_average_raw), light_average_raw * 10080, NA_real_),
    sed_accel_minwk   = if_else(!is.na(sed_average_raw), sed_average_raw * 10080, NA_real_),
    MVPA_accel_minwk  = if_else(!is.na(mod_average_raw), mod_average_raw * 10080, NA_real_)
  )

# Create a total time group from all accel data ie sed, light and mvpa

dat_clean <- dat_clean %>%
  mutate(
    total_time_accel_minwk = light_accel_minwk +
      MVPA_accel_minwk +
      sed_accel_minwk
  )

# Create a total activity group from light and MVPA (not sed)

dat_clean <- dat_clean %>%
  mutate(
    total_activity_accel_minwk = light_accel_minwk + MVPA_accel_minwk
  )


# Rename raw accel catagories for clarity
dat_clean <- dat_clean %>%
  rename(
    proportion_light_accel = light_average_raw,
    proportion_mvpa_accel  = mod_average_raw,
    proportion_sed_accel   = sed_average_raw
  )

# Who guidelines on mins/wk of MVPA

min_cuts <- c(-Inf, 150, 300, Inf)

min_labels <- c("Below_guidelines",
                "Meets_guidelines",
                "Exceeds_guidelines")

# Create WHO MVPA accel group
dat_clean <- dat_clean %>%
  mutate(
    WHO_MVPA_accel_group = cut(
      MVPA_accel_minwk,
      breaks = min_cuts,
      labels = min_labels,
      right = FALSE
    )
  )


##### Load hes and death datasets to create censoring dates

death_clean <- readRDS(file.path(DATA_DERIVED, "death_clean.Rds"))
hes_clean   <- readRDS(file.path(DATA_DERIVED, "hes_clean.Rds"))

# -----------------------------
# 1. Keep dat_clean participant-level
# -----------------------------

death_clean <- death_clean %>%
  mutate(eid = as.character(eid))

hes_clean <- hes_clean %>%
  mutate(eid = as.character(eid))

dat_clean <- dat_clean %>%
  dplyr::left_join(
    death_clean %>% dplyr::select(eid, date_of_death),
    by = "eid"
  )

# -----------------------------
# 2. Create temporary dataset for HES event derivation
# -----------------------------

hes_for_tte <- hes_clean %>%
  dplyr::left_join(
    dat_clean %>%
      dplyr::select(eid, date_assess_A0_raw, date_end_accel_raw),
    by = "eid"
  )


# -----------------------------
# 3. Convert dates
# -----------------------------

hes_for_tte <- hes_for_tte %>%
  mutate(
    date_injury = as.Date(date_injury),
    date_assess_A0_raw = as.Date(date_assess_A0_raw),
    date_end_accel_raw = as.Date(date_end_accel_raw)
  )

# -----------------------------
# 4. Define codes
# -----------------------------


wrist_codes <- c("S52.5", "S52.6", "S62.0")
hip_codes   <- c("S72.0", "S72.1", "S72.2")

Ortho_codes <- c(
  "S02\\.[0-46-9]",
  "S12\\.[0-2]|S12\\.[7-9]",
  "S22\\.[0-5]|S22\\.[8-9]",
  "S32\\.[0-5]|S32\\.[7-8]",
  "S42\\.[0-4]|S42\\.[7-9]",
  "S52\\.[0-9]",
  "S62\\.[0-8]",
  "S72\\.[0-4]|S72\\.[7-9]",
  "S82\\.[0-9]",
  "S92\\.[0-5]"
)

ortho_pattern <- paste(Ortho_codes, collapse = "|")

fragility_codes <-c("S42.2", "S42.3", "S52.2", "S52.3", "S52.4", "S52.5", 
                    "S52.6", "S72.0", "S72.1", "S72.2", "S22.0", "S22.1", "S32.0", "S32.7")

# -----------------------------
# 4b. Reverse causation sensitivity:
# prior fragility fracture within 1 year before A0
# -----------------------------

prior_fragility_1yr <- hes_for_tte %>%
  dplyr::filter(
    !is.na(eid),
    !is.na(date_injury),
    !is.na(date_assess_A0_raw),
    diag_icd10_trimmed %in% fragility_codes,
    date_injury >= date_assess_A0_raw - 365,
    date_injury < date_assess_A0_raw
  ) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    date_prior_fragility_1yr = min(date_injury),
    prior_fragility_event_1yr = 1L,
    .groups = "drop"
  )

dat_clean <- dat_clean %>%
  dplyr::left_join(prior_fragility_1yr, by = "eid") %>%
  dplyr::mutate(
    prior_fragility_event_1yr = dplyr::if_else(
      is.na(prior_fragility_event_1yr),
      0L,
      prior_fragility_event_1yr
    )
  )

table(dat_clean$prior_fragility_event_1yr, useNA = "ifany")

# -----------------------------
# 5. Function to derive first events
# -----------------------------

# get first event after origin (which will be accel or A0)

get_first_event_after_origin <- function(data, event_name, origin_var, codes = NULL, pattern = NULL) {
  
  out_name <- paste0("date_first_", event_name, "_", origin_var)
  origin_col <- if (origin_var == "A0") "date_assess_A0_raw" else "date_end_accel_raw"
  
  df <- data %>%
    dplyr::filter(!is.na(date_injury), !is.na(.data[[origin_col]]))
  
  if (!is.null(codes)) {
    df <- df %>% dplyr::filter(diag_icd10_trimmed %in% codes)
  }
  
  if (!is.null(pattern)) {
    df <- df %>% dplyr::filter(grepl(pattern, diag_icd10_trimmed))
  }
  
  df %>%
    dplyr::filter(date_injury > .data[[origin_col]]) %>%
    dplyr::group_by(eid) %>%
    dplyr::summarise(
      !!out_name := min(date_injury),
      .groups = "drop"
    )
}

# -----------------------------
# 6. Run event derivations
# -----------------------------

event_defs <- list(
  list(event_name = "wrist",   codes = wrist_codes, pattern = NULL),
  list(event_name = "hip",     codes = hip_codes,   pattern = NULL),
  list(event_name = "fragility", codes = fragility_codes, pattern = NULL),
  list(event_name = "allfrac", codes = NULL,        pattern = ortho_pattern)
)

origins <- c("A0", "accel")

event_date_dfs <- list()

for (ev in event_defs) {
  for (org in origins) {
    event_date_dfs[[paste(ev$event_name, org, sep = "_")]] <-
      get_first_event_after_origin(
        data = hes_for_tte,
        event_name = ev$event_name,
        origin_var = org,
        codes = ev$codes,
        pattern = ev$pattern
      )
  }
}

# -----------------------------
# 7. Merge ONE ROW PER EID back
# -----------------------------

dat_clean <- purrr::reduce(
  event_date_dfs,
  ~ dplyr::left_join(.x, .y, by = "eid"),
  .init = dat_clean
)


# Check after merge
sum(duplicated(dat_clean$eid))
dat_clean %>%
  dplyr::count(eid) %>%
  dplyr::filter(n > 1) %>%
  head(10)


# We now add the outcome: time to incident fracture event. Participants can either:

#   - be observed to have a fracture event during their time in the study.
#   - be censored without having had a recorded fracture event.

#[Censoring](https://www.publichealth.columbia.edu/research/population-health-methods/time-event-data-analysis#:~:text=This%20phenomenon%20is%20called%20censoring,participant%20experiences%20a%20different%20event) 
# Censoring may take place at whichever of the following comes first
# (1) date of death
# (2) At the end of records 
# (3) or at the date at which a particular participant was recorded to be lost-to-follow-up.

# Loss-to-follow-up for particular participants is recorded in [field 191]

# [The censoring dates (the end of records) can be obtained from UK Biobank](https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=Data_providers_and_dates). 
# We match participants to the appropriate censoring date based on the country in which they attended the baseline assessment centre.

# The HES data is censored at
# 31 March 2023 England
# 31 August 2022 Scotland
# 31 May 2022 Wales

# Set administrative censoring date for each assessment country

dat_clean <- dat_clean %>%
  mutate(
    date_admin_cens = case_when(
      assessment_country == "England" ~ as.Date("2023-03-31"),
      assessment_country == "Scotland" ~ as.Date("2022-08-31"),
      assessment_country == "Wales" ~ as.Date("2022-05-31")
    )
  )

# Define the final censoring date
dat_clean <- dat_clean %>%
  mutate(
    date_cens = pmin(date_of_death, lost_to_fu_raw, date_admin_cens, na.rm = TRUE)
  )

# We now add a date for end of follow up, which is either the date at which the participant was censored 
# or the date at which they experienced a fracture event in hospital records as long as this occurred before censoring 

# Add follow up variable
# i.e. same as censor date for participants without a hospital-recorded fracture diagnosis,
# event date for participants with hospital-recorded fracture diagnosis that falls within the study period
# Wrist and Hip fracture will be run separately so need separate variables

# Ensure all relevant date variables are Date objects
date_vars <- c(
  "date_cens",
  "date_first_wrist_A0", "date_first_wrist_accel",
  "date_first_hip_A0", "date_first_hip_accel",
  "date_first_fragility_A0", "date_first_fragility_accel",
  "date_first_allfrac_A0", "date_first_allfrac_accel",
  "date_assess_A0_raw", "date_end_accel_raw"
)

dat_clean <- dat_clean %>%
  mutate(across(any_of(date_vars), as.Date))

# Create follow-up variable loop
make_tte_vars <- function(data, event_name, origin_var) {
  
  origin_date <- if (origin_var == "A0") "date_assess_A0_raw" else "date_end_accel_raw"
  suffix <- if (origin_var == "A0") "PA" else "accel"
  
  first_date_col <- paste0("date_first_", event_name, "_", origin_var)
  fu_date_col    <- paste0("date_fu_", event_name, "_", suffix)
  event_col      <- paste0("event_", event_name, "_", suffix)
  time_col       <- paste0("fu_time_", event_name, "_", suffix)
  
  data %>%
    dplyr::mutate(
      
      # Follow-up end
      !!fu_date_col := dplyr::if_else(
        !is.na(.data[[first_date_col]]),
        pmin(.data[[first_date_col]], date_cens),
        date_cens
      ),
      
      # Event indicator (only if under observation)
      !!event_col := dplyr::if_else(
        .data[[origin_date]] <= date_cens,
        as.integer(!is.na(.data[[first_date_col]]) & .data[[first_date_col]] <= date_cens),
        NA_integer_
      ),
      
      # Follow-up time (only valid if under observation)
      !!time_col := dplyr::if_else(
        .data[[origin_date]] <= .data[[fu_date_col]],
        as.numeric(.data[[fu_date_col]] - .data[[origin_date]]),
        NA_real_
      )
    )
}

# apply it

for (ev in c("allfrac", "wrist", "hip", "fragility")) {
  for (org in c("A0", "accel")) {
    dat_clean <- make_tte_vars(dat_clean, ev, org)
  }
}


colnames(dat_clean)
sum(duplicated(dat_clean$eid))
dat_clean %>%
  dplyr::count(eid) %>%
  dplyr::filter(n > 1) %>%
  head(10)



# -----------------------------
# Quick checks
# -----------------------------

# Check event counts
table(dat_clean$event_wrist_PA, useNA = "ifany")
table(dat_clean$event_wrist_accel, useNA = "ifany")
table(dat_clean$event_hip_PA, useNA = "ifany")
table(dat_clean$event_hip_accel, useNA = "ifany")
table(dat_clean$event_fragility_PA, useNA = "ifany")
table(dat_clean$event_fragility_accel, useNA = "ifany")
table(dat_clean$event_allfrac_PA, useNA = "ifany")
table(dat_clean$event_allfrac_accel, useNA = "ifany")

# Alternatively, we might want to analyse the data using age as the timescale, so we add a variable for age at exit in days:

# -----------------------------
# Calculate age at exit for all frac, wrist and hip (both PA and accelerometer)
# -----------------------------
dat_clean <- dat_clean %>%
  mutate(
    # Age at end of follow-up based on all fracture
    age_exit_allfrac_PA     = as.numeric(difftime(date_fu_allfrac_PA, approx_dob_derived, units = "days")) / 365.25,
    age_exit_allfrac_accel  = as.numeric(difftime(date_fu_allfrac_accel, approx_dob_derived, units = "days")) / 365.25,
    
    # Age at end of follow-up based on wrist fracture
    age_exit_wrist_PA       = as.numeric(difftime(date_fu_wrist_PA, approx_dob_derived, units = "days")) / 365.25,
    age_exit_wrist_accel    = as.numeric(difftime(date_fu_wrist_accel, approx_dob_derived, units = "days")) / 365.25,
    
    # Age at end of follow-up based on hip fracture
    age_exit_hip_PA         = as.numeric(difftime(date_fu_hip_PA, approx_dob_derived, units = "days")) / 365.25,
    age_exit_hip_accel      = as.numeric(difftime(date_fu_hip_accel, approx_dob_derived, units = "days")) / 365.25,
    
    # Age at end of follow-up based on fragility fracture
    age_exit_fragility_PA         = as.numeric(difftime(date_fu_fragility_PA, approx_dob_derived, units = "days")) / 365.25,
    age_exit_fragility_accel      = as.numeric(difftime(date_fu_fragility_accel, approx_dob_derived, units = "days")) / 365.25,
    
    # Age at entry
    age_entry_PA            = as.numeric(difftime(date_assess_A0_raw, approx_dob_derived, units = "days")) / 365.25,
    age_entry_accel         = as.numeric(difftime(date_end_accel_raw, approx_dob_derived, units = "days")) / 365.25
  )


colnames(dat_clean)
sum(duplicated(dat_clean$eid))
dat_clean %>%
  dplyr::count(eid) %>%
  dplyr::filter(n > 1) %>%
  head(10)

# We noted that it is well worth inspecting your data to check the code is behaving as expected, especially for some of the logically complex processes in this notebook. 

dat_clean$fu_years_fragility_PA <- dat_clean$fu_time_fragility_PA/365.25
# Follow up time distribution
hist(dat_clean$fu_years_fragility_PA, xlim = c(0,25))

dat_clean$fu_years_fragility_accel <- dat_clean$fu_time_fragility_accel/365.25
# Follow up time distribution
hist(dat_clean$fu_years_fragility_accel, xlim = c(0,25))


dat_clean$fu_years_wrist_PA <- dat_clean$fu_time_wrist_PA/365.25
# Follow up time distribution
hist(dat_clean$fu_years_wrist_PA, xlim = c(0,25))

dat_clean$fu_years_wrist_accel <- dat_clean$fu_time_wrist_accel/365.25
# Follow up time distribution
hist(dat_clean$fu_years_wrist_accel, xlim = c(0,25))

summary(dat_clean$fu_time_allfrac_PA)
sum(dat_clean$fu_time_allfrac_PA < 0, na.rm = TRUE)

summary(dat_clean$fu_time_allfrac_accel)
sum(dat_clean$fu_time_allfrac_accel < 0, na.rm = TRUE)

summary(dat_clean$age_exit_allfrac_accel)
summary(dat_clean$age_exit_allfrac_PA)
summary(dat_clean$age_exit_fragility_accel)
summary(dat_clean$age_exit_fragility_PA)
summary(dat_clean$age_exit_hip_accel)
summary(dat_clean$age_exit_hip_PA)
summary(dat_clean$age_exit_wrist_accel)
summary(dat_clean$age_exit_wrist_PA)




# Save dataset with derivations for Time to Event analysis

saveRDS(
  dat_clean,
  file.path(DATA_DERIVED, "TTE_analysis_dat.rds")
)

## Check number of events of wrist fracture self report and hip fracture HES

dat_clean %>%
  count(Wrist, event_hip_PA)

dat_clean %>%
  count(sex_raw, Wrist, event_hip_PA)

