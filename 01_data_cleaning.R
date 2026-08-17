# This code describes the cleaning steps used on the raw data files for the PA and Fracture study.
# library and file paths are in 00_setup.R

source(here::here("scripts/00_setup.R"))


dat_raw <- fread(
  file.path(DATA_RAW, "participant_data.csv"),
  data.table = FALSE
)

### Remove participants that have requested to withdraw from the study ####

exclude_ids <- fread(
  file.path(DATA_RAW, "w81499_20260310.csv"),
  data.table = FALSE
)

# Assign column name
names(exclude_ids) <- "eid"

exclude_ids$eid <- as.character(exclude_ids$eid)
dat_raw$eid     <- as.character(dat_raw$eid)

dat_raw <- dat_raw %>%
  filter(!eid %in% exclude_ids$eid)

# Check excluded (should = 0)
sum(dat_raw$eid %in% exclude_ids$eid)


#### Data for cross sectional analysis ###################

# ----------------------------
# Helper Functions
# ----------------------------

# Clean character columns (trim, convert NA)
clean_character <- function(x, na_values = c("", "Do not know", "Prefer not to answer")) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% na_values] <- NA
  x
}

# Clean numeric columns (trim, convert NA, keep numeric strings)
clean_numeric <- function(x, na_values = c("", "Do not know", "Prefer not to answer")) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% na_values] <- NA
  numeric_mask <- str_detect(x, "^(\\d+\\.?\\d*|\\d*\\.\\d+)$")
  x[!numeric_mask] <- NA
  as.numeric(x)
}

# Clean age education (special case for "Never went to school")
clean_age_education <- function(x) {
  x <- clean_character(x)
  x[x == "Never went to school"] <- 0
  as.numeric(x)
}

# Apply clean_numeric to a set of columns
clean_columns <- function(data, vars, suffix = "_clean") {
  for (v in vars) {
    new_name <- sub("_raw$", suffix, v)
    data[[new_name]] <- clean_numeric(data[[v]])
  }
  data
}

check_completeness <- function(data, vars) {
  sapply(data[vars], function(x) c("NA" = sum(is.na(x)), "non_NA" = sum(!is.na(x))))
}


# Prepare final clean dataset
get_clean_dataset <- function(data) {
  all_cols <- names(data)
  
  clean_cols <- all_cols[grepl("_clean$", all_cols)]
  raw_cols   <- all_cols[grepl("_raw$", all_cols)]
  
  raw_cols_to_use <- raw_cols[
    !sub("_raw$", "_clean", raw_cols) %in% clean_cols
  ]
  
  other_cols <- all_cols[
    !grepl("_clean$|_raw$", all_cols)
  ]
  
  dplyr::select(data, all_of(c("eid", clean_cols, raw_cols_to_use, other_cols)))
}

# ----------------------------
# Rename columns (shortened for maintainability)
# ----------------------------

rename_map <- c(
  "eid" = "eid",
  "date_assess_A0_raw" = "Date of attending assessment centre | Instance 0",
  "ukb_assess_center_raw" = "UK Biobank assessment centre | Instance 0",
  "age_A0_raw" = "Age when attended assessment centre | Instance 0",
  "year_birth_raw" = "yob",
  "month_birth_raw" = "mob",
  "sex_raw" = "sex",
  "ethnicity_raw" = "Ethnic background | Instance 0",
  "qualification_raw" = "Qualifications | Instance 0",
  "age_education_raw" = "Age completed full time education | Instance 0",
  "tdi_raw" = "townsend_raw",
  "waist_circ_raw" = "Waist circumference | Instance 0",
  "hip_circ_raw" = "Hip circumference | Instance 0",
  "weight_raw" = "Weight | Instance 0",
  "height_raw" = "p50_i0",
  "BMI_raw" = "Body mass index (BMI) | Instance 0",
  "lost_to_fu_raw" = "lost_to_fu",
  "grip_strength_L_raw" = "Hand grip strength (left) | Instance 0",
  "grip_strength_R_raw" = "Hand grip strength (right) | Instance 0",
  "reaction_time_raw" = "Mean time to correctly identify matches | Instance 0",
  "fall_fracture_raw" = "Fracture resulting from simple fall | Instance 0",
  "self_reported_falls_raw" = "Falls in the last year | Instance 0",
  "self_reported_fracture_A0_raw" = "Fractured/broken bones in last 5 years | Instance 0",
  "fracture_location_A0_raw" = "Fractured bone site(s) | Instance 0",
  "duration_mod_raw" = "Duration of moderate activity | Instance 0",
  "duration_vig_raw" = "Duration of vigorous activity | Instance 0",
  "duration_walk_raw" = "Duration of walks | Instance 0",
  "num_days_mod_raw" = "Number of days/week of moderate physical activity 10+ minutes | Instance 0",
  "num_days_vig_raw" = "Number of days/week of vigorous physical activity 10+ minutes | Instance 0",
  "num_day_walk_raw" = "Number of days/week walked 10+ minutes | Instance 0",
  "usual_walking_pace_raw" = "Usual walking pace | Instance 0",
  "WHO_MVPA_raw" = "At or above moderate/vigorous recommendation | Instance 0",
  "WHO_walk_raw" = "At or above moderate/vigorous/walking recommendation | Instance 0",
  "IPAQ_group_raw" = "IPAQ Activity Group | Instance 0",
  "MET_mod_raw" = "MET minutes per week for moderate activity | Instance 0",
  "MET_vig_raw" = "MET minutes per week for vigorous activity | Instance 0",
  "MET_walk_raw" = "MET minutes per week for walking | Instance 0",
  "summed_MET_all_raw" = "Summed MET minutes per week for all activity | Instance 0",
  "summed_days_all_raw" = "Summed days activity | Instance 0",
  "Summed_mins_all_raw" = "Summed minutes activity | Instance 0",
  "light_average_raw" = "Light - Overall average | Instance 0",
  "mod_average_raw" = "Moderate-Vigorous - Overall average | Instance 0",
  "sed_average_raw" = "Sedentary - Overall average | Instance 0",
  "date_end_accel_raw" = "date_end_accel",
  "heel_BMD_raw" = "p3148_i0",
  "risk_taking_raw" = "p2040_i0",
  "drive_fast_raw" = "p1100_i0",
  "smoking_raw" = "p20116_i0",
  "etoh_raw" = "p20117_i0"
)

dat <- dat_raw %>% rename(!!!rename_map)

# ----------------------------
# Clean special character columns
# ----------------------------

# Education
dat$age_education_clean <- clean_age_education(dat$age_education_raw)

# Ethnicity
dat$ethnicity_clean <- clean_character(dat$ethnicity_raw)

# Fractures
dat$self_reported_fracture_A0_clean <- clean_character(dat$self_reported_fracture_A0_raw)
dat$fracture_location_A0_clean <- clean_character(dat$fracture_location_A0_raw)

# Falls
dat$self_reported_falls_clean <- clean_character(dat$self_reported_falls_raw)

# MET & IPAQ
dat$IPAQ_clean     <- clean_character(dat$IPAQ_group_raw)
dat$MET_mod_clean  <- clean_character(dat$MET_mod_raw)
dat$MET_vig_clean  <- clean_character(dat$MET_vig_raw)
dat$MET_walk_clean <- clean_character(dat$MET_walk_raw)
dat$WHO_MVPA_raw  <- clean_character(dat$WHO_MVPA_raw)
dat$WHO_walk_raw  <- clean_character(dat$WHO_walk_raw)

#additional variables
dat$smoking_clean <- clean_character(dat$smoking_raw)
dat$risk_taking_clean <- clean_character(dat$risk_taking_raw)
dat$etoh_clean <-clean_character(dat$etoh_raw)
dat$drive_fast_clean <-clean_character(dat$drive_fast_raw)


# ----------------------------
# Clean numeric columns
# ----------------------------

## Physical activity variables
pa_vars <- c(
  "duration_mod_raw", "num_days_mod_raw",
  "duration_vig_raw", "num_days_vig_raw",
  "duration_walk_raw", "num_day_walk_raw"
)

dat <- clean_columns(dat, pa_vars)

# Function to clean day-duration pairs
clean_pa_pair <- function(days, duration) {
  days_num <- suppressWarnings(as.numeric(days))
  duration_num <- suppressWarnings(as.numeric(duration))
  
  tibble(
    days_clean = case_when(
      is.na(days_num) ~ NA_real_,
      TRUE ~ days_num
    ),
    duration_clean = case_when(
      is.na(days_num) ~ NA_real_,                 # if days missing, both NA
      days_num == 0 ~ 0,                          # if days = 0, duration = 0
      days_num > 0 & is.na(duration_num) ~ NA_real_,  # if days > 0 and duration missing, NA
      TRUE ~ duration_num
    )
  )
}

# Apply once per pair
mod_clean  <- clean_pa_pair(dat$num_days_mod_raw,  dat$duration_mod_raw)
vig_clean  <- clean_pa_pair(dat$num_days_vig_raw,  dat$duration_vig_raw)
walk_clean <- clean_pa_pair(dat$num_day_walk_raw, dat$duration_walk_raw)

dat <- dat %>%
  mutate(
    num_days_mod_clean  = mod_clean$days_clean,
    duration_mod_clean  = mod_clean$duration_clean,
    
    num_days_vig_clean  = vig_clean$days_clean,
    duration_vig_clean  = vig_clean$duration_clean,
    
    num_day_walk_clean  = walk_clean$days_clean,
    duration_walk_clean = walk_clean$duration_clean
  )


# Anthropometrics
anthro_vars <- c("height_raw", "weight_raw", "waist_circ_raw", "BMI_raw", "heel_BMD_raw")
dat <- clean_columns(dat, anthro_vars)

# ----------------------------
# Check completeness (optional)
# ----------------------------
check_completeness(dat, c("duration_mod_clean", "num_days_mod_clean", "MET_mod_raw"))
check_completeness(dat, c("height_clean", "weight_clean", "waist_circ_clean", "BMI_clean"))

# ----------------------------
# Prepare final clean dataset
# ----------------------------
baseline_clean <- get_clean_dataset(dat)

#check values
dplyr::n_distinct(baseline_clean$eid)

# Save cleaned dataset
saveRDS(baseline_clean, file.path(DATA_DERIVED, "baseline_clean.rds"))

# check for duplicated eids
sum(duplicated(baseline_clean$eid))
dplyr::n_distinct(baseline_clean$eid)
nrow(baseline_clean)


################################################################################################################################

# Additional derivations needed for Time To Event (TTE) analysis


dat_death <- fread(
  file.path(DATA_RAW, "death_data.csv"),
  data.table = FALSE
)

# Death data
dat_death$date_death <-
  as.Date(dat_death$date_of_death, format = "%Y-%m-%d")

# A very small number of participants have duplicate records in death data (e.g. perhaps from a second death certificate after post-mortem)
# In this dataset we keep just one record per participant.

dat_death <-
  dat_death[dat_death$ins_index == 0, ]

# Save cleaned outputs
saveRDS(
  dat_death,
  file.path(DATA_DERIVED, "death_clean.rds")
)

##### HES #######


# HES data
# hes is a combination of admissions and diagnosis details merged prior to download. 
# These contain all of the ICD-10 codes for participants in the UKBiobank. 
# We have chosen not to use the ICD-9 codes as these cover dates prior to 1997 in just Wales and Scotland
# Load data

hes <- fread(
  file.path(DATA_RAW, "hes_data.csv"),
  data.table = FALSE
)

# Extract columns eid, epistart, admidate and ICD10_CODE into new df called hes_diag

hes_diag <- hes %>%
  dplyr::select(eid, epistart, admidate,date_injury, diag_icd10)

# Event data
# I have created a 'date_injury' at which the fracture event occured. For this I have used the admissions date (admidate).
# Where this was not available I have used the episode start date (epistart).

# Select out only the relevant orthopaedic codes based on the pre-specified analysis plan.

# Filter rows where ICD10_Code starts with 'S'
hes_filtered <- hes_diag %>%
  filter(
    grepl("^S", diag_icd10) 
  )

#UKB removes the decimal place, for filtering out of non fracture ortho codes we need to put it back in.

# Helper function: insert decimal after 3 characters if longer than 3
insert_icd10_decimal <- function(x) {
  x <- as.character(x)          # ensure character
  is_na <- is.na(x)             # keep track of NAs
  n <- nchar(x)
  
  # Only modify non-NA values longer than 3
  long_idx <- which(!is_na & n > 3)
  x[long_idx] <- paste0(
    substr(x[long_idx], 1, 3),
    ".",
    substr(x[long_idx], 4, n[long_idx])
  )
  
  x
}

# Apply to your filtered HES dataframe
hes_filtered <- hes_filtered %>%
  mutate(
    diag_icd10_full = insert_icd10_decimal(diag_icd10)
  )


hes_trimmed <- hes_filtered %>%
  mutate(
    diag_icd10_trimmed = sub("\\.([0-9])[0-9]+$", ".\\1", diag_icd10_full)
    
  )


# Define the Ortho_codes regex patterns
Ortho_codes <- c(
  "S02\\.[0-9]",  # frac_head not tooth
  "S12\\.[0-2]|S12\\.[7-9]",  # frac_Cx
  "S22\\.[0-5]|S22\\.[8-9]",  # frac_Tx
  "S32\\.[0-5]|S32\\.[7-8]",  # frac_Lx_pelvis
  "S42\\.[0-4]|S42\\.[7-9]",  # frac_Sh_arm
  "S52\\.[0-9]",  # frac_Elb_farm
  "S62\\.[0-8]",  # frac_Hand_carpal
  "S72\\.[0-4]|S72\\.[7-9]",  # frac_Hip_thigh
  "S82\\.[0-9]",  # frac_Knee_LL
  "S92\\.[0-5]|S72\\.[7-9]"  # frac_Ankle_ft
)

# Filter rows based on matching ICD10_Code with regex patterns in Ortho_codes
hes_filtered <- hes_trimmed %>%
  filter(
    sapply(1:nrow(.), function(i) {
      code10 <- diag_icd10_trimmed[i]
      any(sapply(Ortho_codes, function(pattern) grepl(pattern, code10))) 
    })
  )

num_unique_eid <- length(unique(hes_filtered$eid))
print(num_unique_eid)


### find duplicates in diag_icd10 ###

hes_check_duplicates <- hes_filtered %>%
  group_by(eid, date_injury, diag_icd10_trimmed) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1)

# View the result
View(hes_check_duplicates)


#### remove duplicated data i.e exact matches that have been entered more than once for the same date

# Remove duplicates, keeping only the first occurrence
hes_filtered_no_duplicates <- hes_filtered %>%
  group_by(eid, date_injury) %>%
  filter(!duplicated(diag_icd10_trimmed)) %>%
  ungroup()

num_unique_eid <-length(unique(hes_filtered_no_duplicates$eid))
print(num_unique_eid)


# Washout 

# As patients are often re-admitted with the same condition, i.e a hip fracture can be re-admitted with a second fall, 
# the ICD-10 code can appear twice in quick succession. This will have an effect on any incidence calculations. 
# To avoid this it is usual practice to apply a washout period where any repeat of the code is removed. 
# This does run the risk of under reporting as a second fracture is common within a year of a first event but the washout approach is a pragmatic solution. 
# I have applied a one year washout period.

# Convert startdate to Date format (if not already)
HES_clean_ICD10 <- hes_filtered_no_duplicates %>%
  mutate(date_injury = as.Date(date_injury)) %>%
  
  # Group by EID and Ortho_code
  group_by(eid, diag_icd10_trimmed) %>%
  
  # Sort entries by startdate within each group
  arrange(date_injury, .by_group = TRUE) %>%
  
  # Apply washout period logic
  mutate(
    keep_entry = case_when(
      is.na(lag(date_injury)) ~ TRUE,  # Keep the first occurrence
      (date_injury - lag(date_injury)) > 365 ~ TRUE,  # Keep if more than 1 year apart
      TRUE ~ FALSE  # Exclude if within 1 year of the previous entry
    )
  ) %>%
  
  # Keep only the entries marked as TRUE
  filter(keep_entry) %>%
  
  # Drop temporary column
  dplyr::select(-keep_entry) %>%
  
  # Ungroup the data
  ungroup()

# Save cleaned HES

saveRDS(
  HES_clean_ICD10,
  file.path(DATA_DERIVED, "hes_clean.rds")
)


# data check
# To make sure this is the most recent HES set we will now group by injury date from most to least recent
HES_clean <- HES_clean_ICD10 %>%
  arrange(desc(date_injury))

          View(HES_clean)
          



