# This script details the cross-sectional analysis looking at self-reported physical activity 
# and self-reported fracture risk using the data cleaned and derived in 01_data_cleaning.R and 02_data_derivation.R

source(here::here("scripts/00_setup.R"))
source(here("scripts/03_helpers_general.R"))
source(here("scripts/04_helpers_CXA.R"))
source(here("scripts/05_helpers_survival.R"))


# Load processed data from previous scripts
analysis_dat  <- readRDS(file.path(DATA_DERIVED, "analysis_dat.Rds"))

# Keep all variables but rename by dropping the words 'raw' or 'clean'
analysis_dat <- analysis_dat %>%
  rename_with(~ gsub("_clean|_raw", "", .x))

# rename variables for future analysis
analysis_dat <- analysis_dat %>%
  rename(
    SRF = self_reported_fracture_A0,
    WC = waist_circ
  )

# Count EID in full cohort prior to age restriction
analysis_dat %>%
  summarise(n_unique = n_distinct(eid))
    
# Create middle aged cohort A0 for all subsequent analysis

# have had issues with date columns so run this first

analysis_dat$age_A0 <- as.numeric(analysis_dat$age_A0)
analysis_dat <- analysis_dat %>%
  dplyr::mutate(
    date_assess_A0 = data.table::as.IDate(as.Date(date_assess_A0)),
    lost_to_fu     = data.table::as.IDate(as.Date(lost_to_fu))
  )

middle_aged_dat <- analysis_dat %>%
  filter(between(age_A0, 40, 65))



# =========================================================
# Apply exclusions to create complete-case dataset
# =========================================================


# Before working with the data we want to exclude those with missing data in covariates and pre-existing fracture
# We will record how many participants are excluded at each of the steps (e.g. for a flow diagram):

tab_exc <- data.frame("Exclusion" = "Starting cohort (age 40-65)",
                      "Number_excluded" = NA, 
                      "Number_remaining" = nrow(middle_aged_dat)
                      )


# First, we exclude participants with missing fracture data

# Store number before exclusion
nb <- nrow(middle_aged_dat)

# Exclude rows with NA in SRF
middle_aged_dat <- middle_aged_dat[!is.na(middle_aged_dat$SRF), ]

# Add to exclusion table
tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing fracture data",
    Number_excluded = nb - nrow(middle_aged_dat),
    Number_remaining = nrow(middle_aged_dat)
  )
)

# ---- Missing PA  data ----
nb <- nrow(middle_aged_dat)

middle_aged_dat <- middle_aged_dat[
  !is.na(middle_aged_dat$cc_MET_walk_trunc) &
    !is.na(middle_aged_dat$cc_MET_mod_trunc) &
    !is.na(middle_aged_dat$cc_MET_vig_trunc),
]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing activity data",
    Number_excluded = nb - nrow(middle_aged_dat),
    Number_remaining = nrow(middle_aged_dat)
  )
)

# ---- Missing ethnicity ----
nb <- nrow(middle_aged_dat)

middle_aged_dat <- middle_aged_dat[!is.na(middle_aged_dat$ethnicity_derived), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing Ethnicity data",
    Number_excluded = nb - nrow(middle_aged_dat),
    Number_remaining = nrow(middle_aged_dat)
  )
)

# ---- Missing deprivation ----
nb <- nrow(middle_aged_dat)

middle_aged_dat <- middle_aged_dat[!is.na(middle_aged_dat$tdi), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing deprivation data",
    Number_excluded = nb - nrow(middle_aged_dat),
    Number_remaining = nrow(middle_aged_dat)
  )
)

# ---- Missing education ----
nb <- nrow(middle_aged_dat)

middle_aged_dat <- middle_aged_dat[!is.na(middle_aged_dat$education_level), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing qualification data",
    Number_excluded = nb - nrow(middle_aged_dat),
    Number_remaining = nrow(middle_aged_dat)
  )
)

# ---- Missing anthropometrics ----
nb <- nrow(middle_aged_dat)

middle_aged_dat <- middle_aged_dat[
  !is.na(middle_aged_dat$weight) &
    !is.na(middle_aged_dat$height) &
    !is.na(middle_aged_dat$WC) &
    !is.na(middle_aged_dat$BMI),
]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing anthropometric data",
    Number_excluded = nb - nrow(middle_aged_dat),
    Number_remaining = nrow(middle_aged_dat)
  )
)

View(tab_exc)

# Save complete-case analysis dataset
saveRDS(
  middle_aged_dat,
  file.path(DATA_DERIVED, "complete_case_dat.rds")
)

# Save exclusion table for tables/figures script
saveRDS(
  tab_exc,
  file.path(DATA_DERIVED, "tab_exclusions.rds")
)




