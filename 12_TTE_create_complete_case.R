# This script is to be used with the Time To Event (TTE) R markdown file. It uses the same data that was derived
# for the cross sectional analysis along with additional HES and accelerometry data.

library(here)
source(here::here("scripts/00_setup.R"))
source(here("scripts/03_helpers_general.R"))
source(here("scripts/04_helpers_CXA.R"))
source(here("scripts/05_helpers_survival.R"))


# Load processed data from previous scripts
TTE_analysis_dat  <- readRDS(file.path(DATA_DERIVED, "TTE_analysis_dat.Rds"))


check_eid(TTE_analysis_dat)

# have had issues with date columns so run this first

TTE_analysis_dat$age_A0 <- as.numeric(TTE_analysis_dat$age_A0_raw)
TTE_analysis_dat <- TTE_analysis_dat %>%
  dplyr::mutate(
    date_assess_A0_raw = data.table::as.IDate(as.Date(date_assess_A0_raw)),
    lost_to_fu_raw     = data.table::as.IDate(as.Date(lost_to_fu_raw)),
    date_of_death     = data.table::as.IDate(as.Date(date_of_death))
  )

# We will be running 2 separate time to event analysis 
# (1) Using self-report PA (TTE_A0)
# (2) Using accelerometyer PA (TTE_accel)
# Both will use outcome data from HES.

# Each analysis will require its own cohort with exclusions.

####### Analysis one Self-reported PA and HES fracture ############
###################################################################

# Create middle aged cohort A0 for all subsequent analysis of self-reported activity

TTE_A0 <- TTE_analysis_dat %>%
  filter(between(age_A0, 40, 65))

# Ensure sex is a factor
TTE_A0$sex <- factor(TTE_A0$sex)

TTE_A0 <- TTE_A0 %>%
  mutate(
    Wrist = factor(event_wrist_PA, levels = c(0, 1), labels = c("No", "Yes")),
    Hip   = factor(event_hip_PA,   levels = c(0, 1), labels = c("No", "Yes")),
    fragility = factor(event_fragility_PA,   levels = c(0, 1), labels = c("No", "Yes"))
  )

check_eid(TTE_A0)



## Create a complete case analysis data set for primary analysis

## Exclusions 

# Before working with the data we want to exclude those with missing data in covariates and pre-existing fracture
# We will record how many participants are excluded at each of the steps (e.g. for a flow diagram):

tab_exc <- data.frame("Exclusion" = "Starting cohort self-reported PA (age 40-65)",
                      "Number_excluded" = NA, 
                      "Number_remaining" = nrow(TTE_A0)
)

View (tab_exc)


# ---- Missing PA  data ----

nb <- nrow(TTE_A0)

TTE_A0<- TTE_A0[
  !is.na(TTE_A0$cc_MET_walk_trunc) &
    !is.na(TTE_A0$cc_MET_mod_trunc) &
    !is.na (TTE_A0$cc_MET_vig_trunc),
  
]


tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing activity data",
    Number_excluded = nb - nrow(TTE_A0),
    Number_remaining = nrow(TTE_A0)
  )
)

# ---- Missing ethnicity ----
nb <- nrow(TTE_A0)

TTE_A0 <- TTE_A0[!is.na(TTE_A0$ethnicity_derived), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing Ethnicity data",
    Number_excluded = nb - nrow(TTE_A0),
    Number_remaining = nrow(TTE_A0)
  )
)

# ---- Missing deprivation ----
nb <- nrow(TTE_A0)

TTE_A0 <- TTE_A0[!is.na(TTE_A0$tdi_raw), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing deprivation data",
    Number_excluded = nb - nrow(TTE_A0),
    Number_remaining = nrow(TTE_A0)
  )
)

# ---- Missing education ----
nb <- nrow(TTE_A0)

TTE_A0 <- TTE_A0[!is.na(TTE_A0$education_level), ]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing qualification data",
    Number_excluded = nb - nrow(TTE_A0),
    Number_remaining = nrow(TTE_A0)
  )
)

# ---- Missing anthropometrics ----
nb <- nrow(TTE_A0)

TTE_A0 <- TTE_A0[
  !is.na(TTE_A0$weight_clean) &
    !is.na(TTE_A0$height_clean),
   
]

tab_exc <- rbind(
  tab_exc,
  data.frame(
    Exclusion = "Missing anthropometric data",
    Number_excluded = nb - nrow(TTE_A0),
    Number_remaining = nrow(TTE_A0)
  )
)

check_eid(TTE_A0)

# check PA variables complete
sum(is.na(TTE_A0$cc_MET_total_trunc))
sum(is.nan(TTE_A0$cc_MET_total_trunc))
sum(is.infinite(TTE_A0$cc_MET_total_trunc))

# Save exclusion table for tables/figures script
saveRDS(
  tab_exc,
  file.path(DATA_DERIVED, "TTE_tab_exclusions.rds")
)

View(tab_exc)


# Save cleaned outputs
saveRDS(
  TTE_A0,
  file.path(DATA_DERIVED, "TTE_A0_complete_case.rds")
)










