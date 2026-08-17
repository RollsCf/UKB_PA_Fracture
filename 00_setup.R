# This source file specifies the libraries and file paths requried for the project.

# 1. Load packages
packages <- c(
  "dplyr",
  "ggplot2",
  "tidyr",
  "readr",
  "dotenv",
  "here",
  "data.table",
  "ggplot2",
  "stringr",
  "purrr",
  "e1071",
  "gridExtra",
  "officer",
  "broom",
  "emmeans",
  "grid",
  "gtsummary",
  "lubridate",
  "fmsb",
  "gt",
  "survival",
  "epiDisplay",
  "splines",
  "rstatix",
  "patchwork",
  "cmprsk",
  "splines",
  "rms"

)

## nb table1 not loaded here as dependent on gdtools and crashes the script


# Install missing packages if needed
installed <- packages %in% rownames(installed.packages())
if(any(!installed)) install.packages(packages[!installed])

# Load packages
lapply(packages, library, character.only = TRUE)

# 2.  Load environment variables from config.env
dotenv::load_dot_env(here::here("config.env"))

DATA_RAW     <- here(Sys.getenv("DATA_RAW"))
DATA_DERIVED <- here(Sys.getenv("DATA_DERIVED"))
RESULTS_DIR  <- here(Sys.getenv("RESULTS_DIR"))

if (RESULTS_DIR == "") {
  stop("RESULTS_DIR environment variable not set")
}

TABLES_DIR  <- file.path(RESULTS_DIR, "tables")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TTE_TABLES_DIR  <- file.path(RESULTS_DIR, "TTE_tables")
TTE_FIGURES_DIR <- file.path(RESULTS_DIR, "TTE_figures")

dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TTE_TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TTE_FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# =========================
# Word table helper
# =========================
save_table_word <- function(df, table_number, title = NULL, folder_path = TABLES_DIR) {
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  file_path <- file.path(folder_path, paste0("Table_", table_number, ".docx"))
  
  doc <- officer::read_docx()
  
  if (!is.null(title)) {
    doc <- doc %>%
      officer::body_add_par(
        paste0("Table ", table_number, ": ", title),
        style = "Normal"
      ) %>%
      officer::body_add_par("", style = "Normal")
  }
  
  doc <- doc %>% officer::body_add_table(df, style = "table_template")
  print(doc, target = file_path)
  
  message("Saved: ", file_path)
  invisible(file_path)
}