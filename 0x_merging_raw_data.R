dat_raw <- fread("data_raw/participant_data.csv", data.table = FALSE)
additional <-fread("data_raw/additional_variables.csv", data.table = FALSE)

dat_merged <- dat_raw %>%
  left_join(
   additional %>% select(eid, p2040_i0, p1100_i0, p20116_i0, p20117_i0, p3148_i0),
    by = "eid"
  )

write.csv(
  dat_merged,
  "data_raw/participant_data_with_all_additional.csv",
  row.names = FALSE
)

## Manually rename back to participant_data for further analysis ##