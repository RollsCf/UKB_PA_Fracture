# The following helper functions can be used across multiple scripts

# Create a function to check that during process NA eids not introduced
check_eid <- function(data) {
  c(
    total_rows = nrow(data),
    non_missing_eid = sum(!is.na(data$eid)),
    unique_eid = dplyr::n_distinct(data$eid, na.rm = TRUE),
    duplicated_eid = sum(duplicated(data$eid[!is.na(data$eid)])),
    missing_eid = sum(is.na(data$eid))
  )
}


# =========================================================
# Missing activity data comparison
# =========================================================

get_smd <- function(var, data) {
  
  x <- data[[var]]
  g <- data[[group_var]]
  
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- droplevels(g[keep])
  
  if (is.numeric(x)) {
    
    x1 <- x[g == "Has activity data"]
    x2 <- x[g == "Missing activity data"]
    
    smd <- (mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)) /
      sqrt((sd(x1, na.rm = TRUE)^2 + sd(x2, na.rm = TRUE)^2) / 2)
    
    smd <- abs(smd)
    
  } else {
    
    x <- factor(x)
    
    smds <- sapply(levels(x), function(lvl) {
      
      p1 <- mean(x[g == "Has activity data"] == lvl)
      p2 <- mean(x[g == "Missing activity data"] == lvl)
      
      denom <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
      
      if (denom == 0) {
        return(NA_real_)
      } else {
        return(abs((p1 - p2) / denom))
      }
    })
    
    smd <- max(smds, na.rm = TRUE)
  }
  
  tibble(
    variable = var,
    SMD = round(smd, 3)
  )
}
