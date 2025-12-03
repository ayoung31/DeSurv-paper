aggregate_summary_tables <- function(summary_list) {
  if (!length(summary_list)) {
    return(data.frame())
  }
  tables <- lapply(summary_list, `[[`, "summary")
  tables <- tables[lengths(tables) > 0]
  if (!length(tables)) {
    return(data.frame())
  }
  do.call(rbind, tables)
}
