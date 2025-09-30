# Reduce a tibble (one combo’s rows) to best init
select_best_init <- function(df, method_select="surv") {
  
  # only keep initializations with finite loss and no errors
  keep <- df[!df$flag_nan & is.finite(df$loss), ]
  
  # save method for init selection in dataframe
  df$method_select = method_select
  
  if(method_select=="surv"){
    best = df %>% dplyr::slice_max(order_by = sloss, n = 1, with_ties = FALSE)
  }else if(method_select=="nmf"){
    best = df %>% dplyr::slice_min(order_by = nloss, n = 1, with_ties = FALSE)
  }else if(method_select=="loss"){
    best = df %>% dplyr::slice_min(order_by = loss, n = 1, with_ties = FALSE)
  }else{
    stop("This method for selecting the best initialization is not supported")
  }
  
  return(best)
}
