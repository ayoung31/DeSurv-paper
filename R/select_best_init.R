# Reduce a tibble (one combo’s rows) to best init
select_best_init <- function(df, method_select="surv") {
  if(method_select=="surv"){
    df$loss = df$surv_loss
  }else if(method_select=="nmf"){
    df$loss = df$nmf_loss
  }else{
    stop("This method for selecting the best initialization is not supported")
  }
  df$method_select = method_select
  
  keep <- df[!df$flag_nan & is.finite(df$loss), ]
  
  best = df %>% dplyr::slice_min(order_by = loss, n = 1, with_ties = FALSE)
  
  return(best$seed)
}