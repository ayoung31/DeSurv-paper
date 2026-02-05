is_valid_fit = function(path,alpha,ninit){
  if (file.exists(path)) {
    bundle_old = readRDS(path)
    exists = unlist(lapply(bundle_old,function(x) names(x$fits)==as.character(alpha)))
    
    if(all(exists) && (length(bundle_old)==ninit)){
      return(TRUE)
    }
  }
  return(FALSE)
}
