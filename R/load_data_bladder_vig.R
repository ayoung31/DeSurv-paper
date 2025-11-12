load_data_bladder_vig = function(raw_data){

  dat = readRDS(raw_data)
  time = dat$sampInfo$os
  event = abs(1-dat$sampInfo$censOS)
  
  dat$sampInfo$time = time
  dat$sampInfo$event = event
  dat$sampInfo$dataset = "IMvigor"
  dat$sampInfo$ID = dat$sampInfo$sampleID
  
  rownames(dat$ex) = dat$featInfo
  dat$ex = log2(dat$ex+1)
  
  dat$samp_keeps = which(dat$sampInfo$Consensus != "NE-like")
  
  return(dat)
  
}