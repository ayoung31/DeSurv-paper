load_data_bladder_vig = function(raw_data){

  dat = readRDS(raw_data)
  time = dat$sampInfo$os
  event = dat$sampInfo$censOS

  dat$sampInfo$time = time
  dat$sampInfo$event = event
  dat$sampInfo$dataset = "IMvigor"
  dat$sampInfo$ID = dat$sampInfo$sampleID

  rownames(dat$ex) = dat$featInfo
  dat$ex = log2(dat$ex+1)
  dat$dataname = "IMvigor"

  dat$samp_keeps = which(dat$sampInfo$Consensus != "NE-like")

  return(dat)

}

load_data_bladder_vig_nobcg = function(raw_data){

  dat = load_data_bladder_vig(raw_data)

  bcg_flag <- dat$clinical$Intravesical.BCG.administered
  nobcg_idx <- which(bcg_flag == "N")
  dat$samp_keeps = intersect(dat$samp_keeps, nobcg_idx)

  return(dat)

}
