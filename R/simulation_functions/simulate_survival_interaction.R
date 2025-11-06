simulate_survival_interaction <- function(H,
                                          baseline_hazard = 0.01,
                                          censor_rate = 0.005,
                                          pair = c(1,3),
                                          gamma = 1.5) {
  # pair = which programs interact to kill you
  h1 <- H[pair[1], ]
  h2 <- H[pair[2], ]
  
  # nonlinear risk: high only if both high
  linpred <- gamma * (h1 * h2)
  
  event_time <- rexp(length(linpred),
                     rate = baseline_hazard * exp(linpred))
  censor_time <- rexp(length(linpred), rate = censor_rate)
  
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)
  
  data.frame(
    patient = colnames(H),
    time = time,
    status = status,
    linpred_true = linpred  # true risk driver is interaction
  )
}