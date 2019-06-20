
source("gpdfit.R")
source("psipw.R")

### Naive estimator
naive = function(A, Y){
  N_eff = length(A)
  ate = mean(Y[A==1]) - mean(Y[A==0])
  list(ate = ate, N_eff = N_eff)
}


### IPTW estimator
iptw = function(A, Y, w){
  N_eff = n_eff(wt = w)
  ate = mean(w * as.numeric(A==1) * Y) - mean(w * as.numeric(A==0) * Y)
  list(ate = ate, N_eff = N_eff)
}


### Truncate IPTW estimator
# x = sample(x = 1:4, size = 10, replace = TRUE)
# y = x[order(x, decreasing = FALSE)]
# quantile(y, probs = 0.8, names = FALSE)
# tail(x[order(x)], as.integer(0.2 * 10))[1]
iptw_truncate = function(A, Y, w, wcp = 0.2, ceiling = 0){
  if(ceiling == 0){
    S = length(w)
    M = min(as.integer(wcp*S), as.integer(3 * sqrt(S)))
    ceiling = tail(x = w[order(w)], n = M)[1]                           #  truncate largest M weights
    print(paste(c("ceiling = ", ceiling, "#_truncate = ", sum(w>ceiling))))
  }
  
  wtrunc = w
  wtrunc[wtrunc > ceiling] = ceiling
  iptw(A, Y, wtrunc)
}


### normalized IPTW
iptw_norm  = function(A, Y, w){
  denom0 = mean(w * as.numeric(A==0))
  denom1 = mean(w * as.numeric(A==1))
  wnorm = w * as.numeric(A==0) / denom0 + w * as.numeric(A==1) / denom1
    
  N_eff = n_eff(wt = wnorm)
  ate = mean(wnorm * as.numeric(A==1) * Y) - mean(w * as.numeric(A==0) * Y)
  # ate = mean(w * as.numeric(A==1) * Y)/denom1 - mean(w * as.numeric(A==0) * Y)/denom0
  list(ate = ate, N_eff = N_eff)
}


### Pareto-Smoothed Propensity weighting
iptw_ps  = function(A, Y, w, wcp = 0.2, min_cutoff = 20){
  fit = psipw(wt = w, wcp = wcp, MIN_CUTOFF = min_cutoff)
  ws = fit$wt_new
  print(paste(c("N_effective = ", n_eff(wt = ws))))
  
  iptw(A, Y, ws)
}


### Pareto-Smoothed Propensity weighting using IPTW-Norm
iptw_ps_norm  = function(A, Y, w, wcp = 0.2, min_cutoff = 20){
  fit = psipw(wt = w, wcp = wcp, MIN_CUTOFF = min_cutoff)
  ws = fit$wt_new
  
  iptw_norm(A, Y, ws)
}



