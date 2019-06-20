
### Calculate the effective size of a weight vector ------
n_eff = function(wt){
  wt = abs(wt)
  n = length(wt)
  return(n * mean(wt)^2 / mean(wt^2))
}

### 
MAE = function(ate_true, ate_estimate){
  return(abs(ate_true-ate_estimate))
}

### Standardized (mean) difference of each covariate between treated v.s. control ----------
stand_meanDiff = function(X, A){
  treated = as.numeric(A==1)
  control = as.numeric(A==0)
  (mean(X[treated,]) - mean(X[control,])) / sqrt((var(X[treated,]) - var(X[control,]))/2)
}


### Plot the ecdf and weighted ecdf for each pre-treatment covariate ------
#' to demonstrate the effect of covariate balancing and graphical diagnostics
# install.packages("ATE")
# refer ATE::plot()
# plot.ecdf.x = function(X, wt){
#   
# }


### Inverse-CDF of GPD (from Wikipedia) --------
#' p : vector of probabilities
#' mu: location / lower bound
#' k : shape
#' sigma: scale
#' 
# qgpd(p, mu, xi, beta, lower.tail)  in fExtreme
# 
# example:
#     p = (1:20 - 0.5) / 20
#     qgpd(p)
#     fExtremes::qgpd(p = p)
qgpd = function(p, mu = 0, k = 1, sigma = 1, lower.tail = TRUE){
  if (is.nan(sigma) || sigma < 0)
    return(rep(NaN, length(p)))
  if (!lower.tail)
    p = (1-p)
  
  mu + sigma * ((1-p)^(-k) - 1) / k
}


