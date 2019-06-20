#' Generalized Pareto distribution (GPD) parameter estimation
#' p(x|k, sigma)
#' 
#' Estimate the parameters p(x|k, sigma) given a sample (x-mu)
#' 
#' @adaptedFrom https://github.com/stan-dev/loo/blob/master/R/gpdfit.R
#' @references 
#' Zhang, J. & Stephens, M.A. 2009, 'A new and efficient estimation method for 
#' the generalized Pareto distribution', Technometrics, vol. 51, no. 3, pp. 316-25.
gpdfit = function(x){
  N = length(x)
  x = sort.int(x = x, decreasing = FALSE, method = "quick")
  prior = 3                   
  M = 80 + floor(sqrt(N))    
  mseq = seq_len(M)           # 1,2,3, ..., M
  sM = 1 - sqrt(M/(mseq-0.5))
  Nflr = floor(N/4 + 0.5)
  b = 1/x[N] + (sM/prior)/x[Nflr]
  l = N * lx(b, x)
  w = 1 / vapply(X = mseq, FUN = function(j){sum(exp(l - l[j]))}, FUN.VALUE = 0)
  bdotw = sum(b * w)
  k = mean(log1p(-bdotw * x))
  sigma = -k / bdotw
  list(k = k, sigma = sigma)
}


# internal --------
lx = function(a, x){
  a = -a
  k = sapply(X = a, FUN = function(y){mean(log1p(y * x))})
  log(a/k) - k -1
}

