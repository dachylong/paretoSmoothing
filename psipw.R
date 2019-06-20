
source("./gpdfit.R")
source("./utility.R")
#' Pareto smoothed inverse propensity weighting (PSIPW)
#' 
#' @param wt 
#' @param wcp           cut-point percentage. default (0.2) = top 20% largest weights
#' @param MIN_CUTOFF 
psipw = function(wt, wcp = 0.2, MIN_CUTOFF = 20){
  # minimal cutoff value, there must >= 5 ipw larger than this 
  # in order to fit the GPD to the tail
  MIN_TAIL_LENGTH = 5
  N = length(wt)
  x = wt                       
  
  mu = wt_cutpoint(x, wcp, MIN_CUTOFF)  # 
  above_cut = (x > mu)
  x_body = x[!above_cut]                    # [x_min, mu]
  x_tail = x[above_cut]                     # upper tail: (mu, x_max]
  
  tail_len = length(x_tail)
  print(paste0("mu = ", mu, "smoothing length = ", tail_len, "/", N))
  x_new = x
  if (tail_len < MIN_TAIL_LENGTH || all(x_tail == x_tail[1])){    # No smoothing
    if (all(x_tail == x_tail[1]))
      warning("All tail values are the same. ", "Weights are truncated but not smoothed", call. = FALSE)
    else if (tail_len < MIN_TAIL_LENGTH)
      warning("Too few tail samples to fit GPD. \n", "Weights are truncated but not smoothed", call. = FALSE)
    k = Inf
  }
  else {                                                          # Smoothing
    tail_ord = order(x_tail)        # [86, 7, 64, 110, 55, 82, 8, 77, 4, 112, 46]
    fit = gpdfit(x_tail - mu)
    k = fit$k
    sigma = fit$sigma
    print(paste(c("mu = ", mu, "k = ", k, "sigma = ", sigma)))
    qq = qgpd(p = (seq_len(tail_len) - 0.5) /tail_len, mu = mu, k = k, sigma = sigma)
    smoothed_tail = rep.int(0, tail_len)
    smoothed_tail[tail_ord] = qq
    x_new[!above_cut] = x_body
    x_new[above_cut] = smoothed_tail
  }
  
  list(wt_new = x_new, k = k)
}


# internal -------------------------------------------------------
# Decising mu: the cutpoint of the log_weights ---
wt_cutpoint = function(y, wcp, min_cut){          
  if (min_cut > .Machine$double.xmax)
    min_cut = 800
  
  # cp = quantile(y, probs = 1-wcp, names = FALSE)         # 80% quantile point for wcp = 0.2
  S = length(y)
  M = min(as.integer(wcp*S), as.integer(3 * sqrt(S)))      # M = min(0.2S, 3sqrt(S))
  cp = tail(y[order(y)], n = M)[1]                        # 20% largest point
  
  max(cp, min_cut)
}




