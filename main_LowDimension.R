
rm(list = ls())
library(reshape2)
library(ggplot2)
source("estimators.R")

genData = function(N){
  # simulation section 3.1 in Estimating Population Average Causal Effects in the Presence of Non-Overlap: A Bayesian Approach
  extra = 1.15
  n1 = n0 = N/2
  
  A = c(rep(1,n1), rep(0,n0))
  X1 = c(rbinom(n = n1, size = 1, prob = 0.5), rbinom(n = n0, size = 1, prob = 0.4))
  X2 = c(rnorm(n = n1, mean = 1, sd = 1), rnorm(n = n0, mean = -1, sd = 1))
  Y0 = 0.85*X2 + 0.05*X2^2 + 2
  Y1 = 1/(1+exp(-(8*X2 - 1))) + 0.25*X1   ## true potential outcomes & causal effect
  Y = Y1*A + Y0*(1-A)                     ## observed outcome
  data.frame(X1,X2,A,Y,Y1,Y0)
}

# MCMC true ATE   ~= 0.518
set.seed(2017)
ATE_hat = rep(0, 100)
for (N in 1:100){
  data = genData(100000)
  ATE_hat[N] = mean(data$Y1) - mean(data$Y0)
}
ate_true = mean(ATE_hat)
rm(data)
rm(ATE_hat)

set.seed(2018)
Ns = c(100, 200, 500, 1000, 2000, 5000)
Ns = c(50, 100, 200, 300, 500, 1000)
Ns = c(100, 200, 300, 500, 1000, 2000)
reps = 20                                 # no. of replications to be implemented for calculating confidence interval.
MAEs = ESSs = list()

for (i in 1:length(Ns)){
  N = Ns[i]
  print(paste0("Starting experiments with N = ", N, " ---------------"))
  
  # N = Ns[2]
  mae_i = data.frame(NO=rep(N,reps), Naive=rep(0,reps), IPW=rep(0,reps), Trunc=rep(0,reps), Norm=rep(0,reps), PS=rep(0,reps), PSNorm=rep(0,reps))
  ess_i = data.frame(NO=rep(N,reps), Naive=rep(0,reps), IPW=rep(0,reps), Trunc=rep(0,reps), Norm=rep(0,reps), PS=rep(0,reps), PSNorm=rep(0,reps))
  for (k in 1:reps) {
    print(paste0("Replicaiton ", k))
    data = genData(N = N)      ## generate data
    
    ## estimate propensity scores -> IP weights
    emod = glm(formula = A~as.factor(X1)+X2, family = "binomial", data = data)
    ps  = predict(emod, type = "response")
    pAW = ps*as.numeric(data$A==1) + (1-ps)*as.numeric(data$A==0)
    wt  = 1/pAW
    
    ## 
    ate_naive  = naive(A = data$A, Y = data$Y)
    ate_ipw    = iptw(A = data$A, Y = data$Y, w = wt)
    ate_trunc  = iptw_truncate(A = data$A, Y = data$Y, w = wt, ceiling = 0)
    ate_norm   = iptw_norm(A = data$A, Y = data$Y, w = wt)
    ate_ps     = iptw_ps(A = data$A, Y = data$Y, w = wt, wcp = 0.2, min_cutoff = 1)
    ate_psnorm = iptw_ps_norm(A = data$A, Y = data$Y, w = wt, wcp = 0.2, min_cutoff = 1)
    
    mae_naive  = MAE(ate_true = ate_true, ate_estimate = ate_naive$ate)
    mae_ipw    = MAE(ate_true = ate_true, ate_estimate = ate_ipw$ate)
    mae_trunc  = MAE(ate_true = ate_true, ate_estimate = ate_trunc$ate)
    mae_norm   = MAE(ate_true = ate_true, ate_estimate = ate_norm$ate)
    mae_ps     = MAE(ate_true = ate_true, ate_estimate = ate_ps$ate)
    mae_psNorm = MAE(ate_true = ate_true, ate_estimate = ate_psnorm$ate)
    
    mae_i[k, ] = c(N, mae_naive, mae_ipw, mae_trunc, mae_norm, mae_ps, mae_psNorm)
    ess_i[k, ] = c(N, ate_naive$N_eff, ate_ipw$N_eff, ate_trunc$N_eff, ate_norm$N_eff, ate_ps$N_eff, ate_psnorm$N_eff)
    MAEs[[i]] = mae_i
    ESSs[[i]] = ess_i
  }
  
  print(paste0("Experiments with N = ", N, " Done!"))
}


#### Parse results & visualization. ------------
for (i in 1:length(Ns)){
  MAEs[[i]] = melt(MAEs[[i]], id.var = "NO")
  ESSs[[i]] = melt(ESSs[[i]], id.var = "NO")
}

MAEs.melt = MAEs[[1]]
ESSs.melt = ESSs[[1]]
for (i in 2:length(Ns)){
  MAEs.melt = rbind(MAEs.melt, MAEs[[i]])
  ESSs.melt = rbind(ESSs.melt, ESSs[[i]])
}

ESSs.melt[,3] = log(ESSs.melt[,3])
names(MAEs.melt) = c("Sample_Size", "Estimator", "Value")
names(ESSs.melt) = c("Sample_Size", "Estimator", "Value")

### get MAE comparison Table
library(dplyr)

MAEs.melt %>% 
  group_by(Sample_Size, Estimator) %>%
  summarise(mean=mean(Value), sd=sd(Value), n=n()) %>%
  mutate(se = sd / sqrt(n)) ->
  MAEs.comparison 

MAEs.comparison %>% 
  ggplot(aes(x = factor(Sample_Size), y = mean)) +
  geom_bar(aes(fill = Estimator), stat = "identity", position = position_dodge(0.9)) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "MAE") 
  
write.csv(as.data.frame(MAEs.comparison), file = "./result_lowDim.csv")


ESSs.melt %>% 
  group_by(Sample_Size, Estimator) %>%
  summarise(mean=mean(Value), sd=sd(Value), n=n()) %>%
  mutate(se = sd / sqrt(n)) ->
  ESSs.comparison 

ESSs.comparison %>% 
  ggplot(aes(x = factor(Sample_Size), y = mean)) +
  geom_bar(aes(fill = Estimator), stat = "identity", position = position_dodge(0.9)) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "log-ESS") 

ggplot(data = MAEs.melt, aes(x = factor(Sample_Size), y = Value)) + 
  geom_boxplot(aes(fill = factor(Estimator))) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "MAE", limits = c(0, 2.0)) +
  theme(legend.position = 'none')

ggplot(data = ESSs.melt, aes(x = factor(Sample_Size), y = Value)) + 
  geom_boxplot(aes(fill = factor(Estimator))) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "Log Effective Sample Size") + 
  theme(legend.title = element_blank())




