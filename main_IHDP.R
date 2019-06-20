
rm(list = ls())
library(ggplot2)
library(reshape2)
source("./estimators.R")
source("IHDPSimulator.R")

loadDataInCurrentEnvironment()
reps = 20                                 # no. of replications to be implemented for calculating confidence interval.

ATE_A = ATE_B = ATE_C = matrix(0, reps, 7, byrow = TRUE)
colnames(ATE_A) = colnames(ATE_B) = colnames(ATE_C) = c("True", "Naive", "IPW", "Trunc", "Norm", "PS", "PSNorm")
ATEs = MAEs = ESSs = list()

settings = list("A", "B", "C")
for (i in settings) {
  print(paste0("Starting Setting ", i))
  ate_i = data.frame(Setting=rep(i,reps), Truth=rep(0.,reps), Naive=rep(0.,reps), IPW=rep(0.,reps), Trunc=rep(0.,reps), Norm=rep(0.,reps), PS=rep(0.,reps), PSNorm=rep(0.,reps))
  mae_i = data.frame(Setting=rep(i,reps), Naive=rep(0.,reps), IPW=rep(0.,reps), Trunc=rep(0.,reps), Norm=rep(0.,reps), PS=rep(0.,reps), PSNorm=rep(0.,reps))
  ess_i = data.frame(Setting=rep(i,reps), Naive=rep(0.,reps), IPW=rep(0.,reps), Trunc=rep(0.,reps), Norm=rep(0.,reps), PS=rep(0.,reps), PSNorm=rep(0.,reps))
  
  for (k in 1:reps) {
    generateDataForIterInCurrentEnvironment(iter = 100*k, x = x, z = z, w = 0.5, overlap = TRUE, setting = i)
    data = data.frame(A=z.r, Y=y, Y1=y.1, Y0=y.0, u0=mu.0, u1=mu.1)
    
    ## estimate propensity scores -> IP weights
    getPropensityScoreInCurrentEnvironment(x = x.r, z = z.r, p.score = "logistic")
    # df = as.data.frame(x)
    # df$A = data$A
    # m <- glm(z ~ ., data = df, family = binomial())
    # fit = glm(A ~ ., family = "binomial", data = df)
    # ps.z  = fitted(fit)
    pAW = ps.z*as.numeric(data$A==1) + (1- ps.z)*as.numeric(data$A==0)
    wt  = 1/pAW
    
    ##  ate_true   = mean(data$u1 - data$u0)
    ate_true   = 4.0
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
    
    
    ate_i[k, ] = c(i, ate_true, ate_naive$ate, ate_ipw$ate, ate_trunc$ate, ate_norm$ate, ate_ps$ate, ate_psnorm$ate)
    mae_i[k, ] = c(i, mae_naive, mae_ipw, mae_trunc, mae_norm, mae_ps, mae_psNorm)
    ess_i[k, ] = c(i, ate_naive$N_eff, ate_ipw$N_eff, ate_trunc$N_eff, ate_norm$N_eff, ate_ps$N_eff, ate_psnorm$N_eff)
    ATEs[[i]] = ate_i
    MAEs[[i]] = mae_i
    ESSs[[i]] = ess_i
  }
  print(paste0("Experiments with setting = ", i, " Done!"))
}


#### Parse results & visualization. ------------
ate.table = data.frame(setting=c("A","B","C"), truth=rep(NA,3), naive=rep(NA,3),ipw=rep(NA,3),
                                          trunc=rep(NA,3), norm=rep(NA,3), ps=rep(NA,3),ps.norm=rep(NA,3))
for (S in 1:length(settings)){
  ates = colMeans(data.matrix(ATEs[[S]][,-1]))
  ate.table[S,2:8] = ates[1:7]
}
ate.table
#write.csv(x = ate.table, file = "./ateEstIHDP.csv")


for (i in 1:length(settings)){
  MAEs[[i]] = melt(MAEs[[i]], id.var = "Setting")
  ESSs[[i]] = melt(ESSs[[i]], id.var = "Setting")
}

MAEs.melt = MAEs[[1]]
ESSs.melt = ESSs[[1]]
for (i in 2:length(settings)){
  MAEs.melt = rbind(MAEs.melt, MAEs[[i]])
  ESSs.melt = rbind(ESSs.melt, ESSs[[i]])
}

ESSs.melt[,3] = log(as.numeric(ESSs.melt[,3]))
MAEs.melt[,3] = as.numeric(MAEs.melt[,3])
names(MAEs.melt) = c("Setting", "Estimator", "Value")
names(ESSs.melt) = c("Setting", "Estimator", "Value")


ggplot(data = MAEs.melt, aes(x = Setting, y = Value)) + 
  geom_boxplot(aes(fill = factor(Estimator))) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "MAE", limits = c(0,7)) +
  theme(legend.position = "none")

ggplot(data = ESSs.melt, aes(x = Setting, y = Value)) + 
  geom_boxplot(aes(fill = factor(Estimator))) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "Log Effective Sample Size") + 
  theme(legend.title = element_blank())


ATE.A = data.frame(data.matrix(ATEs[[1]][,-1]))
ATE.B = data.frame(data.matrix(ATEs[[2]][,-1]))
ATE.C = data.frame(data.matrix(ATEs[[3]][,-1]))
ATE.A = melt(ATE.A)
ATE.B = melt(ATE.B)
ATE.C = melt(ATE.C)


agg.A = aggregate.data.frame(ATE.A$value, by = list(Estimator=ATE.A$variable), 
                               FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
agg.A = do.call(data.frame, agg.A)
agg.A$se = agg.A$x.sd / sqrt(agg.A$x.n)
colnames(agg.A) = c("Estimator", "mean", "sd", "n", "se")

agg.B = aggregate.data.frame(ATE.B$value, by = list(Estimator=ATE.B$variable), 
                             FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
agg.B = do.call(data.frame, agg.B)
agg.B$se = agg.B$x.sd / sqrt(agg.B$x.n)
colnames(agg.B) = c("Estimator", "mean", "sd", "n", "se")

agg.C = aggregate.data.frame(ATE.C$value, by = list(Estimator=ATE.C$variable), 
                             FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
agg.C = do.call(data.frame, agg.C)
agg.C$se = agg.C$x.sd / sqrt(agg.C$x.n)
colnames(agg.C) = c("Estimator", "mean", "sd", "n", "se")



