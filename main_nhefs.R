
rm(list = ls())
library(ggplot2)
library(reshape2)
source("./estimators.R")

nhefs = as.data.frame(read.csv('./data/nhefs.csv'))       # 1746 x 61
nhefs$cens = as.numeric(is.na(nhefs$wt82))
nhefs$older = as.numeric(nhefs$age > 50 & !is.na(nhefs$age))
nhefs.original = nhefs

nhefs$id = 1:nrow(nhefs)
# group = c(-0.5,8.5, 11.5, 12.5, 15.5, 16.5)         # 0-8, 9-11, 12, 13-15, >=16
# education.code = cut(x = nhefs$school, breaks = group, labels = c(1,2,3,4,5))
# education = cut(x = nhefs$school, breaks = group, labels = c('1. 8th grade or less', '2. HS dropout',
#                                                              '3. HS', '4. College dropout',
#                                                              '5. College or more'))
library(car)
nhefs$education.code <- recode(nhefs$school, " 0:8 = 1; 9:11 = 2; 12 = 3; 13:15 = 4; 16:hi = 5; NA = 6 ")
nhefs$education <- recode(nhefs$education.code, " 1 = '1. 8th grade or less';
                          2 = '2. HS dropout'; 3 = '3. HS';
                          4 = '4. College dropout';  5 = '5. College or more'; 6 = 'Unknown' ")

## nhefs (1746 x 66),  nhefs (1746 x 12)
# nhefs = data.frame(nhefs, education.code, education)
nhefs2 <- nhefs[c("id", "qsmk", "sex", "race", "age", "school",
                  "smokeintensity", "smokeyrs", "exercise", "active", "wt71", "wt82")]
nhefs2 = as.data.frame(na.omit(nhefs2))   # 1566 x 12
nhefs = subset(nhefs, id %in% nhefs2$id)  # 1566 x 66
nhefs0 = subset(x = nhefs, cens==0)       # restricting data for uncensored 1566 x 66

u0 = mean(nhefs0$wt82_71[nhefs$qsmk==0])
se0 = sd(nhefs0$wt82_71[nhefs$qsmk==0]) / sqrt(sum(nhefs0$qsmk==0))
u1 = mean(nhefs0$wt82_71[nhefs$qsmk==1])
se1 = sd(nhefs0$wt82_71[nhefs$qsmk==1]) / sqrt(sum(nhefs0$qsmk==1))

colors = c(rgb(0,0,1, 0.5), rgb(1,0,0, 0.5))
wg = data.frame(group=c("Non-quitter","Quitter"), u=c(u0,u1), se=c(se0,se1))
ggplot(data = wg, aes(x = group, y = u)) + 
  geom_point(size=2.5) +
  geom_errorbar(aes(ymax=u+se, ymin=u-se, width=0.2)) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Weight gain (kg)") +
  theme(legend.position = 'none')


reps = 20
sampleRatio = c(0.5, 0.6, 0.7, 0.8, 0.9)
ratio = sampleRatio[3]
ATEs = ESSs = data.frame(Naive=rep(0,reps), IPW=rep(0,reps), Trunc=rep(0,reps), Norm=rep(0,reps), PS=rep(0,reps), PSNorm=rep(0,reps))

for (k in 1:reps){
  print(paste0("Replicaiton ", k))
  set.seed(100*k)
  # idx = sample(1:nrow(nhefs), size = 1400, replace = FALSE)
  # Sample randomly 80% treated & control group data
  size0 = floor(ratio * sum(nhefs$qsmk==0))                        # 930
  idx0 = sample(nhefs$id[nhefs$qsmk==0], size = size0, replace = FALSE)
  size1 = floor(ratio * sum(nhefs$qsmk==1))                        # 322
  idx1 = sample(nhefs$id[nhefs$qsmk==1], size = size1, replace = FALSE)
  data = data = subset(nhefs0, id %in% c(idx0, idx1))
  data$A = data$qsmk
  data$Y = data$wt82_71
  
  ## estimate propensity scores -> IP weights
  fit <- glm(A ~ as.factor(sex) + as.factor(race) + age + I(age^2) + 
               as.factor(education.code) + smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + 
               as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), 
             family = "binomial", data = data)

  ps  = predict(fit, type = "response")
  pAW = ps*data$A + (1-ps)*(1-data$A)
  wt  = 1/pAW
  
  ## 
  ate_naive  = naive(A = data$A, Y = data$Y)
  ate_ipw    = iptw(A = data$A, Y = data$Y, w = wt)
  ate_trunc  = iptw_truncate(A = data$A, Y = data$Y, w = wt, ceiling = 0)
  ate_norm   = iptw_norm(A = data$A, Y = data$Y, w = wt)
  ate_ps     = iptw_ps(A = data$A, Y = data$Y, w = wt, wcp = 0.2, min_cutoff = 1)
  ate_psnorm = iptw_ps_norm(A = data$A, Y = data$Y, w = wt, wcp = 0.2, min_cutoff = 1)
  
  ATEs[k, ] = c(ate_naive$ate, ate_ipw$ate, ate_trunc$ate, ate_norm$ate, ate_ps$ate, ate_psnorm$ate)
  ESSs[k, ] = c(ate_naive$N_eff, ate_ipw$N_eff, ate_trunc$N_eff, ate_norm$N_eff, ate_ps$N_eff, ate_psnorm$N_eff)
}


ATE.melt = melt(ATEs)
ESS.melt = melt(ESSs)
ESS.melt[,2] = log(ESS.melt[,2])
names(ATE.melt) = names(ESS.melt) = c("Estimator", "Value")


ate.agg = aggregate.data.frame(ATE.melt$Value, by = list(Estimator=ATE.melt$Estimator), 
                               FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
ate.agg = do.call(data.frame, ate.agg)
ate.agg$se = ate.agg$x.sd / sqrt(ate.agg$x.n)
colnames(ate.agg) = c("Estimator", "mean", "sd", "n", "se")

# ess.agg = aggregate.data.frame(ESS.melt$Value, by = list(Estimator=ESS.melt$Estimator), 
#                                FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
# ess.agg = do.call(data.frame, ess.agg)
# ess.agg$se = ess.agg$x.sd / sqrt(ess.agg$x.n)
# colnames(ess.agg) = c("Estimator", "mean", "sd", "n", "se")

ggplot(data = ate.agg, aes(x = Estimator, y = mean, fill=Estimator))+
  geom_bar(stat = "identity", position = position_dodge(0.9)) + 
  geom_errorbar(aes(ymax = mean+se, ymin = mean-se), width=0.5) + 
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "ATE Estimate") 


ggplot(data = ESS.melt, aes(x = Estimator, y = Value)) + 
  geom_boxplot(aes(fill = factor(Estimator))) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "Log Effective Sample Size") + 
  theme(legend.title = element_blank())

ggplot(data = ATE.melt, aes(x = Estimator, y = Value)) + 
  geom_boxplot(aes(fill = factor(Estimator))) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "ATE") + 
  theme(legend.position = 'none', legend.title = element_blank())

ggplot() + 
  geom_histogram(aes(x = wt, group = data$A, fill=factor(data$A))) +
  scale_x_continuous(name = "") + 
  scale_y_continuous(name = "ATE") +
  theme(legend.title = element_text(""))

ggplot() + 
  geom_histogram(aes(x = pAW, group = data$A, fill=factor(data$A)), breaks=seq(0.0,1.0, by = 0.05)) +
  scale_x_continuous(name = "") + 
  scale_y_continuous(name = "ATE") +
  theme(legend.title = element_text(""))


