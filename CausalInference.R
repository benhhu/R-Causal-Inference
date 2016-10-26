#Causal Inference
#install.packages("geepack")
#install.packages("multcomp")
#install.packages("BSagri")
#install.packages("sandwich")
#install.packages("descr")
#install.packages("boot")
library(sas7bdat)
library(car)
library(geepack)
library(multcomp)
library(BSagri)
library(sandwich)
library(descr)
library(boot)

setwd("//file.phhp.ufl.edu/home/huihu/files/RWD/R-Causal-Inference")
nhefs <- read.sas7bdat("nhefs_book.sas7bdat")
nhefs$cens <- as.numeric(is.na(nhefs$wt82))
nhefs$older <- as.numeric(nhefs$age > 50 & !is.na(nhefs$age))

nhefs$education.code <- recode(nhefs$school, " 0:8 = 1; 9:11 = 2; 12 = 3; 13:15 = 4; 16:hi = 5; NA = 6 ")
nhefs$education <- recode(nhefs$education.code, " 1 = '1. 8th grade or less'; 2 = '2. HS dropout'; 3 = '3. HS'; 4 = '4. College dropout'; 5 = '5. College or more'; 6 = 'Unknown' ")

# Analysis restricted to N=1566 
# with non-missing values in the following covariates
nhefs.original <- nhefs # Original data saved for later use
nhefs$id <- 1:nrow(nhefs)
nhefs2 <- nhefs[c("id", "qsmk", "sex", "race", "age", "school","smokeintensity", "smokeyrs", "exercise", "active","wt71", "wt82")]
dim(nhefs2)
nhefs2 <- as.data.frame(na.omit(nhefs2))
dim(nhefs2)
nhefs <- subset(nhefs, id %in% nhefs2$id)
rm(nhefs2)
dim(nhefs)

# mean weight change in those with and without smoking cessation
summary(nhefs$wt82_71[nhefs$cens == 0 & nhefs$qsmk == 0])
sd(nhefs$wt82_71[nhefs$cens == 0 & nhefs$qsmk == 0])
summary(nhefs$wt82_71[nhefs$cens == 0 & nhefs$qsmk == 1])
sd(nhefs$wt82_71[nhefs$cens == 0 & nhefs$qsmk == 1])

summary(glm(wt82_71 ~ qsmk, data = subset(nhefs, cens == 0)))

# restricting data for uncensored
nhefs0 <- subset(nhefs, cens == 0)

# baseline characteristics
years1 <-mean(nhefs0$age[nhefs0$qsmk == 1])
years0 <-mean(nhefs0$age[nhefs0$qsmk == 0])
male1 <-100*mean(nhefs0$sex[nhefs0$qsmk == 1]==0)
male0 <-100*mean(nhefs0$sex[nhefs0$qsmk == 0]==0)
white1 <-100*mean(nhefs0$race[nhefs0$qsmk == 1]==0)
white0 <-100*mean(nhefs0$race[nhefs0$qsmk == 0]==0)
university1 <-100*mean(nhefs0$education.code[nhefs0$qsmk == 1]==5)
university0 <-100*mean(nhefs0$education.code[nhefs0$qsmk == 0]==5)
kg1 <-mean(nhefs0$wt71[nhefs0$qsmk == 1])
kg0 <-mean(nhefs0$wt71[nhefs0$qsmk == 0])
cigs1 <-mean(nhefs0$smokeintensity[nhefs0$qsmk == 1])
cigs0 <-mean(nhefs0$smokeintensity[nhefs0$qsmk == 0])
smoke1 <-mean(nhefs0$smokeyrs[nhefs0$qsmk == 1])
smoke0 <-mean(nhefs0$smokeyrs[nhefs0$qsmk == 0])
noexer1 <-100*mean(nhefs0$exercise[nhefs0$qsmk == 1]==2)
noexer0 <-100*mean(nhefs0$exercise[nhefs0$qsmk == 0]==2)
inactive1 <-100*mean(nhefs0$active[nhefs0$qsmk == 1]==2)
inactive0 <-100*mean(nhefs0$active[nhefs0$qsmk == 0]==2)

baseline.char <- round(matrix(c(years1, years0, male1, male0, white1, white0, university1, university0, kg1, kg0, cigs1, cigs0, smoke1, smoke0, noexer1, noexer0, inactive1, inactive0), ncol = 2, byrow=T), 1) 
dimnames(baseline.char) <- list(c("age, years", "men, %", "white, %", "university, %", "weight, kg", "Chigarettes/day", "year smoking", "little/no exercise, %", "inactive daily life, %"), c("Smoking cessation (A=1)","No smoking cessation (A=0)"))
print(baseline.char)

##############
#IP weighting#
##############

# Estimation of ip weights via a logistic model
fit <- glm(qsmk ~ as.factor(sex) + as.factor(race) + age + I(age^2) + as.factor(education.code) + smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), family = binomial(), data = nhefs0)
summary(fit)

p.qsmk.obs <- ifelse(nhefs0$qsmk == 0, 1 - predict(fit, type = "response"), predict(fit, type = "response"))
nhefs0$w <- 1/p.qsmk.obs
summary(nhefs0$w)
sd(nhefs0$w)

# Estimates from a GLM with cluster option
glm.obj <- glm(wt82_71 ~ qsmk + cluster(seqn), data = nhefs0, weights = w)
beta <- coef(glm.obj)
SE <-sqrt(diag(vcovHC(glm.obj, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
round(cbind(beta, lcl, ucl),1)[2,]

# no association between sex and qsmk in pseudo-population
crosstab(nhefs0$sex, nhefs0$qsmk, weight = nhefs0$w, plot = F, format = "SAS")
# Numbers are approximately double: weight = nhefs0$w/2 makes it close

# "check" for positivity
age.restr <- nhefs0$age[nhefs0$race == 0 & nhefs0$sex == 1]
qsmk.restr <- nhefs0$qsmk[nhefs0$race == 0 & nhefs0$sex == 1]
CrossTable(age.restr, qsmk.restr, expected = F, format = "SAS")

#Stabilized IP Weights
# estimation of denominator of ip weights
denom.fit <- glm(qsmk ~ as.factor(sex) + as.factor(race) + age + I(age^2) + as.factor(education.code) + smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), family = binomial(), data = nhefs0)
denom.p <- predict(denom.fit, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(qsmk~1, family = binomial(), data = nhefs0)
numer.p <- predict(numer.fit, type = "response")

nhefs0$sw <- ifelse(nhefs0$qsmk == 0, ((1-numer.p)/(1-denom.p)), (numer.p/denom.p))

summary(nhefs0$sw)

glm.obj <- glm(wt82_71~ qsmk + cluster(seqn), data = nhefs0, weights = sw)
beta <- coef(glm.obj)
SE <-sqrt(diag(vcovHC(glm.obj, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
round(cbind(beta, lcl, ucl),1)[2,]

# no association between sex and qsmk in pseudo-population
crosstab(nhefs0$sex, nhefs0$qsmk, weight = nhefs0$sw, plot = F, format = "SAS")

#Marginal structural model with a continouous treatment
# Analysis restricted to subjects reporting <=25 cig/day at baseline
nhefs1 <- subset(nhefs0, smokeintensity <=25)
dim(nhefs1)

# estimation of denominator of ip weights
den.fit.obj <- lm(smkintensity82_71 ~ as.factor(sex) + as.factor(race) + age + I(age^2) + as.factor(education.code) + smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71  + I(wt71^2), data = nhefs1)
p.den <- predict(den.fit.obj, type = "response")
dens.den <- dnorm(nhefs1$smkintensity82_71, p.den, summary(den.fit.obj)$sigma)

# estimation of numerator of ip weights
num.fit.obj <- lm(smkintensity82_71 ~ 1, data = nhefs1)
p.num <- predict(num.fit.obj, type = "response")
dens.num <- dnorm(nhefs1$smkintensity82_71, p.num, summary(num.fit.obj)$sigma)

# estimation of Stabilized weights
nhefs1$sw.a = dens.num/dens.den
summary(nhefs1$sw.a)

gee.obj <- geeglm(wt82_71~smkintensity82_71 + I(smkintensity82_71^2), data = nhefs1, std.err = 'san.se', weights = sw.a, id=seqn, corstr="independence")
comp<-glht(gee.obj)
CIadj<-CIGLM(comp,method="Adj")
CIadj

#Marginal structural logistic model
# Estimating the parameters of a marginal structural logistic model
glm.obj <- glm(death ~ qsmk + cluster(seqn), data = nhefs0, weights = sw, family = binomial())
beta <- coef(glm.obj)
SE <-sqrt(diag(vcovHC(glm.obj, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
round(exp(cbind(beta, lcl, ucl)),1)[2,]

#Adjust for censoring
# using the dataset nhefs.original from PROGRAM 12.1
nhefs.original$id <- 1:nrow(nhefs.original)
nhefsx <- nhefs.original[c("id", "qsmk", "sex", "race", "age", "school","smokeintensity", "smokeyrs", "exercise", "active", "wt71")]
# Missing values of wt82 are retained while the others are dropped
dim(nhefsx)
nhefsx <- as.data.frame(na.omit(nhefsx))
dim(nhefsx)
nhefs3 <- subset(nhefs.original, id %in% nhefsx$id)
dim(nhefs3)
# Analysis restricted to N=1629. 
crosstab(nhefs3$cens, nhefs3$qsmk, plot = F, format = "SAS")
# estimation of denominator of ip weights for treatment
denom.fit <- glm(qsmk ~ as.factor(sex) + as.factor(race) + age + I(age^2) + as.factor(education.code) + smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), family = binomial(), data = nhefs3)
denom.p <- predict(denom.fit, type = "response")
# estimation of numerator of ip weights for treatment
numer.fit <- glm(qsmk~1, family = binomial(), data = nhefs3)
numer.p <- predict(numer.fit, type = "response")
# Estimation of stabilized weight for treatment
nhefs3$sw.a <- ifelse(nhefs3$qsmk == 0, ((1-numer.p)/(1-denom.p)),(numer.p/denom.p))
summary(nhefs3$sw.a)
# estimation of denominator of ip weights for not being censored
denom.cens <- glm(cens ~ as.factor(qsmk) + as.factor(sex) + as.factor(race) + age + I(age^2) + as.factor(education.code) + smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), family = binomial(), data = nhefs3)
denom.p.cens <- predict(denom.cens, type = "response")
# estimation of numerator of ip weights for not being censored
numer.cens <- glm(cens~as.factor(qsmk), family = binomial(), data = nhefs3)
numer.p.cens <- predict(numer.cens, type = "response")
# Estimation of stabilized weight for not being censored
nhefs3$sw.c <- ifelse(nhefs3$cens == 0, ((1-numer.p.cens)/(1-denom.p.cens)), 1)
summary(nhefs3$sw.c)
sd(nhefs3$sw.c)
# Estimation of Stabilized Censoring weight (sw)
nhefs3$sw <- nhefs3$sw.a * nhefs3$sw.c
summary(nhefs3$sw)
sd(nhefs3$sw)
# obtaining final estimates
glm.obj <- glm(wt82_71~as.factor(qsmk) + cluster(seqn), data = nhefs3,weights = sw)
beta <- coef(glm.obj)
SE <-sqrt(diag(vcovHC(glm.obj, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
p.value <- 2*(1-pnorm(abs(beta/SE)))
round(cbind(beta, SE, lcl, ucl, p.value),1)[2,]

###############################
#Standardization and G-Formula#
###############################

#Estimating the Mean Outcome via Modeling
# Analysis restricted to N=1566 
# with non-missing values in the following covariates
nhefs.original <- nhefs # Original data saved for later use
nhefs$id <- 1:nrow(nhefs)
nhefs2 <- nhefs[c("id", "qsmk", "sex", "race", "age", "school", "smokeintensity", "smokeyrs", "exercise", "active", "wt71", "wt82")]
dim(nhefs2)
# restricting data for non-missing
nhefs2 <- as.data.frame(na.omit(nhefs2))
dim(nhefs2)
nhefs.id.matched <- subset(nhefs, id %in% nhefs2$id)
dim(nhefs)
# restricting data for uncensored, for comparison with observed outcome
nhefs0 <- subset(nhefs.id.matched, cens == 0)
dim(nhefs0)
table(nhefs0$qsmk) # untreated vs treated
# Estimates
glm.obj <- glm(wt82_71~ as.factor(qsmk) + as.factor(sex) + as.factor(race) + age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2) + I(qsmk*smokeintensity), data = nhefs0)
nhefs0$meanY <- predict(glm.obj, type = "response")
# print variable values corresponding to subject with unique identifier 24770
nhefs0[nhefs0$seqn ==24770,c("meanY", "qsmk", "sex", "race", "age", "education","smokeintensity", "smokeyrs", "exercise", "active","wt71")]
summary(nhefs0$meanY)
summary(nhefs0$wt82_71)

#Standardizing the Mean Outcome to the Confounder Distribution
# 1st copy: equal to original one
nhefs0$interv <- rep(-1, nrow(nhefs0))
dim(nhefs0)
# 2nd copy: treatment set to 0, outcome to missing
nhefs.untr <- nhefs0
nhefs.untr$interv <- rep(0, nrow(nhefs.untr))
nhefs.untr$qsmk <- rep(0, nrow(nhefs.untr))
nhefs.untr$wt82_71 <- rep(NA, nrow(nhefs.untr))
dim(nhefs.untr)
# 3rd copy: treatment set to 1, outcome to missing
nhefs.tr <- nhefs0
nhefs.tr$interv <- rep(1, nrow(nhefs.tr))
nhefs.tr$qsmk <- rep(1, nrow(nhefs.tr))
nhefs.tr$wt82_71 <- rep(NA, nrow(nhefs.tr))
dim(nhefs.tr)
# create a dataset with 3 copies of each subject
onesample <- as.data.frame(rbind(nhefs0, nhefs.untr, nhefs.tr))
dim(onesample)
# Estimates
# linear model to estimate mean outcome conditional on treatment & confounders,
# parameters are estimated using original observations only (interv= -1),
# parameter estimates are used to predict mean outcome for observations 
# with treatment set to 0 (interv=0) and to 1 (innterv=1);
glm.obj <- glm(wt82_71~ as.factor(qsmk) + as.factor(sex) + as.factor(race) + age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), data = onesample)
onesample$meanY <- predict(glm.obj, onesample, type = "response")
# estimate mean outcome in each of the groups interv=0, and interv=1;
# this mean outcome is a weighted average of the mean outcomes 
# in each combination of values of treatment and confounders, 
# that is, the standardized outcome;
with(onesample, tapply(meanY, list(interv), mean))

#Use Bootstrap to Compute the 95% CI
# Compute basic bootstrap confidence interval.
# Below is the function to estimate mean outcome in each of the groups 
# interv=-1, interv=0 & interv=1 conditional on treatment and confounders,
# also include the mean outcome of interv = 1 vs. interv = 0
boot.fun <- function(dat, index, type){
  sampled.data <- dat[index,]
  fit <- glm(formula = wt82_71~ as.factor(qsmk) + as.factor(sex) + as.factor(race) + 
               age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + 
               smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + 
               I(wt71^2), data = sampled.data)
  sampled.data$meanY <- predict(fit, sampled.data, type = "response")
  std.mean <- with(sampled.data, tapply(meanY, list(interv), mean))
  mean.list <- as.numeric(unlist(c(std.mean, diff(std.mean[c(2,3)]))))
  return(mean.list)
}
# Parametric g-formula: Bootstrap results using nboot samples
nboot <- 100
# number of bootstrap samples to run
# creation of bootstarp samples using fixed seed to reproduce results
set.seed(1232)
bootres <- boot(onesample, boot.fun, R = nboot)
bootres
nrow(bootres$t)
boot.stat <- matrix(NA, 4, 4)
for (i in 1:4){
  boot.stat[i,] <- c(bootres$t0[i], sd(bootres$t[,i]), 
                     bootres$t0[i] - 1.96 * sd(bootres$t[,i]), 
                     bootres$t0[i] + 1.96 * sd(bootres$t[,i]))
}
dimnames(boot.stat) <- list(c("observed", "No treatment", 
                              "Treatment", "Difference"), 
                            c("mean", "sd", "LCL", "UCL"))
boot.stat




