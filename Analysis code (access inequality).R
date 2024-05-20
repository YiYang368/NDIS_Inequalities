#################################################################################
#  Analysis code for the paper "Quantifying social inequalities in eligibility   #
#  for, and use of, the Australian National Disability Insurance Scheme"        #
#                                                                               #
#  The code was adapted from the medRCT.R function published by                 #
#  Margarita Moreno-Betancur from the following paper                           # 
#  (https://github.com/moreno-betancur/medRCT/blob/master/medRCT.R):            #
#  Moreno-Betancur et al. "Mediation effects that emulate a target randomised   #
#  trial: Simulation-based evaluation of ill-defined interventions on multiple  #
#  mediators" Statistical Methods in Medical Research 2021. 30(6):1395-1412     #
#                                                                               #
#                                                                               #
#  Analysis: Access Inequality                                                  #
#  Using women/girls (Yes,No) and psychosocial disability as an example         #
#                                                                               #
#  Author: Yi Yang                                                              #
#  Date: 20 May, 2024                                                           #
#################################################################################

#--------------------------------- Function ------------------------------------#  

### Arguments for function:
# dat: A data.frame with the data for analysis
# ind: A vector of indices that define the sample from dat on which to conduct the analysis. Defaults to all rows of dat.
# exposure: inequality group, which must be a binary 0/1 variable (e.g., women or girls)
# outcome: access request outcome, which must be a binary 0/1 variable, with 1=eligible and 0=non-eligible
# confounders: List of confounder(s)
# confounders_factor: confounder(s) and interation term(s) to be included in the model
# mcsim: Number of Monte Carlo simulations to conduct

function_access <-function(dat, ind=1:nrow(dat),exposure, outcome, confounders,confounders_factor, mcsim=200)
{
  #Rename all variables & prepare dataset
  dat$A<-dat[,exposure]
  dat$Y<-dat[,outcome]
  dat<-dat[,c("A","Y",confounders)]

  #Take boostrap sample
  data<-dat[ind,]
  
  #Set flag to capure bootstrap samples to reject
  flag<-FALSE 
  
  #Prepare confounders for formulae
  confs<-paste(confounders_factor,collapse="+")
  
  #Replicate dataset for simulations
  dat2<-data
  dat2[,1:2]<-NA_integer_
  dat2<-coredata(dat2)[rep(seq(nrow(dat2)),mcsim),]
  n<-nrow(dat2)
  
  ## ESTIMATE DISTRIBUTIONS ##
  # Y
  fit<-glm(as.formula(paste("Y~ A+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  ### ESTIMATE OUTCOME EXPECTATION IN EACH ARM ###
  
  #comparator group
  
  dat2$A<-0
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  Y0<-mean(y0)
  
  #inequality group
  
  dat2$A<-1
  y1<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  Y1<-mean(y1)
  
  ### ESTIMATE EFFECTS ###
  Ydiff <- Y1-Y0
  
  res<-c(Y0,Y1,Ydiff)
}

#----------------------------  Access inequality ------------------------------#  

####### Using women/girls (Yes,No) and psychosocial disability as an example 

# load libraries
library(dplyr)
library(boot)

#load data
dat <- NDIS_data %>% filter(disability == "Psychosocial Disability")

# estimand
estimand <- c("Y0","Y1","Ydiff")

# Confounder(s) and interaction term(s)
confounders <-c("PreNDIS_support")
confounders_factor<-c("factor(PreNDIS_support)","A:factor(PreNDIS_support)")

# set seed
set.seed(4082)

# Estimate the effects with the boostrap
bstrap<-boot(data=dat, statistic=function_access,
             exposure="womengirls", outcome="eligibility",
             confounders=confounders, confounders_factor=confounders_factor, 
             mcsim=10, R=500)  # Number of Monte Carlo simulations and bootstrap samples

# Set-up results table
RES<-data.frame(Estimate=bstrap$t0,SE=apply(bstrap$t,2,sd,na.rm=T))
RES$CIlow<-RES$Estimate-1.96*RES$SE
RES$CIupp<-RES$Estimate+1.96*RES$SE

RES_womengirls_access <- cbind(Estimand=estimand,RES)
