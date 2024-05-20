#################################################################################
#  Analysis code for the paper "Quantifying social inequalities in eligibility  #
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
#  Analysis: Plan size and Spending Inequality                                  #
#  Using women/girls (Yes,No) and psychosocial disability as an example         #
#                                                                               #
#  Author: Yi Yang                                                              #
#  Date: 20 May, 2024                                                           #
#################################################################################

#--------------------------------- Function ------------------------------------#  

### Arguments for function:
# dat: A data.frame with the data for analysis
# ind: A vector of indices that define the sample from dat on which to conduct the analysis. Defaults to all rows of dat.
#      This facilitates the use of this function within the boot() function from the boot package
# exposure: inequality group, which must be a binary 0/1 variable (e.g., women or girls)
# outcome: access request outcome, which must be a binary 0/1 variable, with 1=eligible and 0=non-eligible
# confounders: List of confounder(s)
# confounders_factor: confounder(s) and interation term(s) to be included in the model
# mcsim: Number of Monte Carlo simulations to conduct


function_plan_spend <-function(dat, ind=1:nrow(dat),exposure, outcome, mediators, confounders, confounders_factor,mcsim=200)
{

  #Rename all variables & prepare dataset
  dat$A    <-dat[,exposure]
  dat$Y    <-dat[,outcome]
  dat$M1   <-dat[,mediators[1]]
  dat$M1log<-log(dat$M1)
  
  dat<-dat[,c("A","M1","M1log","Y",confounders)]
  
  #Take boostrap sample
  data<-dat[ind,]
  
  #Set flag to capure bootstrap samples to reject
  flag<-FALSE 
  
  #Prepare confounders for formulae
  confs<-paste(confounders_factor,collapse="+")
  
  #Replicate dataset for simulations
  dat2<-data
  dat2[,1:4]<-NA_integer_
  dat2<-coredata(dat2)[rep(seq(nrow(dat2)),mcsim),]
  n<-nrow(dat2)
  
  
  ## ESTIMATE DISTRIBUTIONS ##
  
  # Plan size
  
  fit<-glm(as.formula(paste("M1~A+",confs)),data=data,family=Gamma(link="log"))
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  disp<-sum(resid(fit, type='pear')^2)/fit$df.residual
  
  dat2$A<-0
  mean<- predict(fit,newdata=dat2,type="response")
  var <- mean^2*disp
  sha <- mean ^2 / var   
  sca <- var / mean 
  m1_0_mmmm <- rgamma(n, shape = sha, scale = sca)
  m1_0_mmmm_log <- log(m1_0_mmmm)
  
  dat2$A<-1
  mean <- predict(fit,newdata=dat2,type="response")
  var  <- mean^2*disp
  sha  <- mean ^2 / var  
  sca  <- var / mean 
  m1_1_mmmm <- rgamma(n, shape = sha, scale = sca)
  m1_1_mmmm_log <- log(m1_1_mmmm)     
  
  ### Spending
  
  fit<-glm(as.formula(paste("Y~ A*M1log+",confs)),data=data,family=Gamma(link="log"))
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  disp<-sum(resid(fit, type='pear')^2)/fit$df.residual
  
  dat2$A<-0
  dat2$M1log<-m1_0_mmmm_log
  mean <- predict(fit,newdata=dat2,type="response")
  var <- mean^2*disp
  sha  <- mean ^2 / var   
  sca  <- var / mean 
  y0 <- rgamma(n, shape = sha, scale = sca)
  
  
  dat2$A<-1
  dat2$M1log<-m1_1_mmmm_log
  mean <- predict(fit,newdata=dat2,type="response")
  var <- mean^2*disp
  sha  <- mean ^2 / var  
  sca  <- var / mean 
  y1 <- rgamma(n, shape = sha, scale = sca)
  
  
  dat2$A<-1
  dat2$M1log<-m1_0_mmmm_log
  mean <- predict(fit,newdata=dat2,type="response")
  var <- mean^2*disp
  sha  <- mean ^2 / var  
  sca  <- var / mean 
  y1_m0 <- rgamma(n, shape = sha, scale = sca)
  
  #*#####################*#
  M0<-mean(m1_0_mmmm)
  M1<-mean(m1_1_mmmm)
  Y0<-mean(y0)
  Y1<-mean(y1)
  Y1_m0 <- mean(y1_m0)
  
  TCE_plan <- M1 - M0
  TCE_paid <- Y1 - Y0

  # IIE
  IIE  <-Y1-Y1_m0
  
  # IDE
  IDE <- Y1_m0-Y0
  
  ### Collect and return results
  res<-c(M0,M1,Y0,Y1,Y1_m0,
         TCE_plan,TCE_paid,IIE,IDE) 
}


#-------------------  Plan size and spending inequalities ----------------------#  

####### Using women/girls (Yes,No) and psychosocial disability as an example 
library(dplyr)
library(boot)

#load data
dat <- NDIS_data %>% filter(disability == "Psychosocial Disability")

# Estimand
estimand <- c("M0","M1","Y0","Y1","Y1_m0","TCE_plan","TCE_paid","IIE","IDE")


# Confounder(s) and interaction term(s)
# --- list of confounders
confounders <- c("First_Nations","Age_PlanStart","Not_Major_City","SES_Deciles",
                 "Severity_Score","Other_Supports","previous_NDIS_years","PreNDIS_support") 

# --- confounders and interaction terms
confounders_factor <-c("First_Nations","factor(Age_PlanStart)","Not_Major_City","SES_Deciles",
                       "factor(Severity_Score)","Other_Supports","factor(previous_NDIS_years)","factor(PreNDIS_support)",
                       "A:Severity_Score","A:previous_NDIS_years") 

set.seed(2021) 

bstrap <-boot(data = dat, statistic = function_plan_spend, 
              exposure    = "womengirls", 
              outcome     = "annualised_total_spend", 
              mediators   = "annualised_total_plan",
              confounders = confounders, 
              confounders_factor=confounders_factor,
              mcsim=10, R=500)  # Number of Monte Carlo simulations and bootstrap samples

# Set-up results table
RES<-data.frame(Estimate=bstrap$t0,SE=apply(bstrap$t,2,sd,na.rm=T))
RES$CIlow<-RES$Estimate-1.96*RES$SE
RES$CIupp<-RES$Estimate+1.96*RES$SE

RES_womengirls_plan_spend <- cbind(Estimand=estimand,RES)
