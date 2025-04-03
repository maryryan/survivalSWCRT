#### LIBRARIES ####
library(tidyverse)
library(pracma)
library(survival)
library(stabledist)
library(latex2exp)

#### LOAD KENDALL'S TAU VARIATION SCENARIOS FOR SENSITIVITY ANALYESES (Figure 4; Appendix F) ####
tau_scenarios <- read.table("catheter-scenarios.txt", sep=",")

#### LOAD FUNCTIONS ####
source("functions-biometrics.R")

#### SAMPLE SIZE CALCULATION: 80% power, beta=0.4, m=35 (Section 5) ####
## Step 1: set design parameters ##
m <- 35 # cluster size
J <- 6 # number of time periods
tau_w <- 0.1 # within-period Kendall's tau
tau_b <- 0.05 #between-period Kendall's tau

# set null and alternative effect sizes #
beta0 <- 0
betaA <- 0.4

## Step 2: specify censoring rates and baseline hazard function ##
# set maximum censoring time #
Cp <- 1

# set administrative censoring rate #
pa <- 0.05

# calculate base baseline hazard that results in above administrative censoring rate #
lambdaDET <- function(Cp, pa){#kappa, pa){
  lambda0 <- (1/Cp)*(-log(pa))#^(1/kappa)
  return(lambda0)
}
lambda0_base <- lambdaDET(Cp, pa)

# specify how does the baseline hazard before over time? #
baseline_constant <- 0.05 # specify whether/how baseline hazard changes with time
lambda0 <- lambda0_base+baseline_constant*seq(0,(J-1)) # final baseline hazard

# calculate probability of a cluster being on treatment #
pi_b <- c(0, rep((1/(J-1)), (J-1)))

# set type I error and power #
alpha <- 0.05
power <- 0.8

z_typeI <- qnorm(1-alpha/2)
z_power <- qnorm(power)

## Step 3: Estimate variances and calculate sample size #
# estimate wald variance under the alternative #
design_varA <- sandwich_var(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta=betaA)

# calculate number of clusters needed for wald test #
n_wald <- ( design_varA$var_beta_sqrt * (z_typeI + z_power) )^2/betaA^2

# estimate score variances under null and alternative hypotheses #
design_var0_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta0=beta0, betaA=beta0)
design_varA_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta0=beta0,betaA=betaA)

# calculate number of clusters needed for Self & Mauritsen and Tang score tests #
n_SM <- ( design_varA_score$B * (z_typeI + z_power) )^2/abs(design_varA_score$score)^2
n_Tang <- ( design_varA_score$B * (z_power + z_typeI*(design_var0_score$B/design_varA_score$B)) )^2/abs(design_varA_score$score)^2

#### CATHTAG SENSITIVITY ANALYSIS: n=18, beta=0.4, increasing hazard (Section 5 - Figure 4) #### 
## Step 1: set design parameters ##
m <- 35 # cluster size
n <- 18 # number of clusters
J <- 6 # number of time periods
tau_w <- as.numeric(tau_scenarios[,1]) # within-period Kendall's tau
tau_b <- as.numeric(tau_scenarios[,2]) # between-period Kendall's tau

# set null and alternative effect sizes #
beta0 <- 0
betaA <- 0.4

## Step 2: specify censoring rates and baseline hazard function ##
# set maximum censoring time #
Cp <- 1

# set administrative censoring rate #
pa <- 0.05

# calculate base baseline hazard that results in above administrative censoring rate #
lambdaDET <- function(Cp, pa){
  lambda0 <- (1/Cp)*(-log(pa))
  return(lambda0)
}
lambda0_base <- lambdaDET(Cp, pa)

# specify how does the baseline hazard before over time? #
baseline_constant <- 0.05 # specify whether/how baseline hazard changes with time (increasing by 5% each period here)
lambda0 <- lambda0_base+baseline_constant*seq(0,(J-1)) # final baseline hazard

# calculate probability of a cluster being on treatment #
pi_b <- c(0, rep((1/(J-1)), (J-1)))


## Step 3: loop through Kendall's tau scenarios ##
sensitivity_results <- NULL

for(j in seq(length(tau_w))){
  
  ## Step 3.1: power for Wald t-test ##
  # estimate wald variance under the alternative #
  design_varA <- sandwich_var(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta=betaA)
  
  # calculate predicted power based on your sandwich variance and t test #
  design_power_t <- pt((abs(betaA)/(design_varA$var_beta_sqrt/sqrt(n)) - qt(0.975, n-1)), n-1)
  
  ## Step 3.2: power for score tests ##
  # estimate score variances under null and alternative hypotheses #
  design_var0_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta0=beta0, betaA=beta0)
  design_varA_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta0=beta0,betaA=betaA)
  
  # calculate predicted power based on Tang score test #
  score_power_predict_tang <-  pnorm(abs(design_varA_score$score)/(design_varA_score$B/sqrt(n)) - qnorm(0.975)*(design_var0_score$B/design_varA_score$B))
  # calculate predicted power based on S&M score test #
  score_power_predict_A <-  pnorm(abs(design_varA_score$score)/(design_varA_score$B/sqrt(n)) - qnorm(0.975))
  
  ## Step 3.3: store results ##
  sensitivity_results_temp <- cbind(n, m, J, baseline_constant, tau_w[j], tau_b[j], tau_b[j]/tau_w[j], design_power_t, score_power_predict_A, score_power_predict_tang)
  sensitivity_results <- rbind(sensitivity_results, sensitivity_results_temp)
  
  # progress message #
  print(paste(j, " out of ", length(tau_w), " scenarios done"))
}

colnames(sensitivity_results) <- c("n","m","J","baseline_constant","tau_w", "tau_b", "tau_ratio", "T", "S&M Score", "Tang Score")

## Step 4: plot sensitivity results (Figure 4) ##
sensitivity_colors <- c("#ffffe5", "#f7fcb9", "#d9f0a3","#addd8e","#78c679",
                        "#41ab5d","#238443","#006837","#004529", "#001e12")

sensitivity_results %>%
  as.data.frame() %>%
  pivot_longer(c("T", "S&M Score", "Tang Score"),
               names_to="paradigm", values_to = "power_predict") %>%
  mutate(paradigm=factor(paradigm, levels=c("T", "S&M Score", "Tang Score")),
         power_discrete = factor(case_when(power_predict <= 0.7 ~ "(0,0.7]",
                                           power_predict > 0.7 & power_predict <=0.75 ~ "(0.7,0.75]",
                                           power_predict > 0.75 & power_predict <=0.8 ~ "(0.75,0.8]",
                                           power_predict > 0.8 & power_predict <=0.85 ~ "(0.8,0.85]",
                                           power_predict > 0.85 & power_predict <=0.9 ~ "(0.85,0.9]",
                                           power_predict > 0.9 & power_predict <=0.92 ~ "(0.9,0.92]",
                                           power_predict > 0.92 & power_predict <=0.94 ~ "(0.92,0.94]",
                                           power_predict > 0.94 & power_predict <=0.96 ~ "(0.94,0.96]",
                                           power_predict > 0.96 & power_predict <=0.98 ~ "(0.96,0.98]",
                                           power_predict > 0.98 & power_predict <=1 ~ "(0.98,1]"))) %>%
  ggplot(aes(tau_ratio, tau_w, fill=power_discrete))+
  geom_tile(width=0.05)+
  facet_grid(vars(paradigm))+
  labs(x=TeX("$\\tau_b/\\tau_w$"),
       y=TeX("$\\tau_w$"),
       fill="Predicted Power")+
  scale_fill_manual(values=sensitivity_colors)+
  theme_minimal()+
  theme(strip.background =element_rect(fill="lightgray"))

#### CATHTAG SENSITIVITY ANALYSIS: n=24, beta=0.4, constant hazard (Appendix F - Figure F.1) #### 
## Step 1: set design parameters ##
m <- 35 # cluster size
n <- 24 # number of clusters
J <- 6 # number of time periods
tau_w <- as.numeric(tau_scenarios[,1]) # within-period Kendall's tau
tau_b <- as.numeric(tai_scenarios[,2]) # between-period Kendall's tau

# set null and alternative effect sizes #
beta0 <- 0
betaA <- 0.4

## Step 2: specify censoring rates and baseline hazard function ##
# set maximum censoring time #
Cp <- 1

# calculate base baseline hazard that results in above administrative censoring rate #
lambdaDET <- function(Cp, pa){
  lambda0 <- (1/Cp)*(-log(pa))
  return(lambda0)
}
lambda0_base <- lambdaDET(Cp, pa)

# specify how does the baseline hazard before over time? #
baseline_constant <- 0 # specify whether/how baseline hazard changes with time (no change over time here)
lambda0 <- lambda0_base+baseline_constant*seq(0,(J-1)) # final baseline hazard

# calculate probability of a cluster being on treatment #
pi_b <- c(0, rep((1/(J-1)), (J-1)))

## Step 3: loop through Kendall's tau scenarios ##
sensitivity_results <- matrix(NA,ncol=6, nrow=length(tau_w))

for(j in seq(length(tau_w))){
  
  ## Step 3.1: power for Wald t-test ##
  # estimate wald variance under the alternative #
  design_varA <- sandwich_var(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta=betaA)
  
  # calculate predicted power based on your sandwich variance and t test #
  design_power_t <- pt((abs(betaA)/(design_varA$var_beta_sqrt/sqrt(n)) - qt(0.975, n-1)), n-1)
  
  ## Step 3.2: power for score tests ##
  # estimate score variances under null and alternative hypotheses #
  design_var0_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta0=beta0, betaA=beta0)
  design_varA_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta0=beta0,betaA=betaA)
  
  # calculate predicted power based on Tang score test #
  score_power_predict_tang <-  pnorm(abs(design_varA_score$score)/(design_varA_score$B/sqrt(n)) - qnorm(0.975)*(design_var0_score$B/design_varA_score$B))
  # calculate predicted power based on S&M score test #
  score_power_predict_A <-  pnorm(abs(design_varA_score$score)/(design_varA_score$B/sqrt(n)) - qnorm(0.975))
  
  ## Step 3.3: store results ##
  sensitivity_results <- cbind(n, m, J, baseline_constant, tau_w[j], tau_b[j], tau_b[j]/tau_w[j], design_power_t, score_power_predict_A, score_power_predict_tang)
}

colnames(sensitivity_results) <- c("n","m","J","baseline_constant","tau_w", "tau_b", "tau_ratio", "T", "S&M Score", "Tang Score")

## Step 4: plot results (Appendix F - Figure F.1) ##
sensitivity_colors <- c("#ffffe5", "#f7fcb9", "#d9f0a3","#addd8e","#78c679",
                        "#41ab5d","#238443","#006837","#004529", "#001e12")

sensitivity_results %>%
  as.data.frame() %>%
  pivot_longer(c("T", "S&M Score", "Tang Score"),
               names_to="paradigm", values_to = "power_predict") %>%
  mutate(paradigm=factor(paradigm, levels=c("T", "S&M Score", "Tang Score")),
         power_discrete = factor(case_when(power_predict <= 0.7 ~ "(0,0.7]",
                                           power_predict > 0.7 & power_predict <=0.75 ~ "(0.7,0.75]",
                                           power_predict > 0.75 & power_predict <=0.8 ~ "(0.75,0.8]",
                                           power_predict > 0.8 & power_predict <=0.85 ~ "(0.8,0.85]",
                                           power_predict > 0.85 & power_predict <=0.9 ~ "(0.85,0.9]",
                                           power_predict > 0.9 & power_predict <=0.92 ~ "(0.9,0.92]",
                                           power_predict > 0.92 & power_predict <=0.94 ~ "(0.92,0.94]",
                                           power_predict > 0.94 & power_predict <=0.96 ~ "(0.94,0.96]",
                                           power_predict > 0.96 & power_predict <=0.98 ~ "(0.96,0.98]",
                                           power_predict > 0.98 & power_predict <=1 ~ "(0.98,1]"))) %>%
  ggplot(aes(tau_ratio, tau_w, fill=power_discrete))+
  geom_tile(width=0.05)+
  facet_grid(vars(paradigm))+
  labs(x=TeX("$\\tau_b/\\tau_w$"),
       y=TeX("$\\tau_w$"),
       fill="Predicted Power")+
  scale_fill_manual(values=sensitivity_colors)
theme_minimal()+
  theme(strip.background =element_rect(fill="lightgray"))
  
#### CATHTAG SENSITIVITY ANALYSIS: n=24, beta=0.4, decreasing hazard (Appendix F - Figure F.3) #### 
## Step 1: set design parameters ##
m <- 35 # cluster size
n <- 24 # number of clusters
J <- 6 # number of time periods 
tau_w <- as.numeric(tau_scenarios[,1]) # within-period Kendall's tau
tau_b <- as.numeric(tai_scenarios[,2]) # between-period Kendall's tau

# set null and alternative effect sizes #
beta0 <- 0
betaA <- 0.4

## Step 2: specify censoring rates and baseline hazard function ##
# set maximum censoring time #
Cp <- 1

# calculate base baseline hazard that results in above administrative censoring rate #
lambdaDET <- function(Cp, pa){
  lambda0 <- (1/Cp)*(-log(pa))
  return(lambda0)
}
lambda0_base <- lambdaDET(Cp, pa)

# specify how does the baseline hazard before over time? #
baseline_constant <- 0 # specify whether/how baseline hazard changes with time (decreasing by 5% each time period here)
lambda0 <- lambda0_base+baseline_constant*seq(0,(J-1)) # final baseline hazard

# calculate probability of a cluster being on treatment #
pi_b <- c(0, rep((1/(J-1)), (J-1)))

## Step 3: loop through Kendall's tau scenarios ##
sensitivity_results <- matrix(NA,ncol=6, nrow=length(tau_w))

for(j in seq(length(tau_w))){
  
  ## Step 3.1: power for Wald t-test ##
  # estimate wald variance under the alternative #
  design_varA <- sandwich_var(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta=betaA)
  
  # calculate predicted power based on your sandwich variance and t test #
  design_power_t <- pt((abs(betaA)/(design_varA$var_beta_sqrt/sqrt(n)) - qt(0.975, n-1)), n-1)
  
  ## Step 3.2: power for score tests ##
  # estimate score variances under null and alternative hypotheses #
  design_var0_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta0=beta0, betaA=beta0)
  design_varA_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta0=beta0,betaA=betaA)
  
  # calculate predicted power based on Tang score test #
  score_power_predict_tang <-  pnorm(abs(design_varA_score$score)/(design_varA_score$B/sqrt(n)) - qnorm(0.975)*(design_var0_score$B/design_varA_score$B))
  # calculate predicted power based on S&M score test #
  score_power_predict_A <-  pnorm(abs(design_varA_score$score)/(design_varA_score$B/sqrt(n)) - qnorm(0.975))
  
  ## Step 3.3: store results ##
  sensitivity_results <- cbind(n, m, J, baseline_constant, tau_w[j], tau_b[j], tau_b[j]/tau_w[j], design_power_t, score_power_predict_A, score_power_predict_tang)
}

colnames(sensitivity_results) <- c("n","m","J","baseline_constant","tau_w", "tau_b", "tau_ratio", "T", "S&M Score", "Tang Score")

## Step 4: plot results (Appendix F - Figure F.3) ##
sensitivity_colors <- c("#ffffe5", "#f7fcb9", "#d9f0a3","#addd8e","#78c679",
                        "#41ab5d","#238443","#006837","#004529", "#001e12")

sensitivity_results %>%
  as.data.frame() %>%
  pivot_longer(c("T", "S&M Score", "Tang Score"),
               names_to="paradigm", values_to = "power_predict") %>%
  mutate(paradigm=factor(paradigm, levels=c("T", "S&M Score", "Tang Score")),
         power_discrete = factor(case_when(power_predict <= 0.7 ~ "(0,0.7]",
                                           power_predict > 0.7 & power_predict <=0.75 ~ "(0.7,0.75]",
                                           power_predict > 0.75 & power_predict <=0.8 ~ "(0.75,0.8]",
                                           power_predict > 0.8 & power_predict <=0.85 ~ "(0.8,0.85]",
                                           power_predict > 0.85 & power_predict <=0.9 ~ "(0.85,0.9]",
                                           power_predict > 0.9 & power_predict <=0.92 ~ "(0.9,0.92]",
                                           power_predict > 0.92 & power_predict <=0.94 ~ "(0.92,0.94]",
                                           power_predict > 0.94 & power_predict <=0.96 ~ "(0.94,0.96]",
                                           power_predict > 0.96 & power_predict <=0.98 ~ "(0.96,0.98]",
                                           power_predict > 0.98 & power_predict <=1 ~ "(0.98,1]"))) %>%
  ggplot(aes(tau_ratio, tau_w, fill=power_discrete))+
  geom_tile(width=0.05)+
  facet_grid(vars(paradigm))+
  labs(x=TeX("$\\tau_b/\\tau_w$"),
       y=TeX("$\\tau_w$"),
       fill="Predicted Power")+
  scale_fill_manual(values=sensitivity_colors)
  theme_minimal()+
  theme(strip.background =element_rect(fill="lightgray"))