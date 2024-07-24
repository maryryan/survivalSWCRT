#### LIBRARIES ####
library(tidyverse)
library(pracma)
library(survival)
library(stabledist)
library(latex2exp)
library(gee)
library(geesmv)

#### LOAD FUNCTIONS ####
source("functions-biometrics.R")

#### SIMULATIONS ####
set.seed(8888)

# set simulation scenarios #
simulation_scenarios <- list(betaA=c(rep(0.35,2),rep(0.4,4),rep(0.45,3),rep(0.5,3),0.55,0.6,rep(0.65,3),rep(0.7,3)),
                              n=c(rep(30,2),rep(21,2),rep(30,2),rep(21,2),30,rep(15,2),21,rep(15,2),rep(9,2),15,rep(9,3)),
                              m=c(40,50,40,50,15,25,25,40,15,40,50,15,25,15,40,50,15,25,40,50))
# no cluster size variation #
cv.c <- 0; cv.p <- 0

# number of study periods #
J <-4

# null hypothesis #
beta0 <- 0

# maximum censoring time #
Cp <- 1

# administrative censoring and net censoring rates #
pa <- 0.2; p0 <- 0.3

# how does the baseline hazard before over time? #
baseline_constant <- 0.2

# calculate base baseline hazard that results in above administrative censoring rate #
lambda0_base <- lambdaDET(Cp, pa)
lambda0 <- lambda0_base+baseline_constant*seq(0,(J-1))

# within- and between-period Kendall's tau #
tau_w <- 0.05; tau_b <- 0.01

# number of simulation iterations #
nrep <- 2000

# prepare results matrices #
sim_results_null <- matrix(NA, ncol=52, nrow=length(simulation_scenarios[[1]]))
sim_results_alt <- matrix(NA, ncol=54, nrow=length(simulation_scenarios[[1]]))

colnames(sim_results_null) <- c("n", "m", "periods", "tau_w", "tau_b","beta0", "emp coef",
                                "Emp. Naive SE", "ASE", "ESE",
                                "Design Naive SE", "Design Sandwich SE",
                                "BC MD SE", "BC KC SE", "BC FG SE",
                                "ASE Bias", "ESE Bias",
                                "BC0 Bias", "BC MB Bias", "BC KC Bias", "BC FG Bias",
                                "Test Stat Naive", "Test Stat Sandwich", 
                                "Test Stat MD", "Test Stat KC", "Test Stat FG",
                                "Z Type I Naive", "Z Type I Sandwich",
                                "Z Type I MD", "Z Type I KC", "Z Type I FG",
                                "T Type I Naive", "T Type I Sandwich",
                                "T Type I MD", "T Type I KC", "T Type I FG",
                                "Z type I GEE Indep", "Z type I GEE Exch",
                                "Z type I GEE Indep FG","Z type I GEE Indep KC","Z type I GEE Indep MD",
                                "Z type I GEE Exch FG","Z type I GEE Exch KC","Z type I GEE Exch MD",
                                "T type I GEE Indep", "T type I GEE Exch",
                                "T type I GEE Indep FG","T type I GEE Indep KC","T type I GEE Indep MD",
                                "T type I GEE Exch FG","T type I GEE Exch KC","T type I GEE Exch MD")

colnames(sim_results_alt) <- c("n", "m", "periods", "tau_w", "tau_b","betaA","emp coef",
                               "Emp. Naive SE", "ASE", "ESE",
                               "Design Naive SE", "Design Sandwich SE",
                               "BC MD SE", "BC KC SE", "BC FG SE",
                               "ASE Bias", "ESE Bias",
                               "BC0 Bias", "BC MB Bias", "BC KC Bias", "BC FG Bias",
                               "Test Stat Naive", "Test Stat Sandwich", 
                               "Test Stat MD", "Test Stat KC", "Test Stat FG",
                               "Z Power Design", "T Power Design",
                               "Z Power Naive", "Z Power Sandwich",
                               "Z Power MD", "Z Power KC", "Z Power FG",
                               "T Power Naive", "T Power Sandwich",
                               "T Power MD", "T Power KC", "T Power FG",
                               "Z Power GEE Indep", "Z Power GEE Exch",
                               "Z power GEE Indep FG","Z power GEE Indep KC","Z power GEE Indep MD",
                               "Z power GEE Exch FG","Z power GEE Exch KC","Z power GEE Exch MD",
                               "T power GEE Indep", "T power GEE Exch",
                               "T power GEE Indep FG","T power GEE Indep KC","T power GEE Indep MD",
                               "T power GEE Exch FG","T power GEE Exch KC","T power GEE Exch MD")

# begin simulations #
for(s in seq(length(simulation_scenarios[[1]]))){
  
  # set timing of intervention crossover for each cluster #
  k <- rep( seq(2,J), (n[s]/(J-1)) )
  
  # calculate probability of a cluster being on treatment #
  pi_b <- c(0, rep((1/(J-1)), (J-1)))
  
  # pull simulation parameters #
  n <- simulation_scenarios[["n"]][s]
  m <- simulation_scenarios[["m"]][s]
  betaA <- simulation_scenarios[["betaA"]][s]
  
  # simulate and analyze survival data under null and alternative distributions #
  empirical_var0 <- surSIMULATESW_with_gee(n=n, m=m, cv.c=cv.c, cv.p=cv.p, 
                                           lambda0_c=baseline_constant, 
                                           beta=beta0, pa=pa, p0=p0,
                                           tau_b=tau_b, tau_w=tau_w, Cp=Cp, nrep=nrep, k=k, p.max=J)
  empirical_varA <-surSIMULATESW_with_gee(n=n, m=m, cv.c=cv.c, cv.p=cv.p, 
                                          lambda0_c=baseline_constant, 
                                          beta=betaA, pa=pa, p0=p0,
                                          tau_b=tau_b, tau_w=tau_w, Cp=Cp, nrep=nrep, k=k, p.max=J)
  
  # predict variance under null and alternative distributions#
  design_var0 <- sandwich_var(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta=beta0)
  design_varA <- sandwich_var(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta=betaA)
  
  # calculate predicted power based on your sandwich variance and z test #
  design_power_z <- pnorm(abs(betaA)/design_varA$se_beta - qnorm(0.975))
  
  # calculate predicted power based on your sandwich variance and t test #
  design_power_t <- pt((abs(betaA)/design_varA$se_beta - qt(0.975, n-2)), n-2) 
  
  ## Null results ##
  ASE0 <- mean(empirical_var0[,4])
  ESE0 <- sd(empirical_var0[,1])
  naive0 <- mean(empirical_var0[,3])
  BC0_avg0 <- mean(empirical_var0[,21])
  BC_MD_avg0 <- mean(empirical_var0[,22])
  BC_KC_avg0 <- mean(empirical_var0[,23])
  BC_FG_avg0 <- mean(empirical_var0[,24])
  
  # calculate bias #
  ASE_bias0 <- design_var0$se_beta/ASE0
  ESE_bias0 <- design_var0$se_beta/ESE0
  
  BC0_bias0 <- design_var0$se_beta/BC0_avg0
  BC_MD_bias0 <- design_var0$se_beta/BC_MD_avg0
  BC_KC_bias0 <- design_var0$se_beta/BC_KC_avg0
  BC_FG_bias0 <- design_var0$se_beta/BC_FG_avg0
  
  # calculate test statistics #
  testStat_naive0 <- (empirical_var0[,1]-beta0)/empirical_var0[,3]
  testStat_robust0 <- (empirical_var0[,1]-beta0)/empirical_var0[,4]
  testStat_MD0 <- (empirical_var0[,1]-beta0)/empirical_var0[,22]
  testStat_KC0 <- (empirical_var0[,1]-beta0)/empirical_var0[,23]
  testStat_FG0 <- (empirical_var0[,1]-beta0)/empirical_var0[,24]
  
  testStat_gee_indep <- empirical_var0[,33]/empirical_var0[,36]
  testStat_gee_indep_fg <- empirical_var0[,33]/empirical_var0[,38]
  testStat_gee_indep_kc <- empirical_var0[,33]/empirical_var0[,39]
  testStat_gee_indep_md <- empirical_var0[,33]/empirical_var0[,40]
  testStat_gee_exch <- empirical_var0[,41]/empirical_var0[,44]
  testStat_gee_exch_fg <- empirical_var0[,41]/empirical_var0[,46]
  testStat_gee_exch_kc <- empirical_var0[,41]/empirical_var0[,47]
  testStat_gee_exch_md <- empirical_var0[,41]/empirical_var0[,48]
  
  testStat_results0 <- cbind(testStat_naive0, testStat_robust0,
                             testStat_MD0, testStat_KC0, testStat_FG0)
  
  testStat_gee_results0 <- cbind(testStat_gee_indep, testStat_gee_indep_fg, testStat_gee_indep_kc, testStat_gee_indep_md,
                                 testStat_gee_exch, testStat_gee_exch_fg, testStat_gee_exch_kc, testStat_gee_exch_md)
  
  
  # calculate type I error based on z-test #
  z_test_typeI <- apply(testStat_results0,2, function(x){
    
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qnorm(0.975), 1, 0 )
    }
    return(mean(reject))
  })
  
  z_test_gee_typeI <- apply(testStat_gee_results0,2, function(x){
    
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qnorm(0.975), 1, 0 )
    }
    return(mean(reject))
  })
  
  
  # calculate type I error based on t-test with df=n - 1#
  t_test_typeI <- apply(testStat_results0,2, function(x){
    
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-2), 1, 0 )
    }
    return(mean(reject))
  })
  
  t_test_gee_typeI <- apply(testStat_gee_results0,2, function(x){
    
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-2), 1, 0 )
    }
    return(mean(reject))
  })
  
  # aggregate null simulation results #
  sim_results_null[s,] <- cbind(n,m,J,tau_w,tau_b,beta0,
                                      mean(empirical_var0[,1]),
                                      naive0, ASE0, ESE0,
                                      design_var0$se_naive, design_var0$se_beta,
                                      BC_MD_avg0, BC_KC_avg0, BC_FG_avg0,
                                      ASE_bias0, ESE_bias0,
                                      BC0_bias0, BC_MD_bias0, BC_KC_bias0, BC_FG_bias0,
                                      matrix(apply(testStat_results0,2,mean),nrow=1),
                                      matrix(z_test_typeI, nrow=1), matrix(t_test_typeI, nrow=1),
                                      matrix(z_test_gee_typeI, nrow=1),matrix(t_test_gee_typeI, nrow=1))
  ## Alternative results ##
  ASEA <- mean(empirical_varA[,4])
  ESEA <- sd(empirical_varA[,1])
  naiveA <- mean(empirical_varA[,3])
  BC0_avgA <- mean(empirical_varA[,21])
  BC_MD_avgA <- mean(empirical_varA[,22])
  BC_KC_avgA <- mean(empirical_varA[,23])
  BC_FG_avgA <- mean(empirical_varA[,24])
  
  # calculate bias #
  ASE_biasA <- design_varA$se_beta/ASEA
  ESE_biasA <- design_varA$se_beta/ESEA
  
  BC0_biasA <- design_varA$se_beta/BC0_avgA
  BC_MD_biasA <- design_varA$se_beta/BC_MD_avgA
  BC_KC_biasA <- design_varA$se_beta/BC_KC_avgA
  BC_FG_biasA <- design_varA$se_beta/BC_FG_avgA
  
  # calculate test statistics #
  testStat_naiveA <- (empirical_varA[,1]-beta0)/empirical_varA[,3]
  testStat_robustA <- (empirical_varA[,1]-beta0)/empirical_varA[,4]
  testStat_MDA <- (empirical_varA[,1]-beta0)/empirical_varA[,22]
  testStat_KCA <- (empirical_varA[,1]-beta0)/empirical_varA[,23]
  testStat_FGA <- (empirical_varA[,1]-beta0)/empirical_varA[,24]
  
  testStat_gee_indepA <- empirical_varA[,33]/empirical_varA[,36]
  testStat_gee_indepA_fg <- empirical_varA[,33]/empirical_varA[,38]
  testStat_gee_indepA_kc <- empirical_varA[,33]/empirical_varA[,39]
  testStat_gee_indepA_md <- empirical_varA[,33]/empirical_varA[,40]
  testStat_gee_exchA <- empirical_varA[,41]/empirical_varA[,44]
  testStat_gee_exchA_fg <- empirical_varA[,41]/empirical_varA[,46]
  testStat_gee_exchA_kc <- empirical_varA[,41]/empirical_varA[,47]
  testStat_gee_exchA_md <- empirical_varA[,41]/empirical_varA[,48]
  
  testStat_resultsA <- cbind(testStat_naiveA, testStat_robustA,
                             testStat_MDA, testStat_KCA, testStat_FGA)
  
  testStat_gee_resultsA <- cbind(testStat_gee_indepA, testStat_gee_indepA_fg, testStat_gee_indepA_kc, testStat_gee_indepA_md,
                                 testStat_gee_exchA, testStat_gee_exchA_fg, testStat_gee_exchA_kc, testStat_gee_exchA_md)
  
  
  # calculate empirical power based on z-test #
  z_test_power <- apply(testStat_resultsA,2, function(x){
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qnorm(0.975), 1, 0 )
    }
    mean(reject)
    
  })
  
  z_test_gee_power <- apply(testStat_gee_resultsA,2, function(x){
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qnorm(0.975), 1, 0 )
    }
    mean(reject)
    
  })
  
  # calculate empirical power based on t-test with df=n - 1#
  t_test_power <- apply(testStat_resultsA,2, function(x){
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-2), 1, 0 )
    }
    mean(reject)
    
  })
  
  t_test_gee_power <- apply(testStat_gee_resultsA,2, function(x){
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-2), 1, 0 )
    }
    mean(reject)
    
  })
  
  ## aggregate alternative simulation results #
  sim_results_alt[s,] <- cbind(n,m,J,tau_w,tau_b, betaA,mean(empirical_varA[,1]),
                                     naiveA, ASEA, ESEA,
                                     design_varA$se_naive, design_varA$se_beta,
                                     BC_MD_avgA, BC_KC_avgA, BC_FG_avgA,
                                     ASE_biasA, ESE_biasA,
                                     BC0_biasA, BC_MD_biasA, BC_KC_biasA, BC_FG_biasA,
                                     matrix(apply(testStat_resultsA,2,mean),nrow=1),
                                     design_power_z, design_power_t,
                                     matrix(z_test_power,nrow=1), matrix(t_test_power,nrow=1),
                                     matrix(score_test_power,nrow=1),
                                     matrix(z_test_gee_power,nrow=1), matrix(t_test_gee_power,nrow=1))
  
  print(paste0(Sys.time(), ": Done with beta=", betaA,", J=",J,", n=", n, ", m=", m))
  
}

#### RESULTS ####
## type I error ##
sim_results_null %>% 
  as.data.frame() %>% 
  janitor::clean_names() %>% #glimpse()
  dplyr::select(-c("ase", "ese", "z_type_i_naive","z_type_i_sandwich",
                   "z_type_i_md","z_type_i_fg","z_type_i_kc",
                   "emp_coef","emp_naive_se", "bc_md_se", "bc_kc_se", "bc_fg_se",
                   "test_stat_naive","test_stat_sandwich","test_stat_md",
                   "test_stat_kc","test_stat_fg","beta0"
  )) %>% #glimpse()
  rename(t_type_i_robust = t_type_i_sandwich,
         gee_t_indep_type_i_robust = t_type_i_gee_indep,
         gee_t_indep_type_i_md = t_type_i_gee_indep_md,
         gee_t_indep_type_i_kc = t_type_i_gee_indep_kc,
         gee_t_indep_type_i_fg = t_type_i_gee_indep_fg,
         gee_t_exch_type_i_robust = t_type_i_gee_exch,
         gee_t_exch_type_i_md = t_type_i_gee_exch_md,
         gee_t_exch_type_i_kc = t_type_i_gee_exch_kc,
         gee_t_exch_type_i_fg = t_type_i_gee_exch_fg) %>% 
  pivot_longer(-c("n","m","periods", "tau_w", "tau_b"),
               names_to="decision_method", values_to = "type_i_error") %>%
  mutate(decision_distribution=factor(case_when(
    grepl("^t_", decision_method) ~ "Survival T",
    grepl("^gee_t_indep", decision_method) ~ "GEE T Indep.",
    grepl("^gee_t_exch", decision_method) ~ "GEE T Exch."),
    levels=c("Survival T", "GEE T Indep.", "GEE T Exch.")),
    decision_method = str_extract(decision_method, "(?<=_)[^_]+(?=$)"),#get the word after the last underscore
    decision_method = as.factor(decision_method)) %>% 
  group_by(decision_distribution, decision_method, n,m) %>% 
  summarize(mean(type_i_error)) %>%
  arrange(desc(n),m) %>% 
  print(n=300)

## power ##
sim_results_alt %>% 
  as.data.frame() %>% 
  janitor::clean_names() %>%
  dplyr::select(c("n","m","periods", "tau_w", "tau_b","beta_a",
                  "t_power_sandwich", "t_power_md", "t_power_kc", "t_power_fg",
                  "t_power_gee_indep",
                  "t_power_gee_indep_md", "t_power_gee_indep_kc", "t_power_gee_indep_fg",
                  "t_power_gee_exch",
                  "t_power_gee_exch_md", "t_power_gee_exch_kc", "t_power_gee_exch_fg")) %>% 
  rename(t_power_robust = t_power_sandwich,
    gee_t_indep_power_robust = t_power_gee_indep,
    gee_t_indep_power_md = t_power_gee_indep_md,
    gee_t_indep_power_kc = t_power_gee_indep_kc,
    gee_t_indep_power_fg = t_power_gee_indep_fg,
    gee_t_exch_power_robust = t_power_gee_exch,
    gee_t_exch_power_md = t_power_gee_exch_md,
    gee_t_exch_power_kc = t_power_gee_exch_kc,
    gee_t_exch_power_fg = t_power_gee_exch_fg) %>% 
  pivot_longer(!(c("n","m","periods", "tau_w", "tau_b","beta_a")), names_to="decision_method", values_to="power_error") %>% 
  mutate(decision_distribution=factor(case_when(
    grepl("^t_", decision_method) ~ "Survival T",
    grepl("^gee_t_indep", decision_method) ~ "GEE T Indep.",
    grepl("^gee_t_exch", decision_method) ~ "GEE T Exch."),
    levels=c("Survival T", "GEE T Indep.", "GEE T Exch.")),
    decision_method = str_extract(decision_method, "(?<=_)[^_]+(?=$)"),#get the word after the last underscore
    decision_method = as.factor(decision_method)) %>% 
  na.omit() %>% 
  dplyr::filter(m%in% c(15,25,50)) %>% 
  arrange(desc(n), desc(m)) %>% 
  print(n=150)
