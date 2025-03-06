#### LIBRARIES ####
library(tidyverse)
library(pracma)
library(survival)
library(stabledist)
library(latex2exp)
library(patchwork)

#### LOAD SIMULATION SCENARIOS ####
simulation_scenarios <- read.table("main_sim_alt_scenarios.txt",sep=",",header = T)

#### LOAD FUNCTIONS ####
source("functions-biometrics.R")

## function to determine base baseline hazard ##
lambdaDET <- function(Cp, pa){#kappa, pa){
  lambda0 <- (1/Cp)*(-log(pa))#^(1/kappa)
  return(lambda0)
}

#### SIMULATIONS ####

# no cluster size variation #
cv.c <- 0; cv.p <- 0

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

# number of simulation iterations #
nrep <- 2000

# begin simulations #
for(s in seq(nrow(simulation_scenarios))){
  
  # pull simulation parameters #
  J <- simulation_scenarios[s,"periods"]
  n <- simulation_scenarios[s,"n"]
  m <- simulation_scenarios[s,"m"]
  betaA <- simulation_scenarios[s,"betaA"]
  tau_w <- simulation_scenarios[s,"tau_w"]
  tau_b <- simulation_scenarios[s,"tau_b"]
  
  # calculate period-specific baseline hazards #
  lambda0 <- lambda0_base+baseline_constant*seq(0,(J-1))
  
  # set timing of intervention crossover for each cluster #
  k <- rep( seq(2,J), (simulation_scenarios[[1]][s]/(J-1)) )
  
  # calculate probability of a cluster being on treatment #
  pi_b <- c(0, rep((1/(J-1)), (J-1)))
  
  # simulate data under null and alternative hypotheses #
  empirical_var0 <- surSIMULATESW(n=n, J=m, cv.c=cv.c, cv.p=cv.p, 
                                  lambda0_c=baseline_constant, 
                                  beta=beta0, pa=pa, p0=p0,
                                  tau_b=tau_b, tau_w=tau_w, Cp=Cp, nrep=nrep, k=k, p.max=J)
  empirical_varA <- surSIMULATESW(n=n, J=m, cv.c=cv.c, cv.p=cv.p, 
                                  lambda0_c=baseline_constant, 
                                  beta=betaA, pa=pa, p0=p0,
                                  tau_b=tau_b, tau_w=tau_w, Cp=Cp, nrep=nrep, k=k, p.max=J)
  
  # estimate wald variance under the null and alternative #
  design_var0 <- sandwich_var(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta=beta0)
  design_varA <- sandwich_var(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta=betaA)
  
  # calculate predicted power based on your sandwich variance and t test #
  design_power_t <- pt((abs(betaA)/(design_varA$var_beta_sqrt/sqrt(n)) - qt(0.975, n-2)), n-2)
  
  # estimate score variances under null and alternative hypotheses #
  design_var0_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta0=beta0, betaA=beta0)
  design_varA_score <- sandwich_var_score(m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta0=beta0,betaA=betaA)
  
  # calculate predicted power based on score tests #
  design_power_score_tang_B <- pnorm(abs(design_varA_score$score)/(design_varA_score$B/sqrt(n)) - qnorm(0.975)*(design_var0_score$B/design_varA_score$B))
  design_power_score_SM_B <- pnorm(abs(design_varA_score$score)/(design_varA_score$B/sqrt(n)) - qnorm(0.975))
  designA_A_mod <- (design_varA_score$A/sqrt(n))*sqrt((n-1)/n)
  design0_A_mod <- (design_var0_score$A/sqrt(n))*sqrt((n-1)/n)
  design_power_score_mod <- pnorm(abs(design_varA_score$score)/designA_A_mod- qnorm(0.975)*(design0_A_mod/designA_A_mod))
  
  ## Null results ##
  ASE0 <- mean(empirical_var0[,4])
  ESE0 <- sd(empirical_var0[,1])
  naive0 <- mean(empirical_var0[,3])
  BC0_avg0 <- mean(empirical_var0[,21])
  BC_MD_avg0 <- mean(empirical_var0[,22])
  BC_KC_avg0 <- mean(empirical_var0[,23])
  BC_FG_avg0 <- mean(empirical_var0[,24])
  
  # calculate bias #
  ASE_bias0 <- design_var0$var_beta_sqrt/ASE0
  ESE_bias0 <- design_var0$var_beta_sqrt/ESE0
  
  BC0_bias0 <- design_var0$var_beta_sqrt/BC0_avg0
  BC_MD_bias0 <- design_var0$var_beta_sqrt/BC_MD_avg0
  BC_KC_bias0 <- design_var0$var_beta_sqrt/BC_KC_avg0
  BC_FG_bias0 <- design_var0$var_beta_sqrt/BC_FG_avg0
  
  # calculate test statistics #
  testStat_naive0 <- (empirical_var0[,1]-beta0)/empirical_var0[,3]
  testStat_robust0 <- (empirical_var0[,1]-beta0)/empirical_var0[,4]
  testStat_MD0 <- (empirical_var0[,1]-beta0)/empirical_var0[,22]
  testStat_KC0 <- (empirical_var0[,1]-beta0)/empirical_var0[,23]
  testStat_FG0 <- (empirical_var0[,1]-beta0)/empirical_var0[,24]
  
  testStat_score_robust0 <-empirical_var0[,25]/(empirical_var0[,26])#/sqrt(n))
  testStat_score_robust0_a <- empirical_var0[,27]/(empirical_var0[,26])#/sqrt(n))#uses nom as score and BC-B0 as variance
  testStat_score_robust_mod0 <-empirical_var0[,25]/(empirical_var0[,26]*sqrt((n-1)/n))
  testStat_score_robust_mod0_a <- empirical_var0[,27]/(empirical_var0[,26]*sqrt((n-1)/n))#uses nom as score
  
  testStat_results0 <- cbind(testStat_naive0, testStat_robust0,
                             testStat_MD0, testStat_KC0, testStat_FG0)
  
  testStat_score_results0 <- cbind(
    testStat_score_robust0_a, 
    testStat_score_robust_mod0_a
  ) 
  
  # calculate type I error based on z-test #
  z_test_typeI <- apply(testStat_results0,2, function(x){
    
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qnorm(0.975), 1, 0 )
    }
    return(mean(reject))
  })
  
  # calculate type I error based on t-test with df=n - 2#
  t_test_typeI <- apply(testStat_results0,2, function(x){
    
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-2), 1, 0 )
    }
    return(mean(reject))
  })
  
  # calculate type I error for score tests #
  score_test_typeI <- apply(testStat_score_results0,2, function(x){
    
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qnorm(0.975), 1, 0 )
    }
    return(mean(reject))
  })
  
  sim_results_null <- cbind(n,m,J,tau_w,tau_b,beta0,
                            mean(empirical_var0[,1]),
                            naive0, ASE0, ESE0,
                            design_var0$var_naive, design_var0$var_beta_sqrt,
                            BC_MD_avg0, BC_KC_avg0, BC_FG_avg0,
                            ASE_bias0, ESE_bias0,
                            BC0_bias0, BC_MD_bias0, BC_KC_bias0, BC_FG_bias0,
                            matrix(apply(testStat_results0,2,mean),nrow=1),
                            matrix(z_test_typeI, nrow=1), matrix(t_test_typeI, nrow=1),
                            matrix(apply(testStat_score_results0,2,mean),nrow=1),
                            matrix(score_test_typeI, nrow=1))
  
  ## Alternative results ##
  ASEA <- mean(empirical_varA[,4])
  ESEA <- sd(empirical_varA[,1])
  naiveA <- mean(empirical_varA[,3])
  BC0_avgA <- mean(empirical_varA[,21])
  BC_MD_avgA <- mean(empirical_varA[,22])
  BC_KC_avgA <- mean(empirical_varA[,23])
  BC_FG_avgA <- mean(empirical_varA[,24])
  
  # calculate bias #
  ASE_biasA <- design_varA$var_beta_sqrt/ASEA
  ESE_biasA <- design_varA$var_beta_sqrt/ESEA
  
  BC0_biasA <- design_varA$var_beta_sqrt/BC0_avgA
  BC_MD_biasA <- design_varA$var_beta_sqrt/BC_MD_avgA
  BC_KC_biasA <- design_varA$var_beta_sqrt/BC_KC_avgA
  BC_FG_biasA <- design_varA$var_beta_sqrt/BC_FG_avgA
  
  # calculate test statistics #
  testStat_naiveA <- (empirical_varA[,1]-beta0)/empirical_varA[,3]
  testStat_robustA <- (empirical_varA[,1]-beta0)/empirical_varA[,4]
  testStat_MDA <- (empirical_varA[,1]-beta0)/empirical_varA[,22]
  testStat_KCA <- (empirical_varA[,1]-beta0)/empirical_varA[,23]
  testStat_FGA <- (empirical_varA[,1]-beta0)/empirical_varA[,24]
  
  testStat_score_robustA <-empirical_varA[,25]/(empirical_varA[,26])#/sqrt(n))
  testStat_score_robustA_a <- empirical_varA[,27]/(empirical_varA[,26])#/sqrt(n))#uses nom as score
  testStat_score_robust_modA <-empirical_varA[,25]/(empirical_varA[,26]*sqrt((n-1)/n))
  testStat_score_robust_modA_a <- empirical_varA[,27]/(empirical_varA[,26]*sqrt((n-1)/n))#uses nom as score
  
  testStat_resultsA <- cbind(testStat_naiveA, testStat_robustA,
                             testStat_MDA, testStat_KCA, testStat_FGA)
  
  testStat_score_resultsA <- cbind(
    testStat_score_robustA_a,
    testStat_score_robust_modA_a
  ) 
  
  # calculate empirical power based on z-test #
  z_test_power <- apply(testStat_resultsA,2, function(x){
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qnorm(0.975), 1, 0 )
    }
    mean(reject)
    
  })
  
  # calculate empirical power based on t-test with df=n - 2#
  t_test_power <- apply(testStat_resultsA,2, function(x){
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-2), 1, 0 )
    }
    mean(reject)
    
  })
  
  # calculate empirical power based on score test #
  score_test_power <- apply(testStat_score_resultsA,2, function(x){
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qnorm(0.975), 1, 0 )
    }
    mean(reject)
    
  })
  
  sim_results_alt <- cbind(n,m,J,tau_w,tau_b, betaA,
                           mean(empirical_varA[,1]),
                           naiveA, ASEA, ESEA,
                           design_varA$var_naive, design_varA$var_beta_sqrt,
                           BC_MD_avgA, BC_KC_avgA, BC_FG_avgA,
                           ASE_biasA, ESE_biasA,
                           BC0_biasA, BC_MD_biasA, BC_KC_biasA, BC_FG_biasA,
                           matrix(apply(testStat_resultsA,2,mean),nrow=1),
                           design_power_t,#design_power_t_n2,#design_power_score,
                           matrix(z_test_power,nrow=1), matrix(t_test_power,nrow=1),
                           design_power_score_tang_B,design_power_score_SM_B,
                           design_power_score_mod,
                           matrix(apply(testStat_score_resultsA,2,mean),nrow=1),
                           matrix(score_test_power,nrow=1))
  
  
  print(paste0(Sys.time(), ": Done with beta=", betaA,", J=",J,", n=", n, ", m=", m))
  
  
}

colnames(sim_results_null) <- c("n", "m", "periods", "tau_w", "tau_b","beta0", "emp coef",
                                "Emp. Naive SE", "ASE", "ESE",
                                "Design Naive SE", "Design Sandwich SE",
                                "BC MD SE", "BC KC SE", "BC FG SE",
                                "ASE Bias", "ESE Bias",
                                "BC0 Bias", "BC MD Bias", "BC KC Bias", "BC FG Bias",
                                "Test Stat Naive", "Test Stat Sandwich", 
                                "Test Stat MD", "Test Stat KC", "Test Stat FG",
                                "Z Type I Naive", "Z Type I Sandwich",
                                "Z Type I MD", "Z Type I KC", "Z Type I FG",
                                "T Type I Naive", "T Type I Sandwich",
                                "T Type I MD", "T Type I KC", "T Type I FG",
                                "Score Test Stat Sandwich a",
                                "Score Test Stat Sandwich Mod a",
                                "Score Type I Sandwich a",
                                "Score Type I Sandwich Mod a")

colnames(sim_results_alt) <- c("n", "m", "periods", "tau_w", "tau_b","betaA","emp coef",
                               "Emp. Naive SE", "ASE", "ESE",
                               "Design Naive SE", "Design Sandwich SE",
                               "BC MD SE", "BC KC SE", "BC FG SE",
                               "ASE Bias", "ESE Bias",
                               "BC0 Bias", "BC MB Bias", "BC KC Bias", "BC FG Bias",
                               "Test Stat Naive", "Test Stat Sandwich", 
                               "Test Stat MD", "Test Stat KC", "Test Stat FG",
                               "T Power Design", 
                               "Z Power Naive", "Z Power Sandwich",
                               "Z Power MD", "Z Power KC", "Z Power FG",
                               "T Power Naive", "T Power Sandwich",
                               "T Power MD", "T Power KC", "T Power FG",
                               "Score Power Design Tang B", "Score Power Design SM B",
                               "Score Mod Power Design",
                               "Score Test Stat Sandwich a",
                               "Score Test Stat Sandwich Mod a",
                               "Score Power Sandwich a",
                               "Score Power Sandwich Mod a")

#### RESULTS ####
# type I error plot #
sim_null_results %>% 
  as.data.frame() %>% 
  janitor::clean_names() %>%
  dplyr::select(c("n","m","periods", "tau_w", "tau_b",
                  "t_type_i_sandwich", "t_type_i_md", "t_type_i_kc", "t_type_i_fg",
                  "score_type_i_sandwich_a",
                  "score_type_i_sandwich_mod_a")) %>% 
  rename(score_type_i_sandwich_a_mod = score_type_i_sandwich_mod_a) %>% #,
  pivot_longer(!(c("n","m","periods", "tau_w", "tau_b")), names_to="decision_method", values_to="type_i_error") %>% 
  mutate(decision_distribution=factor(case_when(grepl("^t", decision_method) ~ "T",
                                           grepl("^s", decision_method) ~ "Robust Score"),
                                 levels=c("T", "Robust Score", "Z")),
    decision_method = str_extract(decision_method, "(?<=_)[^_]+(?=$)"),#get the word after the last underscore
    decision_method = as.factor(decision_method)) %>%
  # dplyr::filter(tau_w == tau_scenarios[1,1], tau_b==tau_scenarios[1,2]
  # ) %>% 
  group_by(periods, m, n, decision_distribution, decision_method) %>% 
  summarize(type_i_error=mean(type_i_error, na.rm=T)) %>%
  as.data.frame() %>% 
  ungroup() %>% 
  mutate(decision_method = dplyr::recode(decision_method,
                                         sandwich="Robust SE",
                                         fg="FG Bias Correction",
                                         kc="KC Bias Correction",
                                         md="MD Bias Correction",
                                         a="Non-Modified Score",
                                         mod="(n-1)/n Modified Score"),
         decision_method = factor(decision_method,
                                  levels=c("Robust SE",
                                           "FG Bias Correction",
                                           "KC Bias Correction",
                                           "MD Bias Correction",
                                           "Non-Modified Score",
                                           "(n-1)/n Modified Score")),
         valid_test = case_when(type_i_error < 0.06 ~ 1,
                                type_i_error >= 0.06 ~ 0))%>%
  ggplot(aes(x=n,y=type_i_error))+
  geom_point(aes(shape=decision_method, color=as.factor(m)))+
  ylim(0,0.1)+
  ylab("Empirical Type I Error")+
  xlab("Number of Clusters (n)")+
  geom_hline(yintercept = 0.06, color="dark gray", linetype="dashed", size=0.7)+
  geom_hline(yintercept = 0.05, color="red", linetype="dashed", size=0.7)+
  geom_hline(yintercept = 0.04, color="dark gray", linetype="dashed", size=0.7)+
  theme_bw()+
  guides(color=guide_legend(title="Cluster-Period Size (m)", order=1),
         shape=guide_legend(title="Test Statistic", order=0))+
  facet_grid(vars(decision_distribution),
             vars(periods), labeller=as_labeller(
               c("3"="J=3", "4"="J=4", "5"="J=5", "6"="J=6",
                 "T"="T", "Robust Score"="Robust Score")
             ))+
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c"))

# Differential Power  plot #
m_unique <- sort(unique(score_alt_results$m))
score_wald_alt_plot_diff <- vector(mode="list", length=length(m_unique))

for(i in seq(length(m_unique))){
  score_wald_alt_plot_diff[[i]] <- sim_alt_results %>% 
    as.data.frame() %>% 
    janitor::clean_names() %>% 
    dplyr::select(c("n","m","periods", "tau_w", "tau_b","beta_a",
                    "t_power_design",
                    "score_power_design_new_var_a",
                    "score_power_design_new_tang",
                    "t_power_sandwich", "t_power_md", "t_power_kc", "t_power_fg",
                    "score_power_sandwich_a","score_power_sandwich_mod_a"
    )) %>%
    mutate(Tang_power_sandwich_a = score_power_sandwich_a,
           Tang_power_sandwich_a_mod = score_power_sandwich_mod_a) %>%
    rename(SM_power_sandwich_a = score_power_sandwich_a,
           SM_power_sandwich_a_mod = score_power_sandwich_mod_a,
           SM_score_power_design = score_power_design_new_var_a,
           tang_score_power_design = score_power_design_new_tang
    ) %>%
    pivot_longer(!(c("n","m","periods", "tau_w", "tau_b","beta_a",
                     "t_power_design",
                     "SM_score_power_design",
                     "tang_score_power_design"
    )),
    names_to="decision_method", values_to="power") %>% 
    mutate(decision_distribution=factor(case_when(
      grepl("^t_", decision_method) ~ "T",
      grepl("^SM", decision_method) ~ "S&M Score",
      grepl("^Tang", decision_method) ~ "Tang Score"),
      levels=c("T", "S&M Score", "Tang Score")),
      decision_method = str_extract(decision_method, "(?<=_)[^_]+(?=$)"),#get the word after the last underscore
      decision_method = as.factor(decision_method),
      power_difference = case_when(
        decision_distribution == "T" ~ power - t_power_design,
        decision_distribution == "S&M Score" ~ power - SM_score_power_design,
        decision_distribution == "Tang Score" ~ power -tang_score_power_design
      )) %>% 
    dplyr::filter(tau_w == tau_scenarios[1,1], tau_b==tau_scenarios[1,2],
                  m==m_unique[i]
    ) %>% 
    mutate(decision_method = dplyr::recode(decision_method,
                                           sandwich="Robust SE",
                                           fg="FG Bias Correction",
                                           kc="KC Bias Correction",
                                           md="MD Bias Correction",
                                           a="Non-Modified Score",
                                           mod="(n-1)/n Modified Score"),
           decision_method = factor(decision_method,
                                    levels=c("Robust SE",
                                             "FG Bias Correction",
                                             "KC Bias Correction",
                                             "MD Bias Correction",
                                             "Non-Modified Score",
                                             "(n-1)/n Modified Score")))%>%
    ggplot(aes(x=n,y=power_difference))+
    geom_point(aes(shape=decision_method, color=beta_a))+
    geom_hline(yintercept = 0, color="red", linetype="dashed", linewidth=0.7)+
    geom_hline(yintercept = 0.017, color="dark gray", linetype="dashed", linewidth=0.7)+
    geom_hline(yintercept = -0.017, color="dark gray", linetype="dashed", linewidth=0.7)+
    ylim(-0.075,0.12)+
    theme_bw()+
    facet_grid(vars(decision_distribution),
               vars(periods), labeller=as_labeller(
                 c("3"="J=3", "4"="J=4", "5"="J=5", "6"="J=6",
                   "T"="T", "S&M Score"="S&M Score", "Tang Score"="Tang Score")
               ))+
    ggtitle(paste0("m = ", m_unique[i]))+
    labs(y="Empirical Power - Design Power",
         x="Number of Clusters (n)",
         color="Treatment Effect",
         shape="Test Statistic")+
    scale_color_continuous(limits=c(0.25,0.7))+
    scale_x_continuous(breaks=c(10,20,30))+
    guides(color=guide_colorbar(order=1),
           shape=guide_legend(order=0))+
    theme(plot.title=element_text(hjust=0.5))
  if(!(i %in% c(2,4))){
    score_wald_alt_plot_diff[[i]] <- score_wald_alt_plot_diff[[i]] + theme(legend.position = "none",strip.background.y = element_blank(), strip.text.y = element_blank())
  }
  if(!(i %in% c(1,3))){
    score_wald_alt_plot_diff[[i]] <- score_wald_alt_plot_diff[[i]] + theme(axis.title.y = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
  }
  if(i != 2 ){
    score_wald_alt_plot_diff[[i]] <- score_wald_alt_plot_diff[[i]] + theme(legend.position = "none")
  }
}

(score_wald_alt_plot_diff[[1]] | score_wald_alt_plot_diff[[2]]) / (score_wald_alt_plot_diff[[3]]|score_wald_alt_plot_diff[[4]]) #+ plot_annotation(title = "Cluster-Period Size (m)",theme = theme(plot.title = element_text(hjust = 0.4)))
