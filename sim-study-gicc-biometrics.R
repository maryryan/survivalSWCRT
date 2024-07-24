#### LIBRARIES ####
library(tidyverse)
library(pracma)
library(survival)
library(stabledist)
library(latex2exp)
library(patchwork)

#### LOAD SIMULATION SCENARIOS ####
sim_scenarios_fixed_tauw <- read.table("gicc-tau-scenarios-fixed-tauw.txt",sep=",",header = T)
sim_scenarios_fixed_taub <- read.table("gicc-tau-scenarios-fixed-taub.txt",sep=",",header = T)

#### LOAD FUNCTIONS ####
source("functions-biometrics.R")

#### Fixed tau_w simulation scenarios ####
# set design parameters #
m <- as.numeric(sim_scenarios_fixed_tauw[,1])
n <- as.numeric(sim_scenarios_fixed_tauw[,2])
J <- as.numeric(sim_scenarios_fixed_tauw[,3])
tau_w <- as.numeric(sim_scenarios_fixed_tauw[,4])
tau_b <- as.numeric(sim_scenarios_fixed_tauw[,5])

# set null and alternative effect sizes #
beta0 <- 0
betaA <- 0.4

# maximum censoring time #
Cp <- 1

# administrative censoring and net censoring rates #
pa <- 0.2; p0 <- 0.3

# calculate base baseline hazard that results in above administrative censoring rate #
lambda0_base <- lambdaDET(Cp, pa)

# how does the baseline hazard before over time? #
baseline_constant <- as.numeric(sim_scenarios_fixed_tauw[,6])

# prepare results matrix #
gICC_results_fixed_tauw <- matrix(NA, nrow=length(m), ncol=8)
colnames(gICC_results) <- c("n","m","J","baseline_constant","tau_w", "tau_b", "gICC_w", "gICC_b")

# begin simulations #
for(j in seq(length(tau_w))){
  
  # how does the baseline hazard before over time? #
  lambda0 <- lambda0_base + baseline_constant[j]*seq(0,(J[j]-1))
  
  # set timing of intervention crossover for each cluster #
  k <- rep( seq(2,J[j]), (n[j]/(J[j]-1)) )
  
  # calculate probability of a cluster being on treatment #
  pi_b <- c(0, rep((1/(J[j]-1)), (J[j]-1)))

  # estimate wald variance under the alternative #
  design_varA <- sandwich_var(m=m[j], J=J[j], lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta=betaA)
  
  # collect results #
  gICC_results_fixed_tauw[j,] <- cbind(n[j], m[j], J[j], baseline_constant[j], tau_w[j], tau_b[j], design_varA$icc_w, design_varA$icc_b)

}

#### Fixed tau_b simulation scenarios ####
# set design parameters #
m <- as.numeric(sim_scenarios_fixed_taub[,1])
n <- as.numeric(sim_scenarios_fixed_taub[,2])
J <- as.numeric(sim_scenarios_fixed_taub[,3])
tau_w <- as.numeric(sim_scenarios_fixed_taub[,4])
tau_b <- as.numeric(sim_scenarios_fixed_taub[,5])

# set null and alternative effect sizes #
beta0 <- 0
betaA <- 0.4

# maximum censoring time #
Cp <- 1

# administrative censoring and net censoring rates #
pa <- 0.2; p0 <- 0.3

# calculate base baseline hazard that results in above administrative censoring rate #
lambda0_base <- lambdaDET(Cp, pa)

# how does the baseline hazard before over time? #
baseline_constant <- as.numeric(sim_scenarios_fixed_taub[,6])

# prepare results matrix #
gICC_results_fixed_taub <- matrix(NA, nrow=length(m), ncol=8)
colnames(gICC_results) <- c("n","m","J","baseline_constant","tau_w", "tau_b", "gICC_w", "gICC_b")

# begin simulations #
for(j in seq(length(tau_w))){
  
  # how does the baseline hazard before over time? #
  lambda0 <- lambda0_base+baseline_constant[j]*seq(0,(J[j]-1))
  
  # set timing of intervention crossover for each cluster #
  k <- rep( seq(2,J[j]), (n[j]/(J[j]-1)) )
  
  # calculate probability of a cluster being on treatment #
  pi_b <- c(0, rep((1/(J[j]-1)), (J[j]-1)))
  
  # estimate wald variance under the alternative #
  design_varA <- sandwich_var(m=m[j], J=J[j], lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta=betaA)
  
  # collect results #
  gICC_results_fixed_taub[j,] <- cbind(n[j], m[j], J[j], baseline_constant[j], tau_w[j], tau_b[j], design_varA$icc_w, design_varA$icc_b)
  
}


#### Result plots ####
line_colors <- c("#7bccc4", "#43a2ca", "#0868ac")

fixed_tauw_plot <- ICC_results_fixed_tauw %>% 
  as.data.frame() %>% 
  mutate(J=factor(J, levels=c(3,6,11)),
         tau_w=factor(tau_w, levels=c(0.05, 0.1, 0.2),
                      labels=c('tau[w]*"=0.05"',
                               'tau[w]*"=0.1"',
                               'tau[w]*"=0.2"'))
  ) %>%
  ggplot(aes(tau_b, gICC_b, color=J))+
  geom_line(linewidth=1.05)+
  facet_wrap(vars(tau_w),
             labeller=label_parsed,
             scales="free")+
  labs(x=TeX("$\\tau_b$"),
       y="Between-period g-ICC",
       color="Time Periods")+ 
  ylim(0,0.1)+
  scale_color_manual(values=line_colors)+
  theme_minimal()+
  theme(strip.background =element_rect(fill="lightgray"),
        strip.text=element_text(size=20),
        text=element_text(size=16),
        axis.text=element_text(size=16),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))


fixed_taub_plot <- ICC_results_fixed_taub %>% 
  as.data.frame() %>% 
  mutate(J=factor(J, levels=c(3,6,11)),
         tau_b=factor(tau_b, levels=c(0, 0.04, 0.09),
                      labels=c('tau[b]*"=0"',
                               'tau[b]*"=0.04"',
                               'tau[b]*"=0.09"'))
  ) %>%
  ggplot(aes(tau_w, gICC_w, color=J))+
  geom_line(linewidth=1.05)+
  facet_wrap(vars(tau_b),
             labeller=label_parsed,
             scales="free")+
  labs(x=TeX("$\\tau_w$"),
       y="Within-period g-ICC",
       color="Time Periods")+
  ylim(0,0.25)+
  scale_color_manual(values=line_colors)+
  theme_minimal()+
  theme(strip.background =element_rect(fill="lightgray"),
        strip.text=element_text(size=20),
        text=element_text(size=16),
        axis.text=element_text(size=16),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position="none")

fixed_taub_plot/fixed_tauw_plot