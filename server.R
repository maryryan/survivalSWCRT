#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
#library(plotly)
library(DT)
library(patchwork)
library(latex2exp)
#library(RColorBrewer)
library(shinybrowser)

library(pracma)
library(survival)
library(stabledist)

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  options(shiny.sanitize.errors = FALSE)
  # set maximum between-period kendall's tau to be smaller than within-period tau #
  observe(updateNumericInput(session, "tau_b", max = tau_w() - .Machine$double.eps,))
  
  #### WALD POWER FUNCTIONS ####
  q_0 <- function(s, beta, z_ij, j,lambda0, tau,desmat_prob){
    # censoring distribution
    G_s <- 1-s/tau
    
    # W_js (ratio of s_1j/s_0j)
    pi_b_j <- sum(desmat_prob[,j])
    s_1j <- pi_b_j * exp(beta - exp(beta)*lambda0*s)
    s_0j <- pi_b_j * exp(beta - exp(beta)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
    W_js <- s_1j/s_0j
    
    # Density functions (negative derivative of the survival function)
    # density on treatment arm # 
    f_ij_s <- lambda0*exp(z_ij*beta)*exp(-lambda0*s*exp(z_ij*beta))
    
    
    # final output
    q0 <- G_s*(z_ij-W_js)^2*f_ij_s
    
    return(q0)
  }
  
  q_1 <- function(s, t, beta, z_ij, z_il, j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    
    # transform kendall's tau into correlation parameters
    theta_b <- 1-tau_kb
    theta_w <- 1-tau_kw
    
    # censoring distribution - independent
    G_st <- (1-s/tau)*(1-t/tau)
    
    # Ws
    pi_b_j <- sum(desmat_prob[,j])
    s_1js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s)
    s_0js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
    W_js <- s_1js/s_0js
    
    pi_b_l <- sum(desmat_prob[,l])
    s_1lt <- pi_b_l * exp(beta - exp(beta)*lambda0_l*t)
    s_0lt <- pi_b_l * exp(beta - exp(beta)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
    W_lt <- s_1lt/s_0lt
    
    # Density functions
    # these will be the same if j == l 
    lambda_ij <- lambda0_j*exp(z_ij*beta)
    lambda_il <- lambda0_l*exp(z_il*beta)
    
    # different periods #
    if(j != l){
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^theta_b )
      
      f_ijl_st <- lambda_ij*lambda_il*survivor_jl*(lambda_ij*s)^((1/theta_b)-1)*(lambda_il*t)^((1/theta_b)-1)*( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^(2*theta_b-2)*( 1+((1/theta_b)-1)*((lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b))^(-theta_b) )
      
      
    }else{
      # same period #
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^theta_w )
      
      f_ijl_st <- lambda_ij*lambda_il*survivor_jl*(lambda_ij*s)^((1/theta_w)-1)*(lambda_il*t)^((1/theta_w)-1)*( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^(2*theta_w-2)*( 1+((1/theta_w)-1)*((lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w))^(-theta_w) )
      
    }
    
    q1 <- G_st*(z_ij - W_js)*(z_il - W_lt)*f_ijl_st
    
    
    return(q1)
    
  }
  
  q_2 <- function(s, t, beta, z_ij, z_il, j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    
    # transform kendall's tau into correlation parameters
    theta_b <- 1-tau_kb
    theta_w <- 1-tau_kw
    
    # censoring distribution
    G_st <- (1-s/tau)*(1-t/tau)
    
    # Ws
    pi_b_j <- sum(desmat_prob[,j])
    s_1js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s)
    s_0js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
    W_js <- s_1js/s_0js
    
    pi_b_l <- sum(desmat_prob[,l])#sum(pi_b[1:l])
    s_1lt <- pi_b_l * exp(beta - exp(beta)*lambda0_l*t)
    s_0lt <- pi_b_l * exp(beta - exp(beta)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
    W_lt <- s_1lt/s_0lt
    
    # Density functions
    lambda_ij <- lambda0_j*exp(z_ij*beta)
    lambda_il <- lambda0_l*exp(z_il*beta)
    
    
    # different periods #
    if(j != l){
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^theta_b )
      
      survivor_div <- -(1/s)*(lambda_ij*s)^(1/theta_b)*survivor_jl*((lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b))^(theta_b - 1)
      
    }else{
      # same period #
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^theta_w )
      
      survivor_div <- -(1/s)*(lambda_ij*s)^(1/theta_w)*survivor_jl*((lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w))^(theta_w - 1)
      
    }
    
    # Final output
    q2 <-  G_st*(z_ij-W_js)*(z_il-W_lt)*(-survivor_div)*lambda_il
    
    return(q2)
    
  }
  
  q_3 <- function(s, t, beta, z_ij, z_il, j ,l, lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    
    # transform kendall's tau into correlation parameters
    theta_b <- 1-tau_kb
    theta_w <- 1-tau_kw
    
    # censoring distribution
    G_st <- (1-s/tau)*(1-t/tau)
    
    # Ws
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s)
    s_0js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
    W_js <- s_1js/s_0js
    
    pi_b_l <- sum(desmat_prob[,l])#sum(pi_b[1:l])
    s_1lt <- pi_b_l * exp(beta - exp(beta)*lambda0_l*t)
    s_0lt <- pi_b_l * exp(beta - exp(beta)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
    W_lt <- s_1lt/s_0lt
    
    # Density functions
    lambda_ij <- lambda0_j*exp(z_ij*beta)
    lambda_il <- lambda0_l*exp(z_il*beta)
    
    
    # different periods #
    if(j != l){
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^theta_b )
      
      survivor_div <- -(1/t)*(lambda_il*t)^(1/theta_b)*survivor_jl*((lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b))^(theta_b - 1)
      
    }
    else{
      # same period #
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^theta_w )
      
      survivor_div <- -(1/t)*(lambda_il*t)^(1/theta_w)*survivor_jl*((lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w))^(theta_w - 1)
      
    }
    
    # Final output
    q3 <-  G_st*(z_ij-W_js)*(z_il-W_lt)*(-survivor_div)*lambda_ij
    
    return(q3)
    
  } 
  
  q_4 <- function(s,t, beta, z_ij, z_il, j, l, lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    # transform kendall's tau into correlation parameters
    theta_b <- 1-tau_kb
    theta_w <- 1-tau_kw
    
    # censoring distribution
    G_st <- (1-s/tau)*(1-t/tau)
    
    # Ws
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s)
    s_0js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
    W_js <- s_1js/s_0js
    
    pi_b_l <- sum(desmat_prob[,l])#sum(pi_b[1:l])
    s_1lt <- pi_b_l * exp(beta - exp(beta)*lambda0_l*t)
    s_0lt <- pi_b_l * exp(beta - exp(beta)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
    W_lt <- s_1lt/s_0lt
    
    # Density functions
    lambda_ij <- lambda0_j*exp(z_ij*beta)
    lambda_il <- lambda0_l*exp(z_il*beta)
    
    
    # different periods #
    if(j != l){
      
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^theta_b )
      
    }else{
      # same periods #
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^theta_w )
      
    }
    
    q4 <-  G_st*(z_ij-W_js)*(z_il-W_lt)*survivor_jl*lambda_ij*lambda_il
    
    return(q4)
    
  }
  
  # for the bread A #
  V_0 <- function(s, beta, z_ij, j, lambda0, tau, desmat_prob){
    
    # censoring distribution
    G_s <- 1-s/tau
    
    # W_js (ratio of s_1j/s_0j)
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1j <- s_2j <- pi_b_j * exp(beta - exp(beta)*lambda0*s)
    s_0j <- pi_b_j * exp(beta - exp(beta)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
    
    # Density functions (negative derivative of the survival function)
    f_ij_s <- lambda0*exp(z_ij*beta)*exp(-lambda0*s*exp(z_ij*beta))
    
    v <- G_s*((s_2j/s_0j) - (s_1j^2)/(s_0j^2))*f_ij_s
    
    return(v)
    
  }
  
  
  Q <- function(beta, z_ij, z_il, j, l, lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    # function that integrates and adds together the 4 q components for the non-same subject portion of B
    Q <- dblquad(q_1,0, tau,0,tau,dim=2, beta=beta, z_ij=z_ij, z_il=z_il, j=j, l=l, lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb) -
      dblquad(q_2,0,tau,0,tau, dim=2, beta=beta, z_ij=z_ij, z_il=z_il, j=j, l=l,lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb) -
      dblquad(q_3,0,tau,0,tau,dim=2, beta=beta, z_ij=z_ij, z_il=z_il, j=j, l=l, lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb) +
      dblquad(q_4,0,tau,0,tau, dim=2, beta=beta, z_ij=z_ij, z_il=z_il, j=j, l=l,lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb)
    
    return(Q)
  }
  
  H0 <- function(j, beta, lambda0, tau,desmat_prob){
    # function that integrates and takes the expectation for the same subject, cluster, period portion of B
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    
    H0 <- pi_b_j*pracma::integral(fun=q_0,xmin=0,xmax=tau,beta=beta, z_ij=1,j=j,lambda0=lambda0, tau=tau,desmat_prob=desmat_prob) +
      (1-pi_b_j)*pracma::integral(fun=q_0,xmin=0,xmax=tau, beta=beta, z_ij=0,j=j,lambda0=lambda0, tau=tau,desmat_prob=desmat_prob)
    
    return(H0)
  }
  
  H1 <- function(j, l, beta, lambda0_j, lambda0_l, tau,desmat_prob,tau_kw, tau_kb){
    # function that takes the expectation and adds together the integrated q components for the non-same subject portion of B
    pi_b_min <- ifelse( j < l, sum(desmat_prob[,j]), sum(desmat_prob[,l]) )#ifelse( j < l, sum(pi_b[1:j]), sum(pi_b[1:l]) )
    pi_b_max <- ifelse( j > l, sum(desmat_prob[,j]), sum(desmat_prob[,l]) )#ifelse( j > l, sum(pi_b[1:j]), sum(pi_b[1:l]) )
    if( j > l){
      pi_b_mid <- sum(desmat_prob[,j]) - sum(desmat_prob[,l])#sum(pi_b[(l+1):j])
    }else if(j < l){ 
      pi_b_mid <- sum(desmat_prob[,l]) - sum(desmat_prob[,j])#sum(pi_b[(j+1):l]) 
    }else if(j == l){
      pi_b_mid <- 0
    }
    
    H1 <-  (1-pi_b_max)*Q(beta, 0,0,j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb)
    
    if(j > 1 | l > 1){
      
      if(j > l){
        H1 <- H1 + pi_b_mid*Q(beta, 1,0,j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb)
      }else if(j < l){
        H1 <- H1 + pi_b_mid*Q(beta, 0,1,j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb)
      }
      
      if(j > 1 & l > 1){
        H1 <- H1 + pi_b_min*Q(beta, 1,1,j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb)
      }
    }
    
    return(H1)
  }
  
  H2 <- function(j, beta, lambda0, tau, desmat_prob){
    # function that integrates and takes the expectation of the A components
    
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    
    H2 <- pi_b_j*pracma::integral(fun=V_0, xmin=0, xmax=tau, beta=beta, z_ij=1, j=j, lambda0=lambda0, tau=tau, desmat_prob=desmat_prob) + 
      (1-pi_b_j)*pracma::integral(fun=V_0, xmin=0, xmax=tau, beta=beta, z_ij=0, j=j, lambda0=lambda0, tau=tau, desmat_prob=desmat_prob)
    
    return(H2)
  }
  
  sandwich_var <- function(m,J, lambda0,tau,desmat_prob, tau_kw, tau_kb, beta=beta){
    H0_sum <- H2_sum <- rep(0, J)
    #H_score_sum <- rep(NA,J)
    H1_sum_eq <- H1_sum_uneq <- H0_sum_iccb <- matrix(NA, nrow=J, ncol=J)
    for(j in seq(J)){
      
      H0_sum[j] <- H0(j=j, beta=beta, lambda0=lambda0[j], tau=tau, desmat_prob=desmat_prob)
      H2_sum[j] <- H2(j=j, beta=beta, lambda0=lambda0[j], tau=tau, desmat_prob=desmat_prob)
      
      for(l in seq(J)){
        if(j==l){
          H1_sum_eq[j,l] <- H1(j=j,l=l, beta=beta, lambda0_j=lambda0[j], lambda0_l=lambda0[l], tau=tau, desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb)
        }else{
          H1_sum_uneq[j,l] <- H1(j=j,l=l, beta=beta, lambda0_j=lambda0[j], lambda0_l=lambda0[l], tau=tau, desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb)
          H0_sum_iccb[j,l] <-  H0_sum[j]*H0_sum[l]
        }
        
      }
    }
    
    # meat of sandwich 
    B <- m*sum(H0_sum) + m*(m-1)*sum(H1_sum_eq, na.rm=T) + (m^2)*sum(H1_sum_uneq, na.rm=T)
    B_iccs <-  sum(H0_sum)*( m + m*(m-1)*tau_kw + m^2*(J-1)*tau_kb )
    
    # sandwich bread 
    A <- m*sum(H2_sum)
    
    # MAKE THAT SANDWICH #
    var_beta <- solve(A) %*% B %*% solve(A)
    var_beta_iccs <- solve(A) %*% B_iccs %*% solve(A)
    
    # ICCs #
    icc_w <- sum(H1_sum_eq, na.rm=T)/sum(H0_sum)
    icc_b <- sum(H1_sum_uneq, na.rm=T)/((J-1)*sum(H0_sum))
    
    return(list(var_beta_sqrt=sqrt(var_beta),
                var_naive=sqrt(solve(A)),
                B=sqrt(B),
                A=sqrt(A),
                icc_w=icc_w, icc_b=icc_b))
    
  }
  
  #### SCORE POWER FUNCTIONS ####
  q_0_score <- function(s, beta0, betaA, z_ij, j,lambda0, tau,desmat_prob){
    # censoring distribution
    G_s <- 1-s/tau
    
    # W_js (ratio of s_1j/s_0j)
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1j <- pi_b_j * exp(beta0 - exp(beta0)*lambda0*s)
    s_0j <- pi_b_j * exp(beta0 - exp(beta0)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
    W_js <- s_1j/s_0j
    
    # Density functions (negative derivative of the survival function)
    # density on treatment arm # 
    f_ij_s <- lambda0*exp(z_ij*beta0)*exp(-lambda0*s*exp(z_ij*betaA))
    
    # final output
    q0 <- G_s*(z_ij-W_js)^2*f_ij_s
    
    return(q0)
  }
  
  q_1_score <- function(s, t, beta0, betaA, z_ij, z_il, j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    
    # transform kendall's tau into correlation parameters
    theta_b <- 1-tau_kb
    theta_w <- 1-tau_kw
    
    # censoring distribution - independent
    G_st <- (1-s/tau)*(1-t/tau)
    
    # Ws
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1js <- pi_b_j * exp(beta0 - exp(beta0)*lambda0_j*s)
    s_0js <- pi_b_j * exp(beta0 - exp(beta0)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
    W_js <- s_1js/s_0js
    
    pi_b_l <- sum(desmat_prob[,l])#sum(pi_b[1:l])
    s_1lt <- pi_b_l * exp(beta0 - exp(beta0)*lambda0_l*t)
    s_0lt <- pi_b_l * exp(beta0 - exp(beta0)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
    W_lt <- s_1lt/s_0lt
    
    # Density functions
    # these will be the same if j == l 
    lambda_ij_0 <- lambda0_j*exp(z_ij*beta0)
    lambda_il_0 <- lambda0_l*exp(z_il*beta0)
    lambda_ij <- lambda0_j*exp(z_ij*betaA)
    lambda_il <- lambda0_l*exp(z_il*betaA)
    
    # different periods #
    if(j != l){
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^theta_b )
      
      f_ijl_st <- lambda_ij*lambda_il*survivor_jl*(lambda_ij*s)^((1/theta_b)-1)*(lambda_il*t)^((1/theta_b)-1)*( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^(2*theta_b-2)*( 1+((1/theta_b)-1)*((lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b))^(-theta_b) )
      
    }else{
      # same period #
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^theta_w )
      
      f_ijl_st <- lambda_ij*lambda_il*survivor_jl*(lambda_ij*s)^((1/theta_w)-1)*(lambda_il*t)^((1/theta_w)-1)*( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^(2*theta_w-2)*( 1+((1/theta_w)-1)*((lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w))^(-theta_w) )
    }
    
    q1 <- G_st*(z_ij - W_js)*(z_il - W_lt)*f_ijl_st
    
    return(q1)
    
  }
  
  q_2_score <- function(s, t, beta0, betaA, z_ij, z_il, j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    
    # transform kendall's tau into correlation parameters
    theta_b <- 1-tau_kb
    theta_w <- 1-tau_kw
    
    # censoring distribution
    G_st <- (1-s/tau)*(1-t/tau)
    
    # Ws
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1js <- pi_b_j * exp(beta0 - exp(beta0)*lambda0_j*s)
    s_0js <- pi_b_j * exp(beta0 - exp(beta0)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
    W_js <- s_1js/s_0js
    
    pi_b_l <- sum(desmat_prob[,l])#sum(pi_b[1:l])
    s_1lt <- pi_b_l * exp(beta0 - exp(beta0)*lambda0_l*t)
    s_0lt <- pi_b_l * exp(beta0 - exp(beta0)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
    W_lt <- s_1lt/s_0lt
    
    # Density functions
    lambda_ij_0 <- lambda0_j*exp(z_ij*beta0)
    lambda_il_0 <- lambda0_l*exp(z_il*beta0)
    lambda_ij <- lambda0_j*exp(z_ij*betaA)
    lambda_il <- lambda0_l*exp(z_il*betaA)
    
    # different periods #
    if(j != l){
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^theta_b )
      
      survivor_div <- -(1/s)*(lambda_ij*s)^(1/theta_b)*survivor_jl*((lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b))^(theta_b - 1)
      
    }else{
      # same period #
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^theta_w )
      
      survivor_div <- -(1/s)*(lambda_ij*s)^(1/theta_w)*survivor_jl*((lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w))^(theta_w - 1)
      
    }
    
    # Final output
    q2 <-  G_st*(z_ij-W_js)*(z_il-W_lt)*(-survivor_div)*lambda_il
    
    return(q2)
    
  }
  
  q_3_score <- function(s, t, beta0, betaA, z_ij, z_il, j ,l, lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    
    # transform kendall's tau into correlation parameters
    theta_b <- 1-tau_kb
    theta_w <- 1-tau_kw
    
    # censoring distribution
    G_st <- (1-s/tau)*(1-t/tau)
    
    # Ws
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1js <- pi_b_j * exp(beta0 - exp(beta0)*lambda0_j*s)
    s_0js <- pi_b_j * exp(beta0 - exp(beta0)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
    W_js <- s_1js/s_0js
    
    pi_b_l <- sum(desmat_prob[,l])#sum(pi_b[1:l])
    s_1lt <- pi_b_l * exp(beta0 - exp(beta0)*lambda0_l*t)
    s_0lt <- pi_b_l * exp(beta0 - exp(beta0)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
    W_lt <- s_1lt/s_0lt
    
    # Density functions
    lambda_ij_0 <- lambda0_j*exp(z_ij*beta0)
    lambda_il_0 <- lambda0_l*exp(z_il*beta0)
    lambda_ij <- lambda0_j*exp(z_ij*betaA)
    lambda_il <- lambda0_l*exp(z_il*betaA)
    
    
    # different periods #
    if(j != l){
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^theta_b )
      
      survivor_div <- -(1/t)*(lambda_il*t)^(1/theta_b)*survivor_jl*((lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b))^(theta_b - 1)
      
    }
    else{
      # same period #
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^theta_w )
      
      survivor_div <- -(1/t)*(lambda_il*t)^(1/theta_w)*survivor_jl*((lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w))^(theta_w - 1)
      
    }
    
    # Final output
    q3 <-  G_st*(z_ij-W_js)*(z_il-W_lt)*(-survivor_div)*lambda_ij
    
    return(q3)
    
  } 
  
  q_4_score <- function(s,t, beta0, betaA, z_ij, z_il, j, l, lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    # transform kendall's tau into correlation parameters
    theta_b <- 1-tau_kb
    theta_w <- 1-tau_kw
    
    # censoring distribution
    G_st <- (1-s/tau)*(1-t/tau)
    
    # Ws
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1js <- pi_b_j * exp(beta0 - exp(beta0)*lambda0_j*s)
    s_0js <- pi_b_j * exp(beta0 - exp(beta0)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
    W_js <- s_1js/s_0js
    
    pi_b_l <- sum(desmat_prob[,l])#sum(pi_b[1:l])
    s_1lt <- pi_b_l * exp(beta0 - exp(beta0)*lambda0_l*t)
    s_0lt <- pi_b_l * exp(beta0 - exp(beta0)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
    W_lt <- s_1lt/s_0lt
    
    # Density functions
    lambda_ij_0 <- lambda0_j*exp(z_ij*beta0)
    lambda_il_0 <- lambda0_l*exp(z_il*beta0)
    lambda_ij <- lambda0_j*exp(z_ij*betaA)
    lambda_il <- lambda0_l*exp(z_il*betaA)
    
    
    # different periods #
    if(j != l){
      
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_b) + (lambda_il*t)^(1/theta_b) )^theta_b )
      
    }else{
      # same periods #
      survivor_jl <- exp( -( (lambda_ij*s)^(1/theta_w) + (lambda_il*t)^(1/theta_w) )^theta_w )
      
    }
    
    # Final output
    q4 <-  G_st*(z_ij-W_js)*(z_il-W_lt)*survivor_jl*lambda_ij*lambda_il
    
    return(q4)
    
  }
  
  V_0_score <- function(s, beta0, betaA, z_ij, j, lambda0, tau, desmat_prob){
    
    # censoring distribution
    G_s <- 1-s/tau
    
    # W_js (ratio of s_1j/s_0j)
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1j <- s_2j <- pi_b_j * exp(beta0 - exp(beta0)*lambda0*s)
    s_0j <- pi_b_j * exp(beta0 - exp(beta0)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
    
    # Density functions (negative derivative of the survival function)
    f_ij_s <- lambda0*exp(z_ij*beta0)*exp(-lambda0*s*exp(z_ij*betaA))
    
    v <- G_s*((s_2j/s_0j) - (s_1j^2)/(s_0j^2))*f_ij_s#(W_js - W_js2)*f_ij_s
    
    return(v)
    
  }
  
  Q_score <- function(beta0, betaA, z_ij, z_il, j, l, lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb){
    # function that integrates and adds together the 4 q components for the non-same subject portion of B
    Q <- dblquad(q_1_score,0, tau,0,tau,dim=2, beta0=beta0, betaA=betaA, z_ij=z_ij, z_il=z_il, j=j, l=l, lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb) -
      dblquad(q_2_score,0,tau,0,tau, dim=2,beta0=beta0, betaA=betaA, z_ij=z_ij, z_il=z_il, j=j, l=l,lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb) -
      dblquad(q_3_score,0,tau,0,tau,dim=2, beta0=beta0, betaA=betaA, z_ij=z_ij, z_il=z_il, j=j, l=l, lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb) +
      dblquad(q_4_score,0,tau,0,tau, dim=2, beta0=beta0, betaA=betaA, z_ij=z_ij, z_il=z_il, j=j, l=l,lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb)
    
    return(Q)
  }
  
  H0_score <- function(j, beta0, betaA, lambda0, tau,desmat_prob){
    # function that integrates and takes the expectation for the same subject, cluster, period portion of B
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    
    H0 <- pi_b_j*pracma::integral(fun=q_0_score,xmin=0,xmax=tau,beta0=beta0, betaA=betaA, z_ij=1,j=j,lambda0=lambda0, tau=tau,desmat_prob=desmat_prob) +
      (1-pi_b_j)*pracma::integral(fun=q_0_score,xmin=0,xmax=tau, beta0=beta0, betaA=betaA, z_ij=0,j=j,lambda0=lambda0, tau=tau,desmat_prob=desmat_prob)
    
    return(H0)
  }
  
  H1_score <- function(j, l, beta0, betaA, lambda0_j, lambda0_l, tau,desmat_prob,tau_kw, tau_kb){
    # function that takes the expectation and adds together the integrated q components for the non-same subject portion of B
    pi_b_min <- ifelse( j < l, sum(desmat_prob[,j]), sum(desmat_prob[,l]) )#ifelse( j < l, sum(pi_b[1:j]), sum(pi_b[1:l]) )
    pi_b_max <- ifelse( j > l, sum(desmat_prob[,j]), sum(desmat_prob[,l]) )#ifelse( j > l, sum(pi_b[1:j]), sum(pi_b[1:l]) )
    if( j > l){
      pi_b_mid <- sum(desmat_prob[,j]) - sum(desmat_prob[,l])#sum(pi_b[(l+1):j])
    }else if(j < l){
      pi_b_mid <- sum(desmat_prob[,l]) - sum(desmat_prob[,j])#sum(pi_b[(j+1):l]) 
    }else if(j == l){
      pi_b_mid <- 0
    }
    
    H1 <-  (1-pi_b_max)*Q_score(beta0=beta0, betaA=betaA, 0,0,j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb)
    
    if(j > 1 | l > 1){
      
      if(j > l){
        H1 <- H1 + pi_b_mid*Q_score(beta0=beta0, betaA=betaA, 1,0,j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb)
      }else if(j < l){
        H1 <- H1 + pi_b_mid*Q_score(beta0=beta0, betaA=betaA, 0,1,j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb)
      }
      
      if(j > 1 & l > 1){
        H1 <- H1 + pi_b_min*Q_score(beta0=beta0, betaA=betaA, 1,1,j,l,lambda0_j, lambda0_l, tau,desmat_prob, tau_kw, tau_kb)
      }
    }
    
    return(H1)
  }
  H2_score <- function(j, beta0, betaA, lambda0, tau, desmat_prob){
    # function that integrates and takes the expectation of the A components
    
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    
    H2 <- pi_b_j*pracma::integral(fun=V_0_score, xmin=0, xmax=tau, beta0=beta0, betaA=betaA, z_ij=1, j=j, lambda0=lambda0, tau=tau, desmat_prob=desmat_prob) + 
      (1-pi_b_j)*pracma::integral(fun=V_0_score, xmin=0, xmax=tau, beta0=beta0, betaA=betaA, z_ij=0, j=j, lambda0=lambda0, tau=tau, desmat_prob=desmat_prob)
    
    return(H2)
  }
  
  score_exp <- function(s, beta0, betaA, z_ij, j,lambda0, tau,desmat_prob){
    # censoring distribution
    G_s <- 1-s/tau
    
    # W_js (ratio of s_1j/s_0j)
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    s_1j <- pi_b_j * exp(beta0 - exp(beta0)*lambda0*s)
    s_0j <- pi_b_j * exp(beta0 - exp(beta0)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
    W_js <- s_1j/s_0j
    
    # Density functions (negative derivative of the survival function)
    # density on treatment arm # 
    f_ij_s <- lambda0*exp(z_ij*beta0)*exp(-lambda0*s*exp(z_ij*betaA))
    
    
    # final output
    score_exp <- (G_s-f_ij_s)*(z_ij-W_js)
    return(score_exp)
  }
  
  H_score <- function(j, beta0, betaA, lambda0, tau,desmat_prob){
    # function that integrates and takes the expectation for the same subject, cluster, period portion of B
    pi_b_j <- sum(desmat_prob[,j])#sum(pi_b[1:j])
    
    H_score <- pi_b_j*pracma::integral(fun=score_exp,xmin=0,xmax=tau,beta0=beta0, betaA=betaA, z_ij=1,j=j,lambda0=lambda0, tau=tau,desmat_prob=desmat_prob) +
      (1-pi_b_j)*pracma::integral(fun=score_exp,xmin=0,xmax=tau, beta0=beta0, betaA=betaA, z_ij=0,j=j,lambda0=lambda0, tau=tau,desmat_prob=desmat_prob)
    
    return(H_score)
  }
  
  sandwich_var_score <- function(m,J, lambda0,tau,desmat_prob, tau_kw, tau_kb, beta0=beta0, betaA=betaA){
    H0_sum <- H2_sum <- rep(0, J)
    H_score_sum <- rep(NA,J)
    H1_sum_eq <- H1_sum_uneq <- matrix(NA, nrow=J, ncol=J)
    for(j in seq(J)){
      
      H0_sum[j] <- H0_score(j=j, beta0=beta0,betaA=betaA, lambda0=lambda0[j], tau=tau, desmat_prob=desmat_prob)
      H2_sum[j] <- H2_score(j=j, beta0=beta0,betaA=betaA, lambda0=lambda0[j], tau=tau, desmat_prob=desmat_prob)
      H_score_sum[j] <- H_score(j=j, beta0=beta0, betaA=betaA, lambda0=lambda0[j], tau=tau, desmat_prob=desmat_prob)
      for(l in seq(J)){
        if(j==l){
          H1_sum_eq[j,l] <- H1_score(j=j,l=l, beta0=beta0,betaA=betaA, lambda0_j=lambda0[j], lambda0_l=lambda0[l], tau=tau, desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb)
        }else{
          H1_sum_uneq[j,l] <- H1_score(j=j,l=l, beta0=beta0,betaA=betaA, lambda0_j=lambda0[j], lambda0_l=lambda0[l], tau=tau, desmat_prob=desmat_prob, tau_kw=tau_kw, tau_kb=tau_kb)
        }
        
      }
    }
    
    # meat of sandwich 
    B <- m*sum(H0_sum) + m*(m-1)*sum(H1_sum_eq, na.rm=T) + (m^2)*sum(H1_sum_uneq, na.rm=T)
    
    # sandwich bread 
    A <- m*sum(H2_sum)
    
    # MAKE THAT SANDWICH #
    var_beta <- solve(A) %*% B %*% solve(A)
    
    # also make the score #
    score <- m*sum(H_score_sum)
    score2 <- m*sum(tcrossprod(H_score_sum))
    
    score_var <- sqrt(A)
    
    return(list(var_beta_sqrt=sqrt(var_beta),
                var_naive=sqrt(solve(A)),
                B=sqrt(B), 
                score=score,
                score2=sqrt(score2), 
                A=score_var))
    
  }
  
  
  #### APP FUNCTIONS ####
  designMatrix <- function(design, periods=NULL, clusters=NULL, file=NULL){
    if(design == "balanced"){
      sequences <- periods-1
      
      if( clusters %% sequences != 0 ){
        stop("Must have equal number of clusters per sequence for design to be balanaced")
      }else{
        
        W0 <- matrix(0, ncol=periods, nrow=sequences)
        W0[upper.tri(W0)] <- 1
        W0[rep(1:nrow(W0), each = (clusters/sequences)),]
      }
      
      
    }else if(design == "unbalanced"){
      read.csv(file$datapath, header=FALSE)
    }
  }
  
  lambdaDET <- function(Cp, pa){#kappa, pa){
    if(pa == 0){
      pa <- 0+.Machine$double.eps
    }
    lambda0 <- (1/Cp)*(-log(pa))#^(1/kappa)
    return(lambda0)
  }
  
  #### OUTPUT ####
  ## upload own observes ##
  file1 <- reactive({input$file1})
  
  observeEvent(input$reset, {
    file1 <- NULL
    reset('file1')
  })
  

  ## reacts ##
  n_power <- eventReactive(input$submit, {
    input$n_power
    })
  
  balanced <- eventReactive(input$submit, {
    input$balanced
    })
  
  m <- eventReactive(input$submit, {
    input$m
    })
  
  n_use <- eventReactive(input$submit, {
    if( n_power() == "power" ){
      
      if(balanced() == "balanced"){
        
        input$n
        
      }else if(balanced() == "unbalanced"){
        
        if(is.null(file1())) stop("User needs to upload design matrix before the function can continue")
        desmat <- read.csv(file1()$datapath, header=FALSE)
        
        nrow(desmat)
        
      }
    }else{
      NULL
    }
    
  })
  
  J_use <- eventReactive(input$submit, {
    if( n_power() == "power"){
      if(balanced() == "balanced"){
        
        input$J
        
      }else if(balanced() == "unbalanced"){
        
        if(is.null(file1())) stop("User needs to upload design matrix before the function can continue")
        desmat <- read.csv(file1()$datapath, header=FALSE)
        
        ncol(desmat)
        
      } 
    }else{
      input$J
    }
    
  })
  
  desmat <- eventReactive(input$submit, {
    if( is.null(n_use()) == FALSE){
      if(balanced() == "balanced"){
        
        designMatrix(design = balanced(), periods = J_use(), clusters = n_use())
        
      }else if(balanced() == "unbalanced"){
        
        if(is.null(file1())) stop("User needs to upload design matrix before the function can continue")
        desmat <- read.csv(file1()$datapath, header=FALSE)
        desmat
        
      }
    }else{
      NULL
    }
  
  })
  
  desmat_prob <- reactive({
    if(n_power() == "n"){
      W0 <- matrix( 0, ncol=J_use(), nrow=(J_use()-1) )
      W0[upper.tri(W0)] <- 1
      W0*(1/(J_use()-1))
      #c(0, rep((1/(J_use()-1)), (J_use()-1)))
      
    }else{
      desmat()/n_use() 
    }
  })
  
  power <- eventReactive(input$submit, {
    input$power
    })
  
  tau_w <- reactive({#eventReactive(input$submit, {
    input$tau_w
    })
  tau_b <- reactive({#eventReactive(input$submit, {
    input$tau_b
    })
  
  #icc_w <- reactive({input$icc_w})
  #icc_b <- reactive({input$icc_b})
  
  # cv.c <- 0; cv.p <- 0
  
  betaA <- eventReactive(input$submit, {input$beta})
  pa <- eventReactive(input$submit,{input$admin_censor})
  # p0 <- reactive({input$other_censor})
  constant_baseline <- reactive({#eventReactive(input$submit, {
    input$constant_baseline
    })
  baseline_change <- reactive({#eventReactive(input$submit, {
    input$baseline_change
    })
  
  alpha <- eventReactive(input$submit, {input$sig})
  
  DoF <- eventReactive(input$submit, {
    if(input$DoF == "n1"){
      1
    }else if(input$DoF == "n2"){
      2
    }
  })
  
  ## reactive variance functions ##
  design_varA <- eventReactive(input$submit, {
    Cp <- 1
    beta0 <- 0
    lambda0_base <- lambdaDET(Cp, pa())
    
    if(constant_baseline() == "constant"){
      lambda0 <- lambda0_base+0*seq(0,(J_use()-1))
    }else{
      lambda0 <- lambda0_base+baseline_change()*seq(0,(J_use()-1))
    }
    
    
    #k <- rep( seq(2,J()), (n()/(J()-1)) )
    
    
    # wald stuff #
    sandwich_var(#n=n_use(),
                 m=m(),J=J_use(),
                 lambda0=lambda0,tau=Cp, desmat_prob=desmat_prob(),#pi_b=pi_b,
                 tau_kw=tau_w(), tau_kb=tau_b(), beta=betaA())
  })
  
  
  design_varA_score <- eventReactive(input$submit, {#reactive({
    Cp <- 1
    beta0 <- 0
    lambda0_base <- lambdaDET(Cp, pa())
    
    if(constant_baseline() == "constant"){
      lambda0 <- lambda0_base+0*seq(0,(J_use()-1))
    }else{
      lambda0 <- lambda0_base+baseline_change()*seq(0,(J_use()-1))
    }
    
    #k <- rep( seq(2,J()), (n()/(J()-1)) )
    #pi_b <- c(0, rep((1/(J()-1)), (J()-1)))
    
    sandwich_var_score(#n=n_use(),
                       m=m(),J=J_use(),lambda0=lambda0,tau=Cp, desmat_prob=desmat_prob(),#pi_b=pi_b,
                       tau_kw=tau_w(), tau_kb=tau_b(),
                       beta0=beta0,betaA=betaA())
  })
  
  design_var0_score <- eventReactive(input$submit, {#reactive({
    Cp <- 1
    beta0 <- 0
    lambda0_base <- lambdaDET(Cp, pa())
    
    if(constant_baseline() == "constant"){
      lambda0 <- lambda0_base+0*seq(0,(J_use()-1))
    }else{
      lambda0 <- lambda0_base+baseline_change()*seq(0,(J_use()-1))
    }
    
    #k <- rep( seq(2,J()), (n()/(J()-1)) )
    #pi_b <- c(0, rep((1/(J()-1)), (J()-1)))
    
    # score stuff #
    sandwich_var_score(#n=n_use(),
                       m=m(),J=J_use(),lambda0=lambda0,tau=Cp, desmat_prob=desmat_prob(),#pi_b=pi_b,
                       tau_kw=tau_w(), tau_kb=tau_b(),
                       beta0=beta0, betaA=beta0)
  })
  
  ## Render text first ##
  output$text_context <- renderText({
    if(input$submit == 0){
      paste0("Input your design parameters on the left and press the 'Update View' button to calculate power/number of clusters.")
    }else{
    
      if(n_power()=="power"){
        
          paste0("The predicted power of a ", balanced(), " SW-CRT with J=", J_use(), " periods, n = ", n_use(), " clusters, and m = ", m(), " participants per cluster-period is:") 
      
        }else{
          paste0("For a SW-CRT to obtain at least ", power()*100,"% power with J=", J_use(), " periods and m = ", m(), " participants per cluster-period, the study would need:")
        }
    }
    
  })
  output$text_t <- renderText({
    
    if(n_power()=="power"){
      
      # calculate predicted power based on your sandwich variance and t test #
      
      design_power_t <- pt((abs(betaA())/(design_varA()$var_beta_sqrt/sqrt(n_use())) - qt((1-alpha()/2), n_use()-DoF())), n_use()-DoF())
      
      # output results #
      paste0(round(design_power_t*100, 2), "% under the Wald t-testing paradigm,")
      
    }else{
      
      # calculate number of clusters under given power #
      z_typeI <- qnorm(1-alpha()/2)
      z_power <- qnorm(power())
      #t_typeI <- qt((1-alpha()/2), n_use()-1)
      #t_power <- qt(power(), n_use()-1)
      
      n_result <- ( design_varA()$var_beta_sqrt * (z_typeI + z_power) )^2/betaA()^2
      
      # output results #
      paste0("n=", ceiling(n_result), " clusters under the Wald z-testing paradigm,")
      
      
    }
    
    
  })
  
  output$text_smScore <- renderText({
    
    if(n_power()=="power"){
      
    # calculate predicted power based on score test #
    score_power_predict_A <-  pnorm(abs(design_varA_score()$score)/(design_varA_score()$B/sqrt(n_use())) - qnorm(1-alpha()/2))
    
    # output results #
    paste0(round(score_power_predict_A*100,2), "% under the Self and Mauritsen robust score testing paradigm, and")
    
    }else{
      
      # calculate number of clusters under given power #
      z_typeI <- qnorm(1-alpha()/2)
      z_power <- qnorm(power())
      
      # no need to do anything with n_use here bc score and B are output without that division
      n_result <- ( design_varA_score()$B * (z_typeI + z_power) )^2/abs(design_varA_score()$score)^2#/sqrt(n_use()) * (z_typeI + z_power)^2 )/abs(design_varA_score()$score)^2
      
      # output results #
      paste0("n=", ceiling(n_result), " clusters under the Self and Mauritsen robust score testing paradigm, and")
      
      
    }
    
  })
  
  output$text_tangScore <- renderText({
    
    if(n_power() == "power"){
      
      # calculate predicted power based on score test #
      score_power_predict_tang <-  pnorm(abs(design_varA_score()$score)/(design_varA_score()$B/sqrt(n_use())) - qnorm(1-alpha()/2)*(design_var0_score()$B/design_varA_score()$B))
      
      # output results #
      paste0(round(score_power_predict_tang*100,2), "% under the Tang robust score testing paradigm.")
      
    }else{
      
      # calculate number of clusters under given power #
      z_typeI <- qnorm(1-alpha()/2)
      z_power <- qnorm(power())
      
      n_result <- ( design_varA_score()$B * (z_power + z_typeI*(design_var0_score()$B/design_varA_score()$B)) )^2/abs(design_varA_score()$score)^2#/sqrt(n_use()) * ((z_typeI + z_power)*(1/design_var0_score()$B/design_varA_score()$B))^2 )/abs(design_varA_score()$score)^2
      
      # output results #
      paste0("n=", ceiling(n_result), " clusters under the Tang robust score testing paradigm.")
      
    }
    
    
  })
  
  output$text_ICCs <- renderText({
  
    # output results #
    paste0("The within-period generalized ICC is estimated to be ", round(design_varA()$icc_w,2), " and the between-period generalized ICC is estimated to be ", round(design_varA()$icc_b, 2),".")
    
  })
  
  output$design_matrix <- renderTable({
    if( is.null(n_use()) == TRUE){
      
      desmat_display <- as.data.frame( designMatrix(design = "balanced", periods = J_use(), clusters = (J_use()-1)) )
      colnames(desmat_display) <- paste("Period ", seq(J_use()))
      rownames(desmat_display) <- paste("Sequence ", seq((J_use() - 1)))
    
    }else{
      
      col_names <- paste("Period ", seq(J_use()))
      row_names <- paste("Cluster ", seq(n_use()))
      
      if(balanced() == "balanced"){
        desmat_display <- as.data.frame(desmat())#designMatrix(design = balanced(), periods = J_use(), clusters = n())
        
        colnames(desmat_display) <- col_names
        rownames(desmat_display) <- row_names
        
      }else if(balanced() == "unbalanced"){
        if(is.null(file1())) stop("User needs to upload design matrix before the function can continue")
        desmat_display <- as.data.frame(desmat())#designMatrix(design = balanced(), file=file1())
        
        colnames(desmat_display) <- col_names
        rownames(desmat_display) <- row_names
      }
      
    }
    desmat_display
    
  },digits=0, rownames=TRUE)
  
  
})
