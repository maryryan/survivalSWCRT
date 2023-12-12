#### LIBRARIES ####
library(tidyverse)
library(pracma)
library(survival)
library(stabledist)
library(latex2exp)
library(patchwork)

sim_scenarios_fixed_tauw <- read.table("gicc-tau-scenarios-fixed-tauw.txt",sep=",",header = T)
sim_scenarios_fixed_taub <- read.table("gicc-tau-scenarios-fixed-taub.txt",sep=",",header = T)

#### FUNCTIONS ####
CoxSurBCV <- function(Y,Delta,X,period,ID){
  
  # INPUT
  # Y: Observed time-to-event
  # Delta: Censoring indicator
  # X: Marginal mean covariates (design matrix excluding intercept)
  # period: periods each observation is occurring at
  # ID: Cluster identifier
  
  library(survival)
  library(pracma)
  
  ######################################
  # Step 1: prepare data elements
  ######################################
  # point estimate
  test.cox_cluster <- coxph(Surv(Y,Delta)~X+strata(factor(period))+cluster(ID))
  beta <- as.matrix(coef(test.cox_cluster))
  beta0 <- 0
  Ustar <- Ustar0 <- Ustar_m <- Ustar0_m <-  UUtran <- UUMD <- UUKC <- UUFG <- score_full <- score2_full <- score0_full <- score02_full <- 0
  
  for(p in 1:max(unique(period))){
    
    Y.p <- Y[period==p]
    X.p <- X[period==p]
    Delta.p <- Delta[period==p]
    ID.p <- ID[period==p]
    
    # sort observations by time - WTIHIN PERIOD
    b <- order(Y.p) #give indices of smallest to largest survival times
    Y.p <- sort(Y.p) #rearranges survival times in order of smallest to largest
    X.p <- as.matrix(X.p)[b,,drop=FALSE]
    ID.p <- ID.p[b]
    Delta.p <- Delta.p[b]
    ny <- length(Y.p)#number of observations in period
    nbeta <- dim(as.matrix(X.p))[2]# number of covariates
    UID <- sort(unique(ID.p))
    n <- length(UID)#number of unique IDs/clusters
    
    IDind <- zeros(length(UID), ny)#number of unique IDs by number of observations
    for (i in 1:length(UID)){
      IDind[i, ID.p==UID[i]] <- 1#mark observation as belonging to that unique ID
    }
    
    
    
    
    ######################################################
    # Step 2: Model-based variance estimates
    ######################################################
    # the rate of counting process of event, or dN(t)
    NN <- diag(Delta.p)
    # use the following trick to obtain the at-risk process Y(t)
    # each row is an individual
    # each column is a specific time point (recall the counting process notation)
    IndYY <- (t(repmat(t(Y.p),ny,1))>=repmat(t(Y.p),ny,1))
    Xbeta <- c(X.p%*%beta)
    Xbeta0 <- c(X.p%*%beta0)
    
    # Three S_j matrices for variance calculation
    S0beta <- colSums(IndYY*exp(Xbeta))#sum of Y*e^ZB over all individuals at each time
    S0beta0 <- colSums(IndYY*exp(Xbeta0))
    S1beta <- S1beta0 <- zeros(nbeta,ny)
    for(k in 1:nbeta){
      S1beta[k,] <- colSums(IndYY*exp(Xbeta)*X.p[,k])#sum of Z*Y*e^ZB over all individuals at each time
      S1beta0[k,] <- colSums(IndYY*exp(Xbeta0)*X.p[,k])
    }
    S2beta <- S2beta0 <- array(0, c(ny,1,nbeta,nbeta))
    for(k in 1:nbeta){
      for(s in 1:nbeta){
        S2beta[,,k,s] <- colSums(IndYY*exp(Xbeta)*X.p[,k]*X.p[,s])#sum of Z^2*Y*e^ZB over all individuals at each time - will be same as S1beta if only have treatment indicator
        S2beta0[,,k,s] <- colSums(IndYY*exp(Xbeta0)*X.p[,k]*X.p[,s])
      }
    }
    
    Omega <- Omega0 <- array(0,c(nbeta,nbeta,n))
    score0 <- zeros(nbeta,ny)
    
    # obtain cluster-specific matrices
    for(i in UID){
      # subset observations from each cluster
      S0beta_c <- S0beta[IDind[i,]==1]
      S1beta_c <- S1beta[,IDind[i,]==1,drop=FALSE]
      Delta_c <- Delta.p[IDind[i,]==1]
      
      S0beta0_c <- S0beta0[IDind[i,]==1]
      S1beta0_c <- S1beta0[,IDind[i,]==1,drop=FALSE]
      
      # components for A matrix and score
      for (k in 1:nbeta){
        score0[k,] <- sum(Delta_c*(X.p[,k]-S1beta0_c[k,]/S0beta0_c))
        for (s in 1:nbeta){
          Omega[k,s,i] <- sum(Delta_c*(S2beta[IDind[i,]==1,,k,s]/S0beta_c-
                                         S1beta_c[k,]*S1beta_c[s,]/S0beta_c^2))#this is your A matrix
          Omega0[k,s,i] <- sum(Delta_c*(S2beta0[IDind[i,]==1,,k,s]/S0beta0_c-
                                          S1beta0_c[k,]*S1beta0_c[s,]/S0beta0_c^2))
        }
      }
    }
    
    score0_temp <- sum(score0)
    score02_temp <- score0_temp^2
    score0_full <- score0_full+score0_temp
    score02_full <- score02_full+score02_temp
    Ustar_temp <- apply(Omega,c(1,2),sum) #like Omega_j
    Ustar <- Ustar+Ustar_temp #summing over periods as we go
    Ustar0_temp <- apply(Omega0,c(1,2),sum) #like Omega_j
    Ustar0 <- Ustar0+Ustar0_temp #summing over periods as we go
  }#need to end period look here so we can find naive variance
  naive <- solve(Ustar)#not robust sandwich
  
  ######################################################
  # Step 3: Robust variance estimates and bias corrections
  ######################################################
  nom <- nom2 <- nom0 <- nom02 <- nomMD <- nomFG <- 0
  # start period loop again to find robust/bias-corrected variances, which use naive variance
  for(p in 1:max(unique(period))){
    Y.p <- Y[period==p]
    X.p <- X[period==p]
    Delta.p <- Delta[period==p]
    ID.p <- ID[period==p]
    
    # sort observations by time - WTIHIN PERIOD
    b <- order(Y.p) #give indices of smallest to largest survival times
    Y.p <- sort(Y.p) #rearranges survival times in order of smallest to largest
    X.p <- as.matrix(X.p)[b,,drop=FALSE]
    ID.p <- ID.p[b]
    Delta.p <- Delta.p[b]
    ny <- length(Y.p)#number of observations in period
    nbeta <- dim(as.matrix(X.p))[2]# number of covariates
    UID <- sort(unique(ID.p))
    n <- length(UID)#number of unique IDs/clusters
    
    IDind <- zeros(length(UID), ny)#number of unique IDs by number of observations
    for (i in 1:length(UID)){
      IDind[i, ID.p==UID[i]] <- 1#mark observation as belonging to that unique ID
    }
    
    
    
    
    ######################################################
    # Step 2: Model-based variance estimates
    ######################################################
    # the rate of counting process of event, or dN(t)
    NN <- diag(Delta.p)
    # use the following trick to obtain the at-risk process Y(t)
    # each row is an individual
    # each column is a specific time point (recall the counting process notation)
    IndYY <- (t(repmat(t(Y.p),ny,1))>=repmat(t(Y.p),ny,1))
    Xbeta <- c(X.p%*%beta)
    Xbeta0 <- c(X.p%*%beta0)
    
    # Three S matrices for variance calculation - I would need to do these by period (so will have J S0betas)
    S0beta <- colSums(IndYY*exp(Xbeta))#sum of Y*e^ZB over all individuals at each time
    S1beta <- S1beta0 <- zeros(nbeta,ny)
    S0beta0 <- colSums(IndYY*exp(Xbeta0))
    for(k in 1:nbeta){
      S1beta[k,] <- colSums(IndYY*exp(Xbeta)*X.p[,k])#sum of Z*Y*e^ZB over all individuals at each time
      S1beta0[k,] <- colSums(IndYY*exp(Xbeta0)*X.p[,k])
    }
    S2beta <- S2beta0 <- array(0, c(ny,1,nbeta,nbeta))
    for(k in 1:nbeta){
      for(s in 1:nbeta){
        S2beta[,,k,s] <- colSums(IndYY*exp(Xbeta)*X.p[,k]*X.p[,s])#sum of Z^2*Y*e^ZB over all individuals at each time - will be same as S1beta if only have treatment indicator
        S2beta0[,,k,s] <- colSums(IndYY*exp(Xbeta0)*X.p[,k]*X.p[,s])
      }
    }
    
    
    
    # start of new code #
    # Breslow estimator of baseline hazard
    dHY <- colSums(NN)/c(S0beta)
    HY <- cumsum(dHY)
    dHY0 <- colSums(NN)/c(S0beta0)
    HY0 <- cumsum(dHY0)
    
    # obtain martingale increment: 
    # recall that the martingale is the "residual" in survival context
    epsilon <- NN-IndYY*repmat(t(dHY),ny,1)*exp(Xbeta)#dM
    epsilon0 <- NN-IndYY*repmat(t(dHY0),ny,1)*exp(Xbeta0)#dM
    nom_temp <- nom0_temp <- nom2_temp <- nomMD_temp <- nomFG_temp <- zeros(nbeta,n)
    Omega_m <- Omega0_m <- array(0,c(nbeta,nbeta,n))
    
    # obtain cluster-specific matrices
    for(i in UID){
      # subset observations from each cluster
      X_c <- X.p[IDind[i,]==1,,drop=FALSE]
      epsilon_c <- epsilon[IDind[i,]==1,,drop=FALSE]
      epsilon_c_all <- colSums(epsilon_c)
      S0beta_c <- S0beta[IDind[i,]==1]
      S1beta_c <- S1beta[,IDind[i,]==1,drop=FALSE]
      ny_c <- sum(IDind[i,]==1)
      Delta_c <- Delta.p[IDind[i,]==1]
      IndYY_c <- IndYY[IDind[i,]==1,,drop=FALSE]
      Xbeta_c <- Xbeta[IDind[i,]==1]
      ylxb_c <- IndYY_c*exp(Xbeta_c)*repmat(dHY,ny_c,1)
      
      epsilon0_c <- epsilon0[IDind[i,]==1,,drop=FALSE]
      epsilon0_c_all <- colSums(epsilon0_c)
      S0beta0_c <- S0beta0[IDind[i,]==1]
      S1beta0_c <- S1beta0[,IDind[i,]==1,drop=FALSE]
      Xbeta0_c <- Xbeta0[IDind[i,]==1]
      ylxb0_c <- IndYY_c*exp(Xbeta0_c)*repmat(dHY0,ny_c,1)
      
      # the trick is to loop through the dimension of the coefficients
      # otherwise need to deal with multi-dimensional array, very complex
      for (k in 1:nbeta){
        # components for B matrix
        tempk <- repmat(X_c[,k,drop=FALSE],1,ny)-repmat(S1beta[k,,drop=FALSE]/t(S0beta),ny_c,1)#(Z-W)
        tempk0 <- repmat(X_c[,k,drop=FALSE],1,ny)-repmat(S1beta0[k,,drop=FALSE]/t(S0beta0),ny_c,1)#(Z-W)
        nom_temp[k,i] <- sum(as.matrix(tempk*epsilon_c)%*%repmat(1,ny,1))#(Z-W)dM
        nom0_temp[k,i] <- sum(as.matrix(tempk0*epsilon0_c)%*%repmat(1,ny,1))
        
        # preparation for residual based correction of B matrix - none of this concerns the naive
        for (s in 1:nbeta){
          # true Omega
          Omega_m[k,s,i] <- sum(Delta_c*(S2beta[IDind[i,]==1,,k,s]/S0beta_c-S1beta_c[k,]*S1beta_c[s,]/S0beta_c^2)) -
            sum((repmat(S2beta[,,k,s]/S0beta-S1beta[k,]*S1beta[s,]/S0beta^2,ny_c,1)*ylxb_c)%*%repmat(1,ny,1)) +
            sum((tempk*repmat(X_c[,s,drop=FALSE],1,ny)*ylxb_c)%*%repmat(1,ny,1))
          
          Omega0_m[k,s,i] <- sum(Delta_c*(S2beta0[IDind[i,]==1,,k,s]/S0beta0_c-S1beta0_c[k,]*S1beta0_c[s,]/S0beta0_c^2)) -
            sum((repmat(S2beta0[,,k,s]/S0beta0-S1beta0[k,]*S1beta0[s,]/S0beta0^2,ny_c,1)*ylxb0_c)%*%repmat(1,ny,1)) +
            sum((tempk0*repmat(X_c[,s,drop=FALSE],1,ny)*ylxb0_c)%*%repmat(1,ny,1))
        }
        
      }
      
      # components for MD type correction
      nomMD_temp[,i] <- Ustar%*%solve(Ustar-Omega_m[,,i])%*%nom_temp[,i]
      
      # components for FG type correction
      Hi <- zeros(nbeta,nbeta)
      tempov <- Omega_m[,,i]%*%naive
      diag(Hi) <- 1/sqrt(1-pmin(0.75,c(diag(tempov))))
      nomFG_temp[,i] <- Hi%*%nom_temp[,i]
    }
    
    # check model-based variance (should be same as naive)
    Ustar_m_temp <- apply(Omega_m,c(1,2),sum)
    Ustar_m <- Ustar_m+Ustar_m_temp
    Ustar0_m_temp <- apply(Omega0_m,c(1,2),sum)
    Ustar0_m <- Ustar0_m+Ustar0_m_temp
    
    nom <- nom + nom_temp #summing up results from each period
    nom0 <- nom0 + nom0_temp
    nom2_temp <- tcrossprod(nom_temp)
    nom02_temp <- tcrossprod(nom0_temp)
    nom2 <- nom2 + nom2_temp
    nom02 <- nom02 + nom02_temp
    
    nomMD <- nomMD + nomMD_temp
    
    nomFG <- nomFG + nomFG_temp
  }
  #naive <- solve(Ustar)
  naive_m <- solve(Ustar_m)
  naive0_m <- solve(Ustar0_m)
  
  # robust variance estimator
  UUtran <- tcrossprod(nom) #this is B
  robust <- naive%*%UUtran%*%naive
  
  UUtran0 <- tcrossprod(nom0) #this is B
  
  # variance estimator of MD type
  UUMD <- tcrossprod(nomMD)
  varMD <- naive%*%UUMD%*%naive
  
  # variance estimator of KC type
  UUKC <- tcrossprod(nomMD,nom)
  varKC <- naive%*%(UUKC+t(UUKC))%*%naive/2
  
  # variance estimator of FG type
  UUFG <- tcrossprod(nomFG)
  varFG <- naive%*%UUFG%*%naive
  
  #############################################
  # Output
  # naive: naive or model-based var
  # robust: robust sandwich var
  # varRB: bias-corrected sandwich var due to residual based correction
  # varMD: bias-corrected sandwich var due to MD type
  # varKC: bias-corrected sandwich var due to KC type
  # varFG: bias-corrected sandwich var due to FG type
  #############################################
  bSE <- sqrt(diag(naive))
  bSEBC0 <- sqrt(diag(robust))
  bSEMD <- sqrt(diag(varMD))
  bSEKC <- sqrt(diag(varKC))
  bSEFG <- sqrt(diag(varFG))
  
  outbeta <- cbind(summary(test.cox_cluster)$coefficients,
                   bSE,bSEBC0,#bSERB,
                   bSEMD,bSEKC,bSEFG,
                   score0_full,sqrt(UUtran0),
                   sum(nom0),score02_full,sum(nom02),
                   UUMD,UUKC,UUFG#,bSEMBN,
                   #bSERBMD,bSERBKC,bSERBFG,bSERBMBN
  )
  colnames(outbeta) <- c("coef","exp(coef)","se(coef)","robust se","z","Pr(>|z|)",
                         "MB-stderr","BC0-stderr",
                         "MD-stderr","KC-stderr","FG-stderr",
                         "score0", "BC0-B0",
                         "nom0","score02","nom02",
                         "MD-B","KC-B","FG-B"
  )
  
  return(list(outbeta=outbeta))
}

surGENSW<-function(n, mv, lambda0, beta, tau_b, tau_w, p.max, k){
  
  # INPUT
  # n: Number of clusters
  # mv: Vector (length n) of numbers of observations in each cluster
  # lambda0: baseline hazard function - if this is a vector you can create time-dependent baseline hazards
  # beta: Intervention effect
  # tau_b: Kendall's tau between clusters within a period (must be <= tau_w)
  # tau_w: Kendall's tau between subjects within a cluster-period (must be >= tau_b)
  # k: the sequence each cluster belongs to/time period each cluster begins treatment (vector)
  # p.max: maximum number of periods observed
  
  y <- NULL
  
  # inverse of marginal survival function: inverse Exponential distribution
  t <- function(lambda0, zbeta, S){
    t_S <- (1/lambda0)*(-(log(S))*exp(-zbeta))
    return(t_S)
  }
  
  # function to generate y for one cluster with m observations (this will be all observations across all time periods)
  rnested_gumbel <- function(n, p.max, m, tau_b, tau_w, lambda_ij){
    
    ## A function to sample from a nested Gumbel function
    # n: number of clusters
    # p.max: number of time periods
    # m: number of subjects per cluster-period; can be scalar if constant or a vector if changes by cluster/period
    # tau_b: between-period kendall's tau
    # tau_w: within-period kendall's tau
    # lambda_ij: cluster-period-specific hazard function; this is a matrix where rows are cluster/sequence and columns are period
    
    require(stabledist)
    
    T_ijk <- NULL#matrix(NA, nrow=n, ncol=J*m)
    
    theta_b <- 1-tau_b
    theta_w <- 1-tau_w
    # theta0: parent-level correlation parameter; must be >= 1
    theta0 <- 1/theta_b
    # theta1: child-level correlation parameter; must be > theta0
    theta1 <- theta_b/theta_w
    
    stable0_gamma <- ( cos(pi/(theta0*2)) )
    stable01_gamma <- ( cos(pi/(theta1*2)) )
    
    for(i in seq(n)){
      V0 <- V01 <- Z <- U <- T_ijk_clu <- NULL
      
      # step 1 - pm=1 is the parameterization that makes this positive stable
      V0 <- rstable(1, 1/theta0, 1, stable0_gamma, 0, pm=1)
      
      # step 2
      V01 <- rstable(p.max, 1/theta1, 1, stable01_gamma, 0, pm=1)
      
      # step 3
      Z <- vector(mode="list", length=p.max)
      if( is.null(dim(m)) ){
        
        if(length(m)==1){
          # if there's a constant cluster-period size across clusters & periods #
          for(p in seq(p.max)){
            Z[[p]] <- runif(m, 0, 1)
            
          }
          
        }else if( length(m)==n ){
          # if cluster-period size differs by cluster but not period #
          for(p in seq(p.max)){
            Z[[p]] <- runif(m[i], 0, 1)
            
          }
          
        }else if( length(m) == p.max ){
          # if cluster-period size differs by period but not cluster #
          for(p in seq(p.max)){
            Z[[p]] <- runif(m[p], 0, 1)
          }
          
        }
        
      }else{
        # cluster-period size differs by cluster and period #
        for(p in seq(p.max)){
          Z[[p]] <- runif(sum(m[i,p]), 0, 1)
          
        }
      }
      
      
      # step 4
      U <- vector(mode="list", length=p.max)
      for(p in seq(p.max)){
        U[[p]] <- exp( -( -log(Z[[p]])/V01[p] )^(1/theta1) )
      }
      
      
      # step 5
      if( is.null(dim(lambda_ij)) & n==1 ){
        
        
        if(length(m)==1){
          # if there's a constant cluster-period size across clusters & periods #
          for(p in seq(p.max)){
            T_ijk_clu_temp <- (1/rep(lambda_ij[p], m)) * ( -log(U[[p]])/V0 )^(1/theta0)
            T_ijk_clu <- c(T_ijk_clu, T_ijk_clu_temp)# stacking periods for same cluster on bottom
          }
          
        }else if( length(m) == p.max ){
          # if cluster-period size differs by period but not cluster #
          for(j in seq(p.max)){
            T_ijk_clu_temp <- (1/rep(lambda_ij[j], m[j])) * ( -log(U[[j]])/V0 )^(1/theta0)
            T_ijk_clu <- c(T_ijk_clu, T_ijk_clu_temp)
          }
          
        }
        
        
      }else{
        #lambda_ij is a matrix
        if( is.null(dim(m)) ){
          
          if(length(m)==1){
            # if there's a constant cluster-period size across clusters & periods #
            for(p in seq(p.max)){
              T_ijk_clu_temp <- (1/rep(lambda_ij[i,p], m)) * ( -log(U[[p]])/V0 )^(1/theta0)
              T_ijk_clu <- c(T_ijk_clu, T_ijk_clu_temp)
              
            }
            
          }else if( length(m)==n ){
            # if cluster-period size differs by cluster but not period #
            for(p in seq(p.max)){
              
              T_ijk_clu_temp <- (1/rep(lambda_ij[i,p], m[i])) * ( -log(U[[p]])/V0 )^(1/theta0)
              T_ijk_clu <- c(T_ijk_clu, T_ijk_clu_temp)
              
            }
            
          }else if( length(m) == p.max ){
            # if cluster-period size differs by period but not cluster #
            for(j in seq(p.max)){
              T_ijk_clu_temp <- (1/rep(lambda_ij[i,j], m[j])) * ( -log(U[[j]])/V0 )^(1/theta0)
              T_ijk_clu <- c(T_ijk_clu, T_ijk_clu_temp)
            }
            
          }
          
        }else{
          # cluster-period size differs by cluster and period #
          for(j in seq(p.max)){
            T_ijk_clu_temp <- (1/rep(lambda_ij[i,j], m[i,j])) * ( -log(U[[p]])/V0 )^(1/theta0)
            T_ijk_clu <- c(T_ijk_clu, T_ijk_clu_temp)
          }
          
        }
      }
      
      
      T_ijk <- c(T_ijk, T_ijk_clu)
    }# end n loop
    
    # T_ijk organized as all subjects in cluster 1, period 1, then cluster 1, period 2, etc.
    return(T_ijk)
    
  }
  
  
  
  # generate y for n clusters
  for (i in 1:n){
    # create treatment effect sequence for entire cluster, across time periods #
    zbeta <- c(0,beta)
    lambda_ij <- rep(NA, p.max)
    for(p in seq(p.max)){
      # if we're in the treatment period trt effect is beta[2], otherwise its beta[1] (null)
      zbeta.p <- ifelse(p >= k[i], zbeta[2], zbeta[1])
      lambda_ij[p] <- lambda0[p]*exp(zbeta.p)
    }
    
    yi <- rnested_gumbel(n=1, p.max=p.max, m=mv[i,], tau_b=tau_b, tau_w=tau_w, lambda_ij=lambda_ij)
    y <- c(y, yi) 
  }
  
  return(y)
}

surSIMULATESW <- function(n, J, cv.c, cv.p, lambda0, beta, pa, p0, tau_b, tau_w, Cp, nrep, k, p.max, strat=TRUE, FE=FALSE){
  
  # INPUT
  # n: Number of clusters
  # J: Mean cluster size
  # cv.c: CV of cluster size across clusters
  # cv.p: CV of cluster size across time periods
  # lambda0: baseline hazard function(s) - if this is a vector, you can create period-stratified baseline hazards
  # beta: Intervention effect
  # pa: Desired administrative censoring rate for the control group
  # p0: Net censoring rate in the control arm
  # tau_b: Kendall's tau between clusters within a period (must be <= tau_w)
  # tau_w: Kendall's tau between subjects within a cluster-period (must be >= tau_b)
  # Cp: Administrative censoring time
  # nrep: Number of replications
  # k: the sequence each cluster belongs to/time period each cluster begins treatment (vector)
  # p.max: maximum number of observed periods
  # strat: do you want to fit a stratified cox model on the periods? default is TRUE
  # FE: do you want to fit fixed effects for period in the cox model? default is FALSE
  
  require(survival)
  
  results <- NULL
  
  mv <- matrix(NA, ncol=p.max,nrow=n)
  for (i in 1:nrep){
    # Generate cluster size
    # if cv.c is 0, same cluster size for all clusters and periods
    if (cv.c == 0){
      for(p in 1:p.max){
        
        mv[,p] <- rep(J, n)
        
      }
      
    } else{
      # if cv.c != 0 generate "general" random cluster sizes
      gamma_a <- cv.c^(-2)
      gamma_b <- 1/(J*cv.c^2)
      set.seed(i)
      mv1 <- round(rgamma(n, shape=gamma_a, rate=gamma_b))
      
      # if cv.p is 0, each cluster gets same size across periods
      if(cv.p == 0){
        
        for(p in 1:p.max){
          mv[,p] <- mv1
          mv[mv <= 2,p] <- 2
          
        }
        
      } else{
        # if cv.p != 0 jitter cluster's "general" random size for each period
        mv[,p] <- round(rmvnorm(n, mean=mv1, sigma=diag(cv.p,n)))
        mv[mv <= 2,p] <- 2
        
      }
      
    }
    
    # Create id and Z matrix
    # cluster id #
    id <- rep(1:n, rowSums(mv))#1:n and mv need to be the same length
    # period id for each cluster #
    period <- NULL
    for(j in 1:n){
      for(p in 1:p.max){
        period1 <- rep(p,each=mv[j,p])
        period <- c(period, period1)
      }
      
    }
    
    mt <- sum(mv)# total observations
    Z <- NULL
    for(j in 1:n){
      # treatment indicator for cluster across time periods #
      z1<- c(rep(0,(k[j]-1)), rep(1,(p.max-k[j]+1)))
      
      for(p in 1:p.max){
        # replicate for each subject in each time period #
        z2 <- rep(z1[p],each=mv[j,p])  
        Z <- c(Z,z2)
      }
      
      
    }
    
    # Generate failure times
    set.seed(i)
    y <- surGENSW(n, mv, lambda0, beta, tau_b, tau_w, p.max, k)
    
    # Generate censoring times
    # if admin and net censoring rates are equal #
    if (p0 == pa){
      c <- rep(1, mt)
    } else{
      # if admin and net censoring rates are not equal #
      set.seed(i)
      c <- runif(mt)
    }
    
    # Generate censoring indicator and censored event times
    status <- rep(1, mt)
    # if failure greater than censoring, status == 2 (uniform censoring)
    status[y > c] <- 2
    # if failure greater than 1, status == 0 (admin censoring) #
    status[y > 1] <- 0
    # if y < 1 & < c, status == 1 #
    
    times <- y
    # if failure time greater than censoring time, fill censoring time #
    times[status == 2] <- c[status == 2]
    # if failure time greater than 1, fill 1 (admin censoring time?) #
    times[status == 0] <- 1
    # if y < 1 & < c, time is actual failure time #
    # give those with non-admin censoring status the same status as admin censoring #
    status[status == 2] <- 0
    
    ctrl <- status[Z == 0]
    
    # Cox model
    if(strat==TRUE){
      # won't fit FEs for period if stratifying on period
      if(FE==TRUE){
        return("Pick either strat or FE")
      }else{
        survival <- coxph(Surv(times, status) ~ (Z) + strata(factor(period)) , cluster = id)  
      }
      
    }else{
      if(FE==TRUE){
        survival <- coxph(Surv(times, status) ~ (Z) + factor(period) , cluster = id)  
      }else{
        survival <- coxph(Surv(times, status) ~ (Z) , cluster = id)  
      }
      
    }
    
    BC_results <- CoxSurBCV(times,status,Z,period,id)
    
    results <- rbind(results, c(summary(survival)$coefficients, n, J, mean(mv), beta, tau_b, tau_w, 1-mean(ctrl), BC_results$outbeta))
    if( i %in% c(nrep*seq(0.1,1,0.1)) ) print(paste0(Sys.time(), ": ", (i/nrep)*100, "% done"))
    
  }#end nrep
  colnames(results) <- c("coef", "exp(coef)", "se(coef)", "robust se",
                         "z", "Pr(>|z|)",
                         "n", "J", "empirical J", "beta", "tau_b", "tau_w", "empirical p0",
                         "coef","exp(coef)","se(coef)","robust se","z","Pr(>|z|)",
                         "MB-stderr","BC0-stderr",
                         "MD-stderr","KC-stderr","FG-stderr",
                         "score0", "BC0-B0",
                         "nom0","score02","nom02",
                         "MD-B","KC-B","FG-B"
  )
  
  if (beta == 0){
    name <- "Size"
  } else{
    name <- "Power"}
  cv.c_name <- cv.c*10
  return(results)
}

q_0 <- function(s, beta, z_ij, j,lambda0, tau,pi_b){
  # censoring distribution
  G_s <- 1-s/tau
  
  # W_js (ratio of s_1j/s_0j)
  pi_b_j <- sum(pi_b[1:j])
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

q_1 <- function(s, t, beta, z_ij, z_il, j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  
  # transform kendall's tau into correlation parameters
  theta_b <- 1-tau_kb
  theta_w <- 1-tau_kw
  
  # censoring distribution - independent
  G_st <- (1-s/tau)*(1-t/tau)
  
  # Ws
  pi_b_j <- sum(pi_b[1:j])
  s_1js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s)
  s_0js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
  W_js <- s_1js/s_0js
  
  pi_b_l <- sum(pi_b[1:l])
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

q_2 <- function(s, t, beta, z_ij, z_il, j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  
  # transform kendall's tau into correlation parameters
  theta_b <- 1-tau_kb
  theta_w <- 1-tau_kw
  
  # censoring distribution
  G_st <- (1-s/tau)*(1-t/tau)
  
  # Ws
  pi_b_j <- sum(pi_b[1:j])
  s_1js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s)
  s_0js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
  W_js <- s_1js/s_0js
  
  pi_b_l <- sum(pi_b[1:l])
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

q_3 <- function(s, t, beta, z_ij, z_il, j ,l, lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  
  # transform kendall's tau into correlation parameters
  theta_b <- 1-tau_kb
  theta_w <- 1-tau_kw
  
  # censoring distribution
  G_st <- (1-s/tau)*(1-t/tau)
  
  # Ws
  pi_b_j <- sum(pi_b[1:j])
  s_1js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s)
  s_0js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
  W_js <- s_1js/s_0js
  
  pi_b_l <- sum(pi_b[1:l])
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

q_4 <- function(s,t, beta, z_ij, z_il, j, l, lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  # transform kendall's tau into correlation parameters
  theta_b <- 1-tau_kb
  theta_w <- 1-tau_kw
  
  # censoring distribution
  G_st <- (1-s/tau)*(1-t/tau)
  
  # Ws
  pi_b_j <- sum(pi_b[1:j])
  s_1js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s)
  s_0js <- pi_b_j * exp(beta - exp(beta)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
  W_js <- s_1js/s_0js
  
  pi_b_l <- sum(pi_b[1:l])
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
V_0 <- function(s, beta, z_ij, j, lambda0, tau, pi_b){
  
  # censoring distribution
  G_s <- 1-s/tau
  
  # W_js (ratio of s_1j/s_0j)
  pi_b_j <- sum(pi_b[1:j])
  s_1j <- s_2j <- pi_b_j * exp(beta - exp(beta)*lambda0*s)
  s_0j <- pi_b_j * exp(beta - exp(beta)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
  
  # Density functions (negative derivative of the survival function)
  f_ij_s <- lambda0*exp(z_ij*beta)*exp(-lambda0*s*exp(z_ij*beta))
  
  v <- G_s*((s_2j/s_0j) - (s_1j^2)/(s_0j^2))*f_ij_s
  
  return(v)
  
}


Q <- function(beta, z_ij, z_il, j, l, lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  # function that integrates and adds together the 4 q components for the non-same subject portion of B
  Q <- dblquad(q_1,0, tau,0,tau,dim=2, beta=beta, z_ij=z_ij, z_il=z_il, j=j, l=l, lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb) -
    dblquad(q_2,0,tau,0,tau, dim=2, beta=beta, z_ij=z_ij, z_il=z_il, j=j, l=l,lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb) -
    dblquad(q_3,0,tau,0,tau,dim=2, beta=beta, z_ij=z_ij, z_il=z_il, j=j, l=l, lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb) +
    dblquad(q_4,0,tau,0,tau, dim=2, beta=beta, z_ij=z_ij, z_il=z_il, j=j, l=l,lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb)
  
  return(Q)
}

H0 <- function(j, beta, lambda0, tau,pi_b){
  # function that integrates and takes the expectation for the same subject, cluster, period portion of B
  pi_b_j <- sum(pi_b[1:j])
  
  H0 <- pi_b_j*pracma::integral(fun=q_0,xmin=0,xmax=tau,beta=beta, z_ij=1,j=j,lambda0=lambda0, tau=tau,pi_b=pi_b) +
    (1-pi_b_j)*pracma::integral(fun=q_0,xmin=0,xmax=tau, beta=beta, z_ij=0,j=j,lambda0=lambda0, tau=tau,pi_b=pi_b)
  
  return(H0)
}

H1 <- function(j, l, beta, lambda0_j, lambda0_l, tau,pi_b,tau_kw, tau_kb){
  # function that takes the expectation and adds together the integrated q components for the non-same subject portion of B
  pi_b_min <- ifelse( j < l, sum(pi_b[1:j]), sum(pi_b[1:l]) )
  pi_b_max <- ifelse( j > l, sum(pi_b[1:j]), sum(pi_b[1:l]) )
  if( j > l){
    pi_b_mid <- sum(pi_b[(l+1):j])
  }else if(j < l){
    pi_b_mid <- sum(pi_b[(j+1):l]) 
  }else if(j == l){
    pi_b_mid <- 0
  }
  
  H1 <-  (1-pi_b_max)*Q(beta, 0,0,j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb)
  
  if(j > 1 | l > 1){
    
    if(j > l){
      H1 <- H1 + pi_b_mid*Q(beta, 1,0,j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb)
    }else if(j < l){
      H1 <- H1 + pi_b_mid*Q(beta, 0,1,j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb)
    }
    
    if(j > 1 & l > 1){
      H1 <- H1 + pi_b_min*Q(beta, 1,1,j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb)
    }
  }
  
  return(H1)
}

H2 <- function(j, beta, lambda0, tau, pi_b){
  # function that integrates and takes the expectation of the A components
  
  pi_b_j <- sum(pi_b[1:j])
  
  H2 <- pi_b_j*pracma::integral(fun=V_0, xmin=0, xmax=tau, beta=beta, z_ij=1, j=j, lambda0=lambda0, tau=tau, pi_b=pi_b) + 
    (1-pi_b_j)*pracma::integral(fun=V_0, xmin=0, xmax=tau, beta=beta, z_ij=0, j=j, lambda0=lambda0, tau=tau, pi_b=pi_b)
  
  return(H2)
}

sandwich_var <- function(m,J, lambda0,tau,pi_b, tau_kw, tau_kb, beta=beta){
  H0_sum <- H2_sum <- rep(0, J)
  #H_score_sum <- rep(NA,J)
  H1_sum_eq <- H1_sum_uneq <- H0_sum_iccb <- matrix(NA, nrow=J, ncol=J)
  for(j in seq(J)){
    
    H0_sum[j] <- H0(j=j, beta=beta, lambda0=lambda0[j], tau=tau, pi_b=pi_b)
    H2_sum[j] <- H2(j=j, beta=beta, lambda0=lambda0[j], tau=tau, pi_b=pi_b)
    
    for(l in seq(J)){
      if(j==l){
        H1_sum_eq[j,l] <- H1(j=j,l=l, beta=beta, lambda0_j=lambda0[j], lambda0_l=lambda0[l], tau=tau, pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb)
      }else{
        H1_sum_uneq[j,l] <- H1(j=j,l=l, beta=beta, lambda0_j=lambda0[j], lambda0_l=lambda0[l], tau=tau, pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb)
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

#### SCORE FUNCTIONS ####
q_0_score <- function(s, beta0, betaA, z_ij, j,lambda0, tau,pi_b){
  # censoring distribution
  G_s <- 1-s/tau
  
  # W_js (ratio of s_1j/s_0j)
  pi_b_j <- sum(pi_b[1:j])
  s_1j <- pi_b_j * exp(beta0 - exp(betaA)*lambda0*s)
  s_0j <- pi_b_j * exp(beta0 - exp(betaA)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
  W_js <- s_1j/s_0j
  
  # Density functions (negative derivative of the survival function)
  # density on treatment arm # 
  f_ij_s <- lambda0*exp(z_ij*betaA)*exp(-lambda0*s*exp(z_ij*betaA))
  
  # final output
  q0 <- G_s*(z_ij-W_js)^2*f_ij_s
  
  return(q0)
}

q_1_score <- function(s, t, beta0, betaA, z_ij, z_il, j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  
  # transform kendall's tau into correlation parameters
  theta_b <- 1-tau_kb
  theta_w <- 1-tau_kw
  
  # censoring distribution - independent
  G_st <- (1-s/tau)*(1-t/tau)
  
  # Ws
  pi_b_j <- sum(pi_b[1:j])
  s_1js <- pi_b_j * exp(beta0 - exp(betaA)*lambda0_j*s)
  s_0js <- pi_b_j * exp(beta0 - exp(betaA)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
  W_js <- s_1js/s_0js
  
  pi_b_l <- sum(pi_b[1:l])
  s_1lt <- pi_b_l * exp(beta0 - exp(betaA)*lambda0_l*t)
  s_0lt <- pi_b_l * exp(beta0 - exp(betaA)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
  W_lt <- s_1lt/s_0lt
  
  # Density functions
  # these will be the same if j == l 
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

q_2_score <- function(s, t, beta0, betaA, z_ij, z_il, j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  
  # transform kendall's tau into correlation parameters
  theta_b <- 1-tau_kb
  theta_w <- 1-tau_kw
  
  # censoring distribution
  G_st <- (1-s/tau)*(1-t/tau)
  
  # Ws
  pi_b_j <- sum(pi_b[1:j])
  s_1js <- pi_b_j * exp(beta0 - exp(betaA)*lambda0_j*s)
  s_0js <- pi_b_j * exp(beta0 - exp(betaA)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
  W_js <- s_1js/s_0js
  
  pi_b_l <- sum(pi_b[1:l])
  s_1lt <- pi_b_l * exp(beta0 - exp(betaA)*lambda0_l*t)
  s_0lt <- pi_b_l * exp(beta0 - exp(betaA)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
  W_lt <- s_1lt/s_0lt
  
  # Density functions
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

q_3_score <- function(s, t, beta0, betaA, z_ij, z_il, j ,l, lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  
  # transform kendall's tau into correlation parameters
  theta_b <- 1-tau_kb
  theta_w <- 1-tau_kw
  
  # censoring distribution
  G_st <- (1-s/tau)*(1-t/tau)
  
  # Ws
  pi_b_j <- sum(pi_b[1:j])
  s_1js <- pi_b_j * exp(beta0 - exp(betaA)*lambda0_j*s)
  s_0js <- pi_b_j * exp(beta0 - exp(betaA)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
  W_js <- s_1js/s_0js
  
  pi_b_l <- sum(pi_b[1:l])
  s_1lt <- pi_b_l * exp(beta0 - exp(betaA)*lambda0_l*t)
  s_0lt <- pi_b_l * exp(beta0 - exp(betaA)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
  W_lt <- s_1lt/s_0lt
  
  # Density functions
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

q_4_score <- function(s,t, beta0, betaA, z_ij, z_il, j, l, lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  # transform kendall's tau into correlation parameters
  theta_b <- 1-tau_kb
  theta_w <- 1-tau_kw
  
  # censoring distribution
  G_st <- (1-s/tau)*(1-t/tau)
  
  # Ws
  pi_b_j <- sum(pi_b[1:j])
  s_1js <- pi_b_j * exp(beta0 - exp(betaA)*lambda0_j*s)
  s_0js <- pi_b_j * exp(beta0 - exp(betaA)*lambda0_j*s) + (1-pi_b_j)*exp(-lambda0_j*s)
  W_js <- s_1js/s_0js
  
  pi_b_l <- sum(pi_b[1:l])
  s_1lt <- pi_b_l * exp(beta0 - exp(betaA)*lambda0_l*t)
  s_0lt <- pi_b_l * exp(beta0 - exp(betaA)*lambda0_l*t) + (1-pi_b_l)*exp(-lambda0_l*t)
  W_lt <- s_1lt/s_0lt
  
  # Density functions
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

V_0_score <- function(s, beta0, betaA, z_ij, j, lambda0, tau, pi_b){
  
  # censoring distribution
  G_s <- 1-s/tau
  
  # W_js (ratio of s_1j/s_0j)
  pi_b_j <- sum(pi_b[1:j])
  s_1j <- s_2j <- pi_b_j * exp(beta0 - exp(betaA)*lambda0*s)
  s_0j <- pi_b_j * exp(beta0 - exp(betaA)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
  
  # Density functions (negative derivative of the survival function)
  f_ij_s <- lambda0*exp(z_ij*betaA)*exp(-lambda0*s*exp(z_ij*betaA))
  
  v <- G_s*((s_2j/s_0j) - (s_1j^2)/(s_0j^2))*f_ij_s#(W_js - W_js2)*f_ij_s
  
  return(v)
  
}

Q_score <- function(beta0, betaA, z_ij, z_il, j, l, lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb){
  # function that integrates and adds together the 4 q components for the non-same subject portion of B
  Q <- dblquad(q_1_score,0, tau,0,tau,dim=2, beta0=beta0, betaA=betaA, z_ij=z_ij, z_il=z_il, j=j, l=l, lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb) -
    dblquad(q_2_score,0,tau,0,tau, dim=2,beta0=beta0, betaA=betaA, z_ij=z_ij, z_il=z_il, j=j, l=l,lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb) -
    dblquad(q_3_score,0,tau,0,tau,dim=2, beta0=beta0, betaA=betaA, z_ij=z_ij, z_il=z_il, j=j, l=l, lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb) +
    dblquad(q_4_score,0,tau,0,tau, dim=2, beta0=beta0, betaA=betaA, z_ij=z_ij, z_il=z_il, j=j, l=l,lambda0_j=lambda0_j, lambda0_l=lambda0_l, tau=tau,pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb)
  
  return(Q)
}

H0_score <- function(j, beta0, betaA, lambda0, tau,pi_b){
  # function that integrates and takes the expectation for the same subject, cluster, period portion of B
  pi_b_j <- sum(pi_b[1:j])
  
  H0 <- pi_b_j*pracma::integral(fun=q_0_score,xmin=0,xmax=tau,beta0=beta0, betaA=betaA, z_ij=1,j=j,lambda0=lambda0, tau=tau,pi_b=pi_b) +
    (1-pi_b_j)*pracma::integral(fun=q_0_score,xmin=0,xmax=tau, beta0=beta0, betaA=betaA, z_ij=0,j=j,lambda0=lambda0, tau=tau,pi_b=pi_b)
  
  return(H0)
}

H1_score <- function(j, l, beta0, betaA, lambda0_j, lambda0_l, tau,pi_b,tau_kw, tau_kb){
  # function that takes the expectation and adds together the integrated q components for the non-same subject portion of B
  pi_b_min <- ifelse( j < l, sum(pi_b[1:j]), sum(pi_b[1:l]) )
  pi_b_max <- ifelse( j > l, sum(pi_b[1:j]), sum(pi_b[1:l]) )
  if( j > l){
    pi_b_mid <- sum(pi_b[(l+1):j])
  }else if(j < l){
    pi_b_mid <- sum(pi_b[(j+1):l]) 
  }else if(j == l){
    pi_b_mid <- 0
  }
  
  H1 <-  (1-pi_b_max)*Q_score(beta0=beta0, betaA=betaA, 0,0,j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb)
  
  if(j > 1 | l > 1){
    
    if(j > l){
      H1 <- H1 + pi_b_mid*Q_score(beta0=beta0, betaA=betaA, 1,0,j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb)
    }else if(j < l){
      H1 <- H1 + pi_b_mid*Q_score(beta0=beta0, betaA=betaA, 0,1,j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb)
    }
    
    if(j > 1 & l > 1){
      H1 <- H1 + pi_b_min*Q_score(beta0=beta0, betaA=betaA, 1,1,j,l,lambda0_j, lambda0_l, tau,pi_b, tau_kw, tau_kb)
    }
  }
  
  return(H1)
}
H2_score <- function(j, beta0, betaA, lambda0, tau, pi_b){
  # function that integrates and takes the expectation of the A components
  
  pi_b_j <- sum(pi_b[1:j])
  
  H2 <- pi_b_j*pracma::integral(fun=V_0_score, xmin=0, xmax=tau, beta0=beta0, betaA=betaA, z_ij=1, j=j, lambda0=lambda0, tau=tau, pi_b=pi_b) + 
    (1-pi_b_j)*pracma::integral(fun=V_0_score, xmin=0, xmax=tau, beta0=beta0, betaA=betaA, z_ij=0, j=j, lambda0=lambda0, tau=tau, pi_b=pi_b)
  
  return(H2)
}

score_exp <- function(s, beta0, betaA, z_ij, j,lambda0, tau,pi_b){
  # censoring distribution
  G_s <- 1-s/tau
  
  # W_js (ratio of s_1j/s_0j)
  pi_b_j <- sum(pi_b[1:j])
  s_1j <- pi_b_j * exp(beta0 - exp(betaA)*lambda0*s)
  s_0j <- pi_b_j * exp(beta0 - exp(betaA)*lambda0*s) + (1-pi_b_j)*exp(-lambda0*s)
  W_js <- s_1j/s_0j
  
  # Density functions (negative derivative of the survival function)
  # density on treatment arm # 
  f_ij_s <- lambda0*exp(z_ij*betaA)*exp(-lambda0*s*exp(z_ij*betaA))
  
  
  # final output
  score_exp <- (G_s-f_ij_s)*(z_ij-W_js)
  return(score_exp)
}

H_score <- function(j, beta0, betaA, lambda0, tau,pi_b){
  # function that integrates and takes the expectation for the same subject, cluster, period portion of B
  pi_b_j <- sum(pi_b[1:j])
  
  H_score <- pi_b_j*pracma::integral(fun=score_exp,xmin=0,xmax=tau,beta0=beta0, betaA=betaA, z_ij=1,j=j,lambda0=lambda0, tau=tau,pi_b=pi_b) +
    (1-pi_b_j)*pracma::integral(fun=score_exp,xmin=0,xmax=tau, beta0=beta0, betaA=betaA, z_ij=0,j=j,lambda0=lambda0, tau=tau,pi_b=pi_b)
  
  return(H_score)
}

sandwich_var_score <- function(m,J, lambda0,tau,pi_b, tau_kw, tau_kb, beta0=beta0, betaA=betaA){
  H0_sum <- H2_sum <- rep(0, J)
  H_score_sum <- rep(NA,J)
  H1_sum_eq <- H1_sum_uneq <- matrix(NA, nrow=J, ncol=J)
  for(j in seq(J)){
    
    H0_sum[j] <- H0_score(j=j, beta0=beta0,betaA=betaA, lambda0=lambda0[j], tau=tau, pi_b=pi_b)
    H2_sum[j] <- H2_score(j=j, beta0=beta0,betaA=betaA, lambda0=lambda0[j], tau=tau, pi_b=pi_b)
    H_score_sum[j] <- H_score(j=j, beta0=beta0, betaA=betaA, lambda0=lambda0[j], tau=tau, pi_b=pi_b)
    for(l in seq(J)){
      if(j==l){
        H1_sum_eq[j,l] <- H1_score(j=j,l=l, beta0=beta0,betaA=betaA, lambda0_j=lambda0[j], lambda0_l=lambda0[l], tau=tau, pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb)
      }else{
        H1_sum_uneq[j,l] <- H1_score(j=j,l=l, beta0=beta0,betaA=betaA, lambda0_j=lambda0[j], lambda0_l=lambda0[l], tau=tau, pi_b=pi_b, tau_kw=tau_kw, tau_kb=tau_kb)
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




#### Fixed tau_w simulation scenarios ####
m <- as.numeric(sim_scenarios_fixed_tauw[,1])
n <- as.numeric(sim_scenarios_fixed_tauw[,2])
J <- as.numeric(sim_scenarios_fixed_tauw[,3])
tau_w <- as.numeric(sim_scenarios_fixed_tauw[,4])
tau_b <- as.numeric(sim_scenarios_fixed_tauw[,5])
beta0 <- 0
betaA <- 0.4
Cp <- 1
baseline_constant <- as.numeric(sim_scenarios_fixed_tauw[,6])

gICC_results_fixed_tauw <- matrix(NA, nrow=length(m), ncol=8)
colnames(gICC_results) <- c("n","m","J","baseline_constant","tau_w", "tau_b", "gICC_w", "gICC_b")

for(j in seq(length(tau_w))){
  
  lambda0 <- 1+baseline_constant[j]*seq(0,(J[j]-1))
  k <- rep( seq(2,J[j]), (n[j]/(J[j]-1)) )
  pi_b <- c(0, rep((1/(J[j]-1)), (J[j]-1)))

  design_varA <- sandwich_var(n=n[j], m=m[j], J=J[j], lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta=betaA)
  
  gICC_results_fixed_tauw[j,] <- cbind(n[j], m[j], J[j], baseline_constant[j], tau_w[j], tau_b[j], design_varA$icc_w, design_varA$icc_b)

}

#### Fixed tau_b simulation scenarios ####
m <- as.numeric(sim_scenarios_fixed_taub[,1])
n <- as.numeric(sim_scenarios_fixed_taub[,2])
J <- as.numeric(sim_scenarios_fixed_taub[,3])
tau_w <- as.numeric(sim_scenarios_fixed_taub[,4])
tau_b <- as.numeric(sim_scenarios_fixed_taub[,5])
beta0 <- 0
betaA <- 0.4
Cp <- 1
baseline_constant <- as.numeric(sim_scenarios_fixed_taub[,6])

gICC_results_fixed_taub <- matrix(NA, nrow=length(m), ncol=8)
colnames(gICC_results) <- c("n","m","J","baseline_constant","tau_w", "tau_b", "gICC_w", "gICC_b")

for(j in seq(length(tau_w))){
  
  lambda0 <- 1+baseline_constant[j]*seq(0,(J[j]-1))
  k <- rep( seq(2,J[j]), (n[j]/(J[j]-1)) )
  pi_b <- c(0, rep((1/(J[j]-1)), (J[j]-1)))
  
  design_varA <- sandwich_var(n=n[j], m=m[j], J=J[j], lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w[j], tau_kb=tau_b[j], beta=betaA)
  
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