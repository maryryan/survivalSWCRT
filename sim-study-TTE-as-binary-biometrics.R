#### LIBRARIES ####
library(tidyverse)
library(pracma)
library(survival)
library(stabledist)
library(latex2exp)
library(gee)
library(geesmv)

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
        gee_ind <- gee(status ~ Z + factor(period), id=id, corstr="independence", family=binomial)
        gee_exch <- gee(status ~ Z + factor(period), id=id, corstr="exch", family=binomial)
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
    results_gee_ind <- rbind(results_gee_ind, summary(gee_ind)$coefficients)
    results_gee_exch <- rbind(results_gee_exch, summary(gee_ex)$coefficients)
    
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
  return(list(survival=results,
              gee_ind=results_gee_ind,
              gee_exch=results_gee_exch))
}


surSIMULATESW_with_gee <- function(n, J, cv.c, cv.p, lambda0, beta, pa, p0, tau_b, tau_w, Cp, nrep, k, p.max, strat=TRUE, FE=FALSE){
  
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
  
  # Function to calculate lambda0 based on the desired administrative censoring rate
  #                                        for the control group pa
  # lambdaDET <- function(Cp, kappa, pa){
  #   lambda0 <- (1/Cp)*(-log(pa))^(1/kappa)
  #   return(lambda0)
  # }
  
  # Function to calculate theta based on Kendall’s tau
  # thetaDET <- function(tau){
  #   theta <- 0.5*(1/tau - 1)
  #   return(theta)
  # }
  
  # lambda0 <- lambdaDET(Cp, kappa, pa)
  #theta <- thetaDET(tau)
  results <- results_surv <- results_gee <- NULL
  
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
    
    #mz <- c(sum(mv[1:(floor(n/2))]), sum(mv[(floor(n/2)+1):n])) # subjects in each arm
    mt <- sum(mv)#sum(mv*max(k)) # total observations
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
    #Z <- rep(c(0, 1), mz) 
    
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
    
    # testStat_naive <- (coef-beta)/summary(survival)$coefficients[3]
    # testStat_robust <- (coef-beta)/summary(survival)$coefficients[4]
    # testStat_MD <- (coef-beta)/BC_results[9]
    # testStat_KC <- (coef-beta)/BC_results[10]
    # testStat_FG <- (coef-beta)/BC_results[11]
    # 
    # testStat_results <- c(testStat_naive, testStat_robust,
    #                       testStat_MD, testStat_KC, testStat_FG)
    
    results_surv <- rbind(results_surv, c(summary(survival)$coefficients, n, J, mean(mv), beta, tau_b, tau_w, 1-mean(ctrl), BC_results$outbeta))#, testStat_results)
    
    # binarizing results #
    gee_indep <- gee(status ~ Z + factor(period), id=id, corstr="independence", family=binomial)
    gee_exch <- gee(status ~ Z + factor(period), id=id, corstr="exchangeable", family=binomial)
    gee_data <- cbind(id,status,Z,period)
    colnames(gee_data) <- c("id","status","Z","period")
    gee_data <- gee_data %>% as.data.frame() %>% mutate(period = as.factor(period)) %>%as.data.frame()
    # fg correction #
    gee_indep_fg <- GEE.var.fg(status ~ Z + period, id="id", data=gee_data, corstr="independence", family=binomial)$cov.beta[2]
    gee_exch_fg <- GEE.var.fg(status ~ Z + period, id="id", data=gee_data, corstr="exchangeable", family=binomial)$cov.beta[2]
    # kc correction #
    gee_indep_kc <- GEE.var.kc(status ~ Z + period, id="id", data=gee_data, corstr="independence", family=binomial)$cov.beta[2]
    gee_exch_kc <- GEE.var.kc(status ~ Z + period, id="id", data=gee_data, corstr="exchangeable", family=binomial)$cov.beta[2]
    # md correction #
    gee_indep_md <- GEE.var.md(status ~ Z + period, id="id", data=gee_data, corstr="independence", family=binomial)$cov.beta[2]
    gee_exch_md <- GEE.var.md(status ~ Z + period, id="id", data=gee_data, corstr="exchangeable", family=binomial)$cov.beta[2]
    
    results_gee <- rbind(results_gee, c(summary(gee_indep)$coefficients[2,],sqrt(gee_indep_fg),sqrt(gee_indep_kc),sqrt(gee_indep_md), 
                                        summary(gee_exch)$coefficients[2,],sqrt(gee_exch_fg),sqrt(gee_exch_kc),sqrt(gee_exch_md)))
    
    if( i %in% c(nrep*seq(0.1,1,0.1)) ) print(paste0(Sys.time(), ": ", (i/nrep)*100, "% done"))
    
  }#end nrep
  results <- cbind(results_surv, results_gee)
  colnames(results) <- c("coef", "exp(coef)", "se(coef)", "robust se",
                         "z", "Pr(>|z|)",
                         "n", "J", "empirical J", "beta", "tau_b", "tau_w", "empirical p0",
                         "coef","exp(coef)","se(coef)","robust se","z","Pr(>|z|)",
                         "MB-stderr","BC0-stderr",#"score",#"RB-stderr",
                         "MD-stderr","KC-stderr","FG-stderr",
                         "score0", "BC0-B0",
                         "nom0","score02","nom02",
                         "MD-B","KC-B","FG-B",#MBN-stderr",
                         #"RBMD-stderr","RBKC-stderr","RBFG-stderr","RBMBN-stderr"
                         "gee indep coef", "gee indep naive se", "gee indep naive z", "gee indep robust se", "gee indep robust z",
                         "gee indep fg se","gee indep kc se","gee indep md se",
                         "gee exch coef", "gee exch naive se", "gee exch naive z", "gee exch robust se", "gee exch robust z",
                         "gee exch fg se","gee exch kc se","gee exch md se"
  )
  
  if (beta == 0){
    name <- "Size"
  } else{
    name <- "Power"}
  cv.c_name <- cv.c*10
  #save(results, file = paste0(name, "_beta_", cv.c_name, ".RData"))
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

#### SIMULATIONS ####
set.seed(8888)
simulation_scenarios <- list(betaA=c(rep(0.35,2),rep(0.4,4),rep(0.45,3),rep(0.5,3),0.55,0.6,rep(0.65,3),rep(0.7,3)),
                              n=c(rep(30,2),rep(21,2),rep(30,2),rep(21,2),30,rep(15,2),21,rep(15,2),rep(9,2),15,rep(9,3)),
                              m=c(40,50,40,50,15,25,25,40,15,40,50,15,25,15,40,50,15,25,40,50))

cv.c <- 0; cv.p <- 0
J <-4
beta0 <- 0

pa <- 0.2; p0 <- 0.3
tau_w <- 0.05; tau_b <- 0.01
Cp <- 1

nrep <- 2000

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

for(s in seq(length(simulation_scenarios[[1]]))){
  
  lambda0 <- 1+0.2*seq(0,(J-1))
  k <- rep( seq(2,J), (n[s]/(J-1)) )
  pi_b <- c(0, rep((1/(J-1)), (J-1)))
  
  n <- simulation_scenarios[["n"]][s]
  m <- simulation_scenarios[["m"]][s]
  betaA <- simulation_scenarios[["betaA"]][s]
  
  # simulate and analyze survival data under null and alternative distributions #
  empirical_var0 <- surSIMULATESW_with_gee(n=n, J=m, cv.c=cv.c, cv.p=cv.p, 
                                  lambda0=lambda0, beta=beta0, pa=pa, p0=p0,
                                  tau_b=tau_b, tau_w=tau_w, Cp=Cp, nrep=nrep, k=k, p.max=J)
  empirical_varA <- surSIMULATESW_with_gee(n=n, J=m, cv.c=cv.c, cv.p=cv.p, 
                                  lambda0=lambda0, beta=betaA, pa=pa, p0=p0,
                                  tau_b=tau_b, tau_w=tau_w, Cp=Cp, nrep=nrep, k=k, p.max=J)
  
  # predict variance under null and alternative distributions#
  design_var0 <- sandwich_var(n=n,m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta=beta0)
  design_varA <- sandwich_var(n=n,m=m,J=J, lambda0=lambda0,tau=Cp, pi_b=pi_b, tau_kw=tau_w, tau_kb=tau_b, beta=betaA)
  
  # calculate predicted power based on your sandwich variance and z test #
  design_power_z <- pnorm(abs(betaA)/design_varA$se_beta - qnorm(0.975))
  
  # calculate predicted power based on your sandwich variance and t test #
  design_power_t <- pt((abs(betaA)/design_varA$se_beta - qt(0.975, n-1)), n-1) 
  
  ## Null results ##
  ASE0 <- mean(empirical_var0[,4])
  ESE0 <- sd(empirical_var0[,1])
  naive0 <- mean(empirical_var0[,3])
  BC0_avg0 <- mean(empirical_var0[,21])
  BC_MD_avg0 <- mean(empirical_var0[,22])
  BC_KC_avg0 <- mean(empirical_var0[,23])
  BC_FG_avg0 <- mean(empirical_var0[,24])
  
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
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-1), 1, 0 )
    }
    return(mean(reject))
  })
  
  t_test_gee_typeI <- apply(testStat_gee_results0,2, function(x){
    
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-J), 1, 0 )#dof is n-(1+(J-1)) since using period fixed effects
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
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-1), 1, 0 )
    }
    mean(reject)
    
  })
  
  t_test_gee_power <- apply(testStat_gee_resultsA,2, function(x){
    reject <- rep(NA,length(x))
    for(i in seq(length(x))){
      reject[i] <- ifelse( abs(x[i]) >= qt(0.975, n-J), 1, 0 )
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
