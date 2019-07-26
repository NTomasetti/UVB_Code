rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
source('Functions/scoreVB.R')
source('Functions/metropolisHastings.R')
sourceCpp('Functions/ARP.cpp')

# This is set to each of 1:500
rep <- 1 

# Store results, set counter
forecast <- list()
counter <- 1

# Set up updating
T <- 500
initialBatch <- 100
updateSize <- 25
batches <- (T - initialBatch) / updateSize + 1
data <- c(100, seq(100 + initialBatch, T+100, length.out = batches))

# Other setup
lags <- 3
MCsamples <- 2000 # To estimate forecast densities as a Monte Carlo Estimate
ISsamples <- 100 # Number of theta draws for UVB-IS

dim <- 2 + lags
priorMean <- rep(0, dim)
priorVar <- diag(10, dim)
priorVarInv <- solve(priorVar)
priorLinv <- solve(t(chol(priorVar)))

set.seed(rep)

# Simulate Data

mu <- rnorm(1)

stationary <- FALSE
sd <- 1
while(!stationary){
  phi <- rnorm(lags, 0, sd)
  roots <- polyroot(c(1, -phi))
  if(all(Mod(roots) > 1)){
    stationary <- TRUE
  } else {
    sd <- sd * 0.95
  }
}

sigmaSq <- 1/rgamma(1, 5, 5)
x <- rnorm(lags, 0, sqrt(sigmaSq))
for(t in (lags+1):(T+100+1)){
  x <- c(x, mu + sum(phi * (x[t-(1:lags)] -mu)) + rnorm(1, 0, sqrt(sigmaSq)))
}

# Forecast support
support <- seq(min(x)-sd(x), max(x)+sd(x), length.out = 1000)

# MCMC for exact forecast densities
for(t in 1:batches){
  
  # Fit MCMC
  MCMCfit <- metropolisHastings(start = rep(0, dim),
                                iter = 15000,
                                model = ARjointDensSingle,
                                burnIn = 10000,
                                x = x[(data[1]+1 - lags):data[t+1]],
                                lags = lags,
                                mu = priorMean,
                                varInv = priorVarInv)
  
  # Set up forecast density
  densMCMC <- rep(0, 1000)
  # Forecast
  for(i in 1:MCsamples){
    sigSq <- exp(MCMCfit[i, 1])
    mu <- MCMCfit[i, 2]
    phi <- MCMCfit[i, 3:(2 + lags)]
    mean <- mu + phi %*% (x[data[t+1]+1-(1:lags)] - mu)
    densMCMC <- densMCMC + dnorm(support, mean, sqrt(sigSq)) / MCsamples
  }
  lower <- max(which(support < x[data[t+1]+1]))
  upper <- min(which(support > x[data[t+1]+1]))
  
  if(lower == -Inf){
    lsMCMC <- log(densMCMC[upper])
  } else if(upper == Inf) {
    lSMCMC <- log(densMCMC[lower])
  } else {
    lsMCMC <- log(linearInterpolate(support[lower], support[upper], densMCMC[lower], densMCMC[upper], x[data[t+1]+1]))
  }
  
  forecast[[counter]] <- data.frame(t = data[t + 1],
                                    ls = lsMCMC,
                                    inference = 'MCMC',
                                    K = 1:3,
                                    runTime = NA,
                                    ESS = NA,
                                    ELBO = NA,
                                    id = rep)
  counter <- counter + 1
}

## Apply VB approximations with K = 1, 2, or 3 component mixture approximations
for(K in 1:3){
  # Total runtime
 VBTotal <- UVBTotal <- ISTotal <- 0

  # Initial Lambda: posterior means, then log of posterior standard dev, then un-normalised mixture weights (z, where mixture weights = exp(z) / sum(exp(z)))
  # Add a small random peturbation
  lambda <- c(rep(c(-1, rep(0, lags+1)), K) + rnorm(dim*K, 0, 0.1),
              rnorm(dim*K, -1, 0.1),
              rep(1, K))
  
  # Store lambda values
  VBfit <- UVBfit <- ISfit <- matrix(0, length(lambda), batches)
  for(t in 1:batches){
   
    if(t == 1){
      
      # Fit SVB
      
      startVB <- Sys.time()
      VB <- scoreVB(data = x[(data[1]+1 - lags):data[t+1]],
                    lambda = lambda,
                    model = AR3VB,
                    dimTheta = 5,
                    mix = K,
                    S = 25,
                    mean = matrix(priorMean, ncol = 1),
                    SigInv = array(priorVarInv, dim = c(5, 5, 1)),
                    weights = 1,
                    dets = sqrt(det(priorVarInv)))
      endVB <- Sys.time() - startVB
      VBfit[,t] <- UVBfit[,t] <- ISfit[,t] <- VB$lambda
      ELBO <- rep(VB$LB[VB$iter - 1], 3)

      
    } else {
      
      # Fit SVB
      startVB <- Sys.time()
      VB <- scoreVB(data = x[(data[1]+1 - lags):data[t+1]],
                    lambda = lambda,
                    model = AR3VB,
                    dimTheta = 5,
                    mix = K,
                    S = 25,
                    mean = matrix(priorMean, ncol = 1),
                    SigInv = array(priorVarInv, dim = c(5, 5, 1)),
                    weights = 1,
                    dets = sqrt(det(priorVarInv)))
      endVB <- Sys.time() - startVB
      VBfit[,t] <- VB$lambda
      ELBO <- VB$LB[VB$iter - 1]
 
      # UVB
      
      # Set up priors from previous fit, also grab variance determinant component of the MVN prior
      startUVB <- Sys.time()
      updateMean <- matrix(UVBfit[1:(dim*K), t-1], ncol = K)
      updateVarInv <- array(0, dim = c(dim, dim, K))
      dets <- rep(0, K)
      for(k in 1:K){
        sd <- exp(UVBfit[dim*K + dim*(k-1) + 1:dim, t-1])
        updateVarInv[,,k] <- diag(1/sd^2)
        dets[k] <- 1 / prod(sd)
      }
      updateZ <- UVBfit[dim*K*2 + 1:K, t-1] 
      # Normalise weights
      updateWeight <- exp(updateZ) / sum(exp(updateZ))
      # Fit model
      UVB <- scoreVB(data = x[(data[t]+1-lags):data[t+1]],
                     lambda = UVBfit[,t-1],
                     model = AR3VB,
                     dimTheta = 5,
                     mix = K,
                     S = 25,
                     mean = updateMean,
                     SigInv = updateVarInv,
                     dets = dets,
                     weights = updateWeight)
      
      endUVB <- Sys.time() - startUVB
      UVBfit[,t] <- UVB$lambda
      ELBO <- c(ELBO, UVB$LB[UVB$iter - 1])
        
      #UVB-IS
      startIS <- Sys.time()
      #Extract means, variance, mixture weights of prior distribution
      
      ISMean <- matrix(ISfit[1:(dim*K), t-1], ncol = K)
      ISSd <- matrix(exp(ISfit[dim*K + 1:(dim*K), t-1]), ncol = K)
      ISZ <- ISfit[2*dim*K + 1:K, t-1]
      ISWeight <- exp(ISZ) / sum(exp(ISZ))
      # Draw samples from the prior distribution
      ISdraw <- matrix(0, ISsamples, dim)
      for(i in 1:ISsamples){
        group <- sample(1:K, 1, prob = ISWeight)
        ISdraw[i, ] <- mvtnorm::rmvnorm(1, ISMean[,group], diag(ISSd[,group]^2))
      }
      
      # Calcualte the density of these draws under the prior distributoin
      qIS <- rep(0, ISsamples)
      for(k in 1:K){
        qIS <- qIS + ISWeight[k] * mvtnorm::dmvnorm(ISdraw, ISMean[,k], diag(ISSd[,k]^2))
      }
      # Calculate the log-joint density of these samples using C++ functions
      if(K == 1){
	# The prior is a multivariate normal with parameters ISMean, Inverse Variance = diag(1/SD^2)
        pIS <- ARjointDensSingleMatrix(x[(data[t]+1 - lags):(data[t+1])], ISdraw, ISMean, diag(c(1 / ISSd^2)), lags)
      } else {
	# Alternatively the prior is a mixture of multivariate normals, and we calculate the inverse variance for each component.
        varInv <- array(0, dim = c(2 + lags, 2 + lags, K))
        dets <- rep(0, K)
        for(k in 1:K){
          varInv[,,k] <- diag(c(1 / ISSd[,k]^2))
          dets[k] <- prod(1 / ISSd[,k])
        }
        pIS <- ARjointDensMixMatrix(x[(data[t]+1 - lags):(data[t+1])], ISdraw, ISMean, varInv, dets, ISWeight, lags)
      }
      # Run VB Update
      VBIS <- UVBIS(lambda = ISfit[,t-1],
                    qScore = AR3Score,
                    samples = ISdraw,
                    dSamples = qIS,
                    logjoint = c(log(pIS)),
                    maxIter = 2000,
                    mix = K)

      endIS <-  Sys.time() - startIS
      ISfit[,t] <- VBIS$lambda
      ELBO <- c(ELBO, VBIS$LB[VBIS$iter - 1])

    }
    
        
    if(attr(endVB, 'units') == 'mins'){
      endVB <- endVB * 60
    }
    if(t == 1){
      endIS <- endVB
      endUVB <- endVB
    } else {
      if(attr(endUVB, 'units') == 'mins'){
        endUVB <- endUVB * 60
      }
      if(attr(endIS, 'units') == 'mins'){
        endIS <- endIS * 60
      }
    }
   
    VBTotal <- VBTotal + as.numeric(endVB)
    UVBTotal <- UVBTotal + as.numeric(endUVB)
    ISTotal <- ISTotal + as.numeric(endIS)
    
    # Draw samples for forecasting
    VBMean <- matrix(VBfit[1:(dim*K), t], ncol = K)
    VBSd <- matrix(exp(VBfit[dim*K + 1:(dim*K), t]), ncol = K)
    VBZ <- VBfit[2*dim*K + 1:K, t]
    VBWeight <- exp(VBZ) / sum(exp(VBZ))
    VBdraw <- matrix(0, MCsamples, dim)
    for(i in 1:MCsamples){
      group <- sample(1:K, 1, prob = VBWeight)
      VBdraw[i, ] <- mvtnorm::rmvnorm(1, VBMean[,group], diag(VBSd[,group]^2))
    }
    
    if(t > 1){
      # Also draw UVB and UVBIS for T_2 onwards
      UVBMean <- matrix(UVBfit[1:(dim*K), t], ncol = K)
      UVBSd <- matrix(exp(UVBfit[dim*K + 1:(dim*K), t]), ncol = K)
      UVBZ <- UVBfit[2*dim*K + 1:K, t]
      UVBWeight <- exp(UVBZ) / sum(exp(UVBZ))
      UVBdraw <- matrix(0, MCsamples, dim)
      for(i in 1:MCsamples){
        group <- sample(1:K, 1, prob = UVBWeight)
        UVBdraw[i, ] <- mvtnorm::rmvnorm(1, UVBMean[,group], diag(UVBSd[,group]^2))
      }
    
      ISMean <- matrix(ISfit[1:(dim*K), t], ncol = K)
      ISSd <- matrix(exp(ISfit[dim*K + 1:(dim*K), t]), ncol = K)
      ISZ <- ISfit[2*dim*K + 1:K, t]
      ISWeight <- exp(ISZ) / sum(exp(ISZ))
      ISdraw <- matrix(0, MCsamples, dim)
      for(i in 1:MCsamples){
        group <- sample(1:K, 1, prob = ISWeight)
        ISdraw[i, ] <- mvtnorm::rmvnorm(1, ISMean[,group], diag(ISSd[,group]^2))
      }
    }
   
    # Forecast
    # Set up forecast densities
    densVB <- densUVB <- densIS <- rep(0, 1000)
 
    
    for(i in 1:MCsamples){
      # SVB forecast density
      sigSqVB <- exp(VBdraw[i, 1])
      muVB <- VBdraw[i, 2]
      phiVB <- VBdraw[i, 3:(2 + lags)]
      meanVB <- muVB + phiVB %*% (x[(data[t+1]+1-(1:lags))] - muVB)
      densVB <- densVB + dnorm(support, meanVB, sqrt(sigSqVB)) / MCsamples
      
      if(t > 1){
        # UVB forecast density
        sigSqUVB <- exp(UVBdraw[i, 1])
        muUVB <- UVBdraw[i, 2]
        phiUVB <- UVBdraw[i, 3:(2 + lags)]
        meanUVB <- muUVB + phiUVB %*% (x[(data[t+1]+1-(1:lags))] - muUVB)
        densUVB <- densUVB + dnorm(support, meanUVB, sqrt(sigSqUVB)) / MCsamples
        
        #UVBIS forecast density
        sigSqIS <- exp(ISdraw[i, 1])
        muIS <- ISdraw[i, 2]
        phiIS <- ISdraw[i, 3:(2 + lags)]
        meanIS <- muIS + phiIS %*% (x[(data[t+1]+1-(1:lags))] - muIS)
        densIS <- densIS + dnorm(support, meanIS, sqrt(sigSqIS)) / MCsamples
        
       
      } 
    }
    
    if(t == 1){
      densUVB <- densIS <- densVB
    }
    lower <- max(which(support < x[data[t+1]+1]))
    upper <- min(which(support > x[data[t+1]+1]))
    
    if(lower == -Inf){
      lsVB <- log(densVB[upper])
      lsUVB <- log(dens(UVB[upper]))
      lsIS <- log(densIS[upper])
    } else if(upper == Inf) {
      lsVB <- log(densVB[lower])
      lsUVB <- log(densUVB[lower])
      lsIS <- log(densIS[lower])
    } else {
      # Linear interpolation of forecast density at observed value
      lsVB <- log(linearInterpolate(support[lower], support[upper], densVB[lower], densVB[upper], x[data[t+1]+1]))
      lsUVB <- log(linearInterpolate(support[lower], support[upper], densUVB[lower], densUVB[upper], x[data[t+1]+1]))
      lsIS <- log(linearInterpolate(support[lower], support[upper], densIS[lower], densIS[upper], x[data[t+1]+1]))
    }
    
    # Save results
    forecast[[counter]] <- data.frame(t = data[t+1],
                                      ls = c(lsVB, lsUVB, lsIS),
                                      inference = c('VB', 'UVB', 'UVB-IS'),
                                      K = K,
                                      runTime = c(VBTotal, UVBTotal, ISTotal),
                                      ELBO = ELBO,
                                      id = rep)
    counter <- counter + 1
  }
}

forecast <- do.call(rbind.data.frame, forecast)

write.csv(forecast, paste0('AR3/iter_', rep, '.csv'), row.names = FALSE)

