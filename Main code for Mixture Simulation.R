rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
source('Functions/mixture normal funs.R')
sourceCpp('Functions/mixtureNormal.cpp')

# Set from 1 to 500
rep <- 1

set.seed(rep)

# Initialise everything
N <- 100
T <- 100
initialBatch <- 10
updateSize <- 10
batches <- (T - initialBatch) / updateSize + 1
data <- seq(initialBatch, T, length.out = batches)
MCsamples <- 2000
ISsamples <- 100

# Simulate data
mu <- rnorm(2, 0, 0.25)
sigmaSq <- runif(2, 1, 2)
# True group labels
group <- sample(0:1, N, replace = TRUE)
y <- matrix(0, T, N)
results <- list()
counter <- 1

for(i in 1:N){
  y[,i] <- rnorm(T, mu[group[i]+1], sqrt(sigmaSq[group[i]+1]))
}

# Prior components
priorMean <- rep(0, 4)
priorVar <- diag(10, 4)
priorLinv <- solve(t(chol(priorVar)))
priorVarInv <- solve(priorVar)

# MCMC
for(t in 1:batches){
  # Fit MCMC
  MCMCfit <- MixNormMCMC(y = y[1:data[t], ],
                         reps = 15000,
                         drawT = c(0.5, 0.5, 0, 0),
                         drawK = sample(0:1, N, replace = TRUE),
                         hyper = list(mean = priorMean, varInv = priorVarInv),
                         thin = 1,
                         stepsize = 0.01)
  keep <- floor(seq(10001, 15000, length.out = MCsamples))
  
  # Predict class labels via draws of K taken in MCMC
  class <- rep(0, N)
  # Forecast
  for(j in 1:N){
    class1 <- 0
    for(i in 1:MCsamples){
      class1 <- class1 + MCMCfit$K[keep[i], j] / MCsamples
    }
    class[j] <- class1 > 0.5
  }
  score <- max(sum(class == group), sum(class != group)) / N
  
  results[[counter]] <- data.frame(t = data[t],
                              score = score,
                              inference = 'MCMC',
                              K = 1:3,
                              runTime = NA,
                              ELBO = NA,
                              id = rep)
  counter <- counter + 1
}

# Fit VB for K = 1, 2, 3 mixture approximations
for(K in 1:3){
  # Start Time
  VBTotal <- 0
  UVBTotal <- 0 
  ISTotal <- 0
    
  # Initial Lambda, 4 theta values & K mixture components
  # First elements: Posterior Means, t1 mix1, t2 mix1, ..., t2 mix2 etc
  # Next: Posterior log sd
  # Last: unnormalised mixture weights
  lambda <- c(rnorm(4*K, 0, 0.1), rnorm(4*K, -1, 0.5), rnorm(K, 0, 0.25))
    
  # Store converged lambda
  VBfit <- UVBfit <- ISfit <- matrix(0, length(lambda), batches)
    
  for(t in 1:batches){

    
    if(t == 1){
      # Fit SVB
      startVB <- Sys.time()
      VB <- fitVBScoreMN(data = y[1:data[t],],
                         S = 25,
                         lambda = lambda,
                         mixPrior = 1,
                         mixApprox = K,
                         mean = matrix(priorMean, ncol = 1),
                         SigInv = array(priorVarInv, dim = c(4, 4, 1)),
                         probK = rep(0.5, N),
                         weights = 1,
                         dets = det(priorVarInv) ^ 0.5)
      endVB <- Sys.time() - startVB
      VBfit[,t] <- UVBfit[,t] <- ISfit[,t]  <- VB$lambda
      ELBO <- rep(VB$LB[VB$iter], 3)
    } else {
      # SVB
      startVB <- Sys.time()
      VB <- fitVBScoreMN(data = y[1:data[t],],
                         S = 25,
                         lambda = lambda,
                         mixPrior = 1,
                         mixApprox = K,
                         mean = matrix(priorMean, ncol = 1),
                         SigInv = array(priorVarInv, dim = c(4, 4, 1)),
                         probK = rep(0.5, N),
                         weights = 1,
                         dets = det(priorVarInv) ^ 0.5)
      endVB <- Sys.time() - startVB
      VBfit[,t] <- VB$lambda
      ELBO <- VB$LB[VB$iter - 1]
      
      # UVB
      startUVB <- Sys.time()
      # Set up prior
      updateMean <- matrix(UVBfit[1:(4*K), t-1], ncol = K)
      updateSd <- matrix(exp(UVBfit[4*K + 1:(4*K), t-1]), ncol = K)
      updateVarInv <- array(0, dim = c(4, 4, K))
      dets <- rep(0, K)
      for(k in 1:K){
        updateVarInv[,,k] <- diag(1/updateSd[,k]^2)
        dets[k] <- 1 / prod(updateSd[,k])
      }
      updateZ <- UVBfit[8*K + 1:K, t-1]
      updateWeight <- exp(updateZ) / sum(exp(updateZ))
      # Update lambda
      UVB <- fitVBScoreMN(data = y[(data[t-1]+1):data[t],],
                         S = 25,
                         lambda = UVBfit[,t-1],
                         mixPrior = K,
                         mixApprox = K,
                         mean = updateMean,
                         SigInv = updateVarInv,
                         probK = probKUVB,
                         weights = updateWeight,
                         dets = dets)
      endUVB <- Sys.time() - startUVB
      UVBfit[,t] <- UVB$lambda
      ELBO <- c(ELBO, UVB$LB[UVB$iter])
      
      # UVB-IS
      
      startIS <- Sys.time()
      
      # Set up prior
      ISMean <- matrix(ISfit[1:(4*K), t-1], ncol = K)
      ISSd <- matrix(exp(ISfit[4*K + 1:(4*K), t-1]), ncol = K)
      ISZ <- ISfit[2*4*K + 1:K, t-1]
      ISWeight <- exp(ISZ) / sum(exp(ISZ))
      # Sample from prior (old approx)
      ISdraw <- matrix(0, ISsamples, 4)
      for(i in 1:ISsamples[isrep]){
        groupIS <- sample(1:K, 1, prob = ISWeight)
        ISdraw[i, ] <- rnorm(4, ISMean[,groupIS], ISSd[,groupIS])
      }
      # Calcualte density under the old distribution
      qIS <- rep(0, ISsamples)
      for(k in 1:K){
        qIS <- qIS + ISWeight[k] * mvtnorm::dmvnorm(ISdraw, ISMean[,k], diag(ISSd[,k]^2))
      }
      # Calculate the log-joint distribution through the UVB prior recursion
      # log liklihood + log(p(k)) was calculated already for each group, subtract log(p(k)) and reconstruct with weights
      loglik <- rep(0, ISsamples)
      for(j in 1:N){
        ll <- matrix(0, ISsamples, 2)
        for(i in 1:ISsamples){
          ll[i, ] <- postK(y[(data[t-1]+1):data[t],j], ISdraw[i,], c(0, 0))
        }
        loglik <- loglik + log((1 - probKIS[j]) * exp(ll[,1]) +
                                 probKIS[j] * exp(ll[,2]))
      }
      loglik <- loglik - max(loglik)
      # Prior is qIS
      logjoint <- loglik + log(qIS)  
      # Run VB Update
      IS <- ISUVB(lambda = ISfit[,t-1, isrep],
                  qScore = mixNormScore,
                  samples = ISdraw,
                  dSamples = qIS,
                  logjoint = logjoint,
                  maxIter = 2000,
                  mix = K)
      endIS <- Sys.time() - startIS
    
      ISfit[,t] <- IS$lambda
      ELBO <- c(ELBO, IS$LB[IS$iter])
    
    }
      
    if(attr(endVB, 'units') == 'mins'){
      endVB <- endVB * 60
    }
    if(t == 1){
      endUVB <- endIS <- endVB
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
    
    # Draw samples from VB
    VBMean <- matrix(VBfit[1:(4*K), t], ncol = K)
    VBSd <- matrix(exp(VBfit[4*K + 1:(4*K), t]), ncol = K)
    VBZ <- VBfit[2*4*K + 1:K, t]
    VBWeight <- exp(VBZ) / sum(exp(VBZ))
    VBdraw <- matrix(0, MCsamples, 4)
    for(i in 1:MCsamples){
      groupVB <- sample(1:K, 1, prob = VBWeight)
      VBdraw[i, ] <- mvtnorm::rmvnorm(1, VBMean[,groupVB], diag(VBSd[,groupVB]^2))
    }
    
    if(t > 1){
      # Also draw from UVB
      UVBMean <- matrix(UVBfit[1:(4*K), t], ncol = K)
      UVBSd <- matrix(exp(UVBfit[4*K + 1:(4*K), t]), ncol = K)
      UVBZ <- UVBfit[2*4*K + 1:K, t]
      UVBWeight <- exp(UVBZ) / sum(exp(UVBZ))
      UVBdraw <- matrix(0, MCsamples, 4)
      for(i in 1:MCsamples){
        groupUVB <- sample(1:K, 1, prob = UVBWeight)
        UVBdraw[i, ] <- mvtnorm::rmvnorm(1, UVBMean[,groupUVB], diag(UVBSd[,groupUVB]^2))
      }
      # And UVB-IS
      ISMean <- matrix(ISfit[1:(4*K), t], ncol = K)
      ISSd <- matrix(exp(ISfit[4*K + 1:(4*K), t]), ncol = K)
      ISZ <- ISfit[2*4*K + 1:K, t]
      ISWeight <- exp(ISZ) / sum(exp(ISZ))
      ISdraw <- matrix(MCsamples, 4)
      for(i in 1:MCsamples){
        groupIS <- sample(1:K, 1, prob = ISWeight)
        ISdraw[i, ] <- mvtnorm::rmvnorm(1, ISMean[,groupIS], diag(ISSd[,groupIS]^2))
      }
    }
     
    # Calculate p(k = 1) using each set of samples
    # These will store the log posterior for k for each unit i, value of k, and sampled value of theta
    probKVB <- probKUVB <- probKIS <- array(0, dim = c(MCsamples, 2, N))

    for(i in 1:MCsamples){
      for(j in 1:N){
        probKVB[i, , j] <- postK(y[1:data[t],j], VBdraw[i,], c(log(0.5), log(0.5)))
        if(t > 1){
          probKUVB[i, , j] <- postK(y[1:data[t],j], UVBdraw[i,], c(log(0.5), log(0.5)))
          probKIS[i, , j] <- postK(y[1:data[t],j], ISdraw[i,], c(log(0.5), log(0.5)))
        }
      }
    }
    if(t == 1){
      probKUVB <- probKVB
      probKIS <- probKVB
    }
    
    # Average the log posteriors across different theta samples
    probKVB <- apply(probKVB, 2:3, mean)
    probKUVB <- apply(probKUVB, 2:3, mean)
    probKIS <- apply(probKIS, 2:3, mean)
    
    # Normalise log posterior
    probKVB <- apply(probKVB, 2, function(x) {y = x - max(x); exp(y[2]) / sum(exp(y))})
    probKUVB <- apply(probKUVB, 2, function(x) {y = x - max(x); exp(y[2]) / sum(exp(y))})
    probKIS <- apply(probKIS, 2, function(x) {y = x - max(x); exp(y[2]) / sum(exp(y))})
 
    # p-hat k == 1 if posterior prob > 0.5
    classVB <- probKVB > 0.5
    classUVB <- probKUVB > 0.5
    classIS <- probKIS > 0.5
    scoreVB <- max(sum(classVB == group), sum(classVB != group)) / N
    scoreUVB <- max(sum(classUVB == group), sum(classUVB != group)) / N
    scoreIS <-  max(sum(x == group), sum(x != group)) / N
      
    results[[counter]] <-  data.frame(t = data[t],
                                      score = c(scoreVB, scoreUVB, scoreIS),
                                      inference = c('VB', 'UVB', 'UVB-IS'),
                                      runTime = c(VBTotal, UVBTotal, ISTotal), 
                                      ELBO = ELBO,
                                      id = rep)
    counter <- counter + 1
  }
}
results <- do.call(rbind.data.frame, results)
write.csv(results, paste0('mixNorm/iter_', rep, '.csv'), row.names = FALSE)
