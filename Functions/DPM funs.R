#' Score based VB for the Dirichlet Process Mixture
#' 
#' This function runs a VB optimisation for the vehicle dirichlet process mixture model defined in the UVB paper.
#' It may be applied to either SVB or UVB with a suitable choice of the prior.
#' This relies on functions defined in dpm.cpp to calculate the gradients. It also requires the mvtnorm package.
#' @param data A matrix of observations, where columns correspond to different units (i) and rows to time period (t)
#' @param lambda An 5 x N matrix of parameter values. Each column corresponds to the parameters of a single unique theta cluster. 
#' The first row is the mean of log sigma, the second row is the mean of mu,
#' while the next three rows correspond to the three entires of the lower triangle (cholesky) of the variance matrix
#' @param logPrior A N x N matrix of log prior probabilities for the indicator variables, used for updates. This is the log of the entires of the q(k | y) distribution, one column per vehicle.
#' @param alphaDPM The mass parameter of the Dirichlet Process
#' @param S The number of monte carlo samples per iteration
#' @param priorMean A 2 x N matrix of the prior mean parameter for theta star: log variance (first row) and mu (second row). 
#' Set this to the base measure for all columns for the initial fit, or to the previous VB approximation for updates. This defaults to the zero vector.
#' @param priorSig A 2 x 2 x N array of the variance matrix for each theta star cluster, either the base measure repeated for each slice or the previous VB fit.
#' This defaults to the identity matrix.
#' @param maxIter The maximum number of gradient iterations
#' @param alpha Adam control parameter
#' @param beta1 Adam control parameter
#' @param beta2 Adam control parameter
#' @param threshold ELBO convergence threshold
#' @param suppressProgress if TRUE, suppress printed messages
fitVBDPM <- function(data, lambda, logPrior = NULL, alphaDPM = 1, S = 50,
                     priorMean = matrix(0, 2, N), priorSig = array(rep(c(diag(2)), ncol(data)),c(2, 2, ncol(data))),
                     maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01, suppressProgress = FALSE){
  
  N <- ncol(data)
  T <- nrow(data)
  
  # Control parameters for the optimisation. 

  iter <- 1
  M <- V <- matrix(0, nrow(lambda), ncol(lambda))
  e <- 1e-8
  
  # Store ELBO values, and mean / previous mean of a set of 5 consecutive ELBO values
  LB <- numeric(maxIter)
  meanLB <- 0
  oldMeanLB <- 0
  # Difference between mean LB values
  diff <- threshold + 1
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    
    # Store gradient and control variate h function
    grad <- h <- array(0, dim = c(nrow(lambda), ncol(lambda), S))

    eval <- numeric(S)
    # S many monte carlo samples
    for(s in 1:S){
      
      # Theta Row 1 : Log Standard Deviation
      # Theta Row 2 : Means
      theta <- matrix(0, 2, N)
      prior <- 0
      for(i in 1:N){
        mean <- lambda[1:2, i]
        l <- matrix(c(lambda[3:4, i], 0, lambda[5, i]), 2)
        sigma <- l %*% t(l)
        theta[,i] <- mvtnorm::rmvnorm(1, mean, sigma)
        prior <- prior + mvtnorm::dmvnorm(theta[,i], priorMean[,i], priorSig[, , i], log = TRUE)
      }
      
      # Initial k_1 = 1, 
      k <- 1
      # nk = number of times each value of k is repeated
      nk <- c(1, rep(0, N-1))
      # Active clusters = unique values of k
      active <- 1
      
      # likelihood of y | k, theta
      py <- sum(dnorm(data[,1], theta[2, 1], exp(theta[1, 1]), log = TRUE))
      
      for(i in 2:N){
        # Calculate likelihood of current data under each active cluster
        probY <- vapply(c(active, i), function(x) sum(dnorm(data[,i], theta[2, x], exp(theta[1, x]), log = TRUE)), runif(1))
        
        # add previous log prior values if provided, else use CRP probabilities
        if(!is.null(logPrior)){
          probY <- probY + logPrior[c(active, i), i]
          postProb <- probY - max(probY)
          postProb <- exp(postProb) / sum(exp(postProb))
        } else {
          probY <- probY + log(c(nk[active], alphaDPM) / (alphaDPM + i - 1))
          postProb <- probY - max(probY)
          postProb <- exp(postProb) / sum(exp(postProb))
        }
        # sample K from the posterior probabilty, update log likelihood conditioned on sampled value of k
        newK <- sample(c(active, i), 1, prob = postProb)
        py <- py + probY[which(c(active, i) == newK)] 
        k <- c(k, newK)
        nk[newK] <- nk[newK] + 1
        active <- unique(k)
      }
      # score gradient derivative function
      derivs <- mixtureNormalDPM(theta = theta,
                                 lambdaIn = matrix(lambda, ncol = 1),
                                 logp = py + prior)
      eval[s] <- derivs$val
      
      grad[,,s] <- matrix(derivs$grad[1:(5 * N * mixApprox)], nrow = N, ncol = 5)
      h[,,s] <- grad[,,s] / derivs$val
      
    }
    # transform to vector and estimate control variate via vapply
    gL <- matrix(grad, ncol = S)
    hL <- matrix(h, ncol = S)
    a <- vapply(1:nrow(gL), function(x) cov(gL[x, ], hL[x,], use = 'complete.obs') / var(hL[x,], na.rm = TRUE), 1)
    a[is.na(a)] <- 0
    # estimate gradient as a mean, rearrange resulting vector back to matrix format
    gradient <- vapply(1:nrow(gL), function(x) mean(gL[x, ] - a[x] * hL[x, ], na.rm = TRUE), runif(1))
    gradient <- matix(gradient, N, 5)
    
    gradientSq <- vapply(1:nrow(gL), function(x) mean((gL[x, ] - a[x] * hL[x, ])^2, na.rm = TRUE), runif(1))
    gradientSq <- matrix(gradientSq, N, 5)
    
    LB[iter] <- mean(eval[is.finite(eval)], na.rm = TRUE)
    
    # Apply Adam update
    
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    update <- alpha * Mst / sqrt(Vst + e)
    if(any(is.na(update))){
      print('Break - NA gradient value')
      break
    } else {
      lambda <- lambda + update
    }
    
    # calculate mean elbo over previous five iterations to reduce noise in esimtation and assess convergence
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 50 == 0 & !suppressProgress){
      print(paste0('iter: ', iter, ' ELBO: ', LB[iter]))
    }
    iter <- iter + 1
  }
  iter <- iter - 1
  if(!suppressProgress){
    print(paste0('iter: ', iter, ' ELBO: ', LB[iter]))
  }
 
  return(list(lambda=lambda, LB = LB[1:iter], iter = iter))
}

#' Sample k and theta from the DPM
#'
#' Given the optimal lambda values, sample from the approximation
#' @param data A matrix of observations, where columns correspond to different units (i) and rows to time period (t)
#' @param lambda An 5 x N matrix of parameter values. Each column corresponds to the parameters of a single theta cluster. 
#' The first row is the mean of log sigma, the second row is the mean of mu,
#' while the next three rows correspond to the three entires of the lower (cholesky) triangle of the variance matrix
#' @param logPrior A N x N matrix of log prior probabilities for the indicator variables, used for updates. This is the log of the entires of the q(k | y) distribution
#' @param alphaDPM The mass parameter of the Dirichlet Process
#' @param S The number of samples to generate
VBDPMsample <- function(data, lambda, logPrior = NULL, alphaDPM = 1,  S = 50){
  
  # The sampled values of k
  select <- matrix(0, N, S)
  # The sampled values of theta
  thetaFull <- array(0, dim = c(2, N, S))
  
  for(s in 1:S){
    
    component <- 1
    
    theta <- matrix(0, 2, N)
    for(i in 1:N){
      mean <- lambda[1:2, i]
      l <- matrix(c(lambda[3:4, i], 0, lambda[5, i]), 2)
      sigma <- l %*% t(l)
      theta[,i] <- mvtnorm::rmvnorm(1, mean, sigma)
    }
    thetaFull[,,s] <- theta
    
    k <- 1
    nk <- c(1, rep(0, N-1))
    active <- 1
    select[1, s] <- 1
    for(i in 2:N){
    
      probY <- vapply(c(active, i), function(x) sum(dnorm(data[,i], theta[2, x], exp(theta[1, x]), log = TRUE)), runif(1))
      
      if(!is.null(logPrior)){
        probY <- probY + logPrior[c(active, i), i]
        postProb <- probY - max(probY)
        postProb <- exp(postProb) / sum(exp(postProb))
      } else {
        probY <- probY + log(c(nk[active], alphaDPM) / (alphaDPM + i - 1))
        postProb <- probY - max(probY)
        postProb <- exp(postProb) / sum(exp(postProb))
      }
      
      newK <- sample(c(active, i), 1, prob = postProb)
      k <- c(k, newK)
      nk[newK] <- nk[newK] + 1
      active <- unique(k)
      select[i, s] <- newK
    }
  }
  list(k = select, theta = thetaFull)
}

#' Marginalise the probabilties for k over different values of theta
#' 
#' This function calcules the q(k | y) value as a monte carlo average of S many posterior draws
#' @param data A matrix of observations, where columns correspond to different units (i) and rows to time period (t)
#' @param lambda An 5 x N matrix of parameter values. Each column corresponds to the parameters of a single theta cluster. 
#' The first row is the mean of log sigma, the second row is the mean of mu,
#' while the next three rows correspond to the three entires of the lower (cholesky) triangle of the variance matrix
#' @param logPrior A N x N matrix of log prior probabilities for the indicator variables, used for updates. This is the log of the entires of the q(k | y) distribution
#' @param alphaDPM The mass parameter of the Dirichlet Process
#' @param S The number of samples to average over
VBDPMprobs <- function(data, lambda, logPrior = NULL, alphaDPM = 1, S = 50){
  
  logProb <- matrix(0, N, N)
  mixApprox <- length(weightsZ)
  
  for(s in 1:S){
  
    theta <- matrix(0, 2, N)
    for(i in 1:N){
      mean <- lambda[1:2, i]
      l <- matrix(c(lambda[3:4. i], 0, lambda[5, i]), 2)
      sigma <- l %*% t(l)
      theta[,i] <- mvtnorm::rmvnorm(1, mean, sigma)
    }
    
    
    for(i in 1:N){
      
      probY <- vapply(1:i, function(x) sum(dnorm(data[,i], theta[2, x], exp(theta[1, x]), log = TRUE)), runif(1))
      if(i < N){
        probY <- c(probY, rep(-Inf, N-i))
      }
      
      
      if(!is.null(logPrior)){
        probY <- probY + logPrior[, i]
      }
      logProb[,i] <- logProb[,i] + 1/S * probY
      
    }
  }
  logProb
}


#' Calculate the predictive density of each vehicle given the DPM approximation
#' 
#' This will estimate (via Monte Carlo) the predictive density implied by each theta cluster, then apply this to each vehicle via sampled values of k
#' @param data A matrix of observations, where columns correspond to different units (i) and rows to time period (t)
#' @param lambda An 5 x N matrix of parameter values. Each column corresponds to the parameters of a single theta cluster. 
#' The first row is the mean of log sigma, the second row is the mean of mu,
#' while the next three rows correspond to the three entires of the lower (cholesky) triangle of the variance matrix
#' @param logPrior A N x N matrix of log prior probabilities for the indicator variables, used for updates. This is the log of the entires of the q(k | y) distribution
#' @param alphaDPM The mass parameter of the Dirichlet Process
#' @param S The number of samples to average over
#' @param support A vector of values to calculate the density over.
DPMindivDens <- function(data, lambda, logPrior = NULL, alphaDPM = 1, S = 100, support = seq(-6, 6, length.out = 1000)){
  
  # One column per unique theta value
  groupDens <- matrix(0, length(support), ncol(data))
  # One column per vehicle, as a weigthed average (according to k samples) of groupDens columns
  indivDens <- matrix(0, length(support), ncol(data))
  
  for(s in 1:S){
    
    theta <- matrix(0, 2, N)
    for(i in 1:N){
      mean <- lambda[1:2, i]
      l <- matrix(c(lambda[3:4, i], 0, lambda[5, i]), 2)
      sigma <- l %*% t(l)
      theta[,i] <- mvtnorm::rmvnorm(1, mean, sigma)
    }
    k <- 1
    nk <- c(1, rep(0, N-1))
    active <- 1
    
    for(i in 2:N){
      probY <- vapply(c(active, i), function(x) sum(dnorm(data[,i], theta[2, x], exp(theta[1, x]), log = TRUE)), runif(1))
      if(!is.null(logPrior)){
        probY <- probY + logPrior[c(active, i), i]
        postProb <- probY - max(probY)
        postProb <- exp(postProb) / sum(exp(postProb))
      } else {
        probY <- probY + log(c(nk[active], alphaDPM) / (alphaDPM + i - 1))
        postProb <- probY - max(probY)
        postProb <- exp(postProb) / sum(exp(postProb))
      }
      newK <- sample(c(active, i), 1, prob = postProb)
      k <- c(k, newK)
      nk[newK] <- nk[newK] + 1
      active <- unique(k)
    }
    
    for(i in 1:N){
      groupDens[,i] <- dnorm(support, theta[2, i], exp(theta[1, i]))
    }
    
    for(i in 1:N){
      indivDens[,i] <- indivDens[,i] + groupDens[,k[i]] / S
    }
  }
  indivDens
}

#' Coordinate Ascent Variational Inference for the DPM model
#' 
#' This function implements the MFVB approximation to the DPM.
#' @param y Matrix of data
#' @param priorMu Base distribution prior mean
#' @param priorSig Base distribution prior variance
#' @param alphaDPM DPM mass parameter
#' @param K number of DPM clusters
#' @param iter Maximum number of coordinate ascent iterations
CAVI <- function(y, priorMu, priorSig, alphaDPM, K, iter = 100){
  N <- ncol(y)
  T <- nrow(y)
  
  muNaught <- priorMu
  lambdaNaught <- priorSig
  
  sample <- rlnorm(100000, priorMu, sqrt(priorSig))
  sigmaMean <- mean(sample)
  sigmaVar <- var(sample)
  
  alphaNaught <- sigmaMean^2 / sigmaVar + 2
  C <- -log(sum(1/sample)) - 1/100000 * sum(log(sample)) 
  for(i in 1:100){
    alphaNaught <- distr::igamma(log(100000 * alphaNaught) + C)
  }
  kappaNaught <- 100000 * alphaNaught / sum(1/sample)
  
  rho <- matrix(1/K, nrow = N, ncol = K)
  alphaHat <- rep(5, K)
  betaHat <- rep(5, K)
  muHat <- colMeans(y[,1:K])
  lambdaHat <- rep(0.5, N)
  ahat <- runif(K)
  bhat <- runif(K)
  
  ySums <- matrix(0, 2, N)
  ySums[1, ] <- colSums(y)
  ySums[2, ] <- colSums(y^2)
  
  for(i in 1:iter){
    
    oldParams <- list(rho = rho, alphaHat = alphaHat, betaHat = betaHat, muHat = muHat, lambdaHat = lambdaHat, ahat = ahat, bhat = bhat)

    # Beta
    for(k in 1:K){
      ahat[k] <- 1 + sum(rho[,k])
      if(k == K){
        bhat[k] <- alphaDPM
      } else {
        bhat[k] <- alphaDPM + sum(rho[,(k+1):K])
      }
    }
    
    # K
    elogsig <- rep(0, K)
    elogbeta <- rep(0, K)
    cumulative <- 0
    for(k in 1:K){
      simsig <- 1/rgamma(1000, alphaHat[k], betaHat[k])
      elogsig[k] <- mean(log(simsig))
      if(k == 1){
        elogbeta[k] <- digamma(ahat[k]) - digamma(ahat[k] + bhat[k])
        cumulative <- cumulative + digamma(bhat[k]) - digamma(ahat[k] + bhat[k])
      } else {
        elogbeta[k] <- digamma(ahat[k]) - digamma(ahat[k] + bhat[k]) + cumulative
        cumulative <- cumulative + digamma(bhat[k]) - digamma(ahat[k] + bhat[k])
      }
    }
    for(j in 1:N){
      loglik <- rep(0, K)
      for(k in 1:K){
        loglik[k] <- -T/2 * elogsig[k] - alphaHat[k] / (2 * betaHat[k]) * (ySums[2, j] - 2 * muHat[k] * ySums[1, j] + T * (muHat[k]^2 + lambdaHat[k]))
      }
      loglik <- loglik + elogbeta
      loglik <- loglik - max(loglik)
      rho[j,] <- exp(loglik) / sum(exp(loglik))
    }
    
    # Mu Stars
    for(k in 1:K){
      meanNumer <- lambdaNaught * sum(rho[,k] * alphaHat[k] / betaHat[k] * ySums[1,]) + muNaught
      meanDenom <- lambdaNaught * T * alphaHat[k] / betaHat[k] * sum(rho[,k]) + 1
      muHat[k] <- meanNumer / meanDenom
      lambdaHat[k] <- lambdaNaught / meanDenom
    }
    
    # Sigma Stars
    for(k in 1:K){
      alphaHat[k] <- alphaNaught + T/2 * sum(rho[,k])
      betaHat[k] <- kappaNaught + 1/2 * sum(rho[,k] * (ySums[2, ] + T * (muHat[k]^2 + lambdaHat[k]) - 2 * muHat[k] * ySums[1, ]))
    }
    
    diffRho <- mean(rho - oldParams$rho)
    diffAlpha <- mean(alphaHat - oldParams$alphaHat)
    diffBeta <- mean(betaHat - oldParams$betaHat)
    diffMu <- mean(muHat - oldParams$muHat)
    diffLambda <- mean(lambdaHat - oldParams$lambdaHat)
    diffA <- mean(ahat - oldParams$ahat)
    diffB <- mean(bhat - oldParams$bhat)
    
    if(abs(diffMu) < 1e-5 & abs(diffLambda) < 1e-5 & abs(diffA) < 0.01 & abs(diffB) < 0.1){
      break
    }
  }
  if(i == iter){
    print(paste('Failed to converge after', iter, 'iterations'))
  } else {
    print(paste('Converged at iteration', i))
  }
  list(rho = rho, alphaHat = alphaHat, betaHat = betaHat, muHat = muHat, lambdaHat = lambdaHat, ahat = ahat, bhat = bhat)
}

