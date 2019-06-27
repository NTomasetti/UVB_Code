#MCMC for the Mixture Normal Clustering Model
# Approach uses metropolis hastings for the log variance parameters, and gibbs sampling for the indicator variables and mean parameters
MixNormMCMC <- function(y, reps, drawT, drawK, hyper, thin, stepsize, mix = 2, alpha = rep(1, mix), suppressProgress = FALSE){
  N <- ncol(y)
  T <- nrow(y)
  
  
  accept <- 0
  nSave <- floor(reps/thin)
  saveTheta <- matrix(0, nSave, 2*mix)
  saveK <- matrix(0, nSave, N)
  
  stepsizeCons <- 0.44 * (1 - 0.44)
  
  for(i in 1:reps){
    
    # Runtime timing
    if(i == 50){
      startTime <- Sys.time()
    } else if(i == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    }
    # print progress
    if(!suppressProgress & i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else if(mins > 1){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      } else {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', ceiling(mins * 60), ' seconds.'))
      }
    }
    
    candidate <- drawT
    candidate[1:mix] <- candidate[1:mix] + stepsize * rnorm(mix)
    
    canDens <- MNLogDens(y, candidate, drawK, hyper$mean, hyper$varInv, mix)
    oldDens <- MNLogDens(y, drawT, drawK, hyper$mean, hyper$varInv, mix)
    
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      drawT <- candidate
      oldDens <- canDens
      stepsize <- stepsize + c * (1 - 0.44) / (18 + i)
    } else {
      stepsize <- stepsize - c * 0.44 / (18 + i)
    }
    
    Nj <- table(c(0:(mix-1), drawK)) - 1
    ySums <- rep(0, mix)
    for(j in 1:N){
      ySums[drawK[j]] <- ySums[drawK[j]] + sum(y[,j])
    }
    meanNumer <- exp(drawT[1:mix]) * hyper$mean[mix + 1:mix] + diag(hyper$var)[mix + 1:mix]  * ySums
    meanDenom <- exp(drawT[1:mix]) + Nj * T * diag(hyper$var)[mix + 1:mix]
    var <- diag(hyper$var)[mix + 1:mix] * exp(drawT[1:mix]) / meanDenom
    drawT[mix + 1:mix] <- rnorm(mix, meanNumer / meanDenom, sqrt(var))
    
    # Sample K for this draw
    for(j in 1:N){
      pk <- postK(y[,j], drawT, alpha)
      pk <- pk - max(pk)
      pk <- exp(pk) / sum(exp(pk))
      drawK[j] <- which(rmultinom(1, 1, pk) == 1) - 1
    }
    if(i %% thin == 0){
      saveTheta[i/thin, ] <- drawT
      saveK[i/thin, ] <- drawK
    }
  }
  return(list(theta = saveTheta,
              K = saveK,
              acceptRate = accept/reps,
              stepsize = stepsize))
  
  
}

# Fit the score gradient to the mixture normal model
fitVBScoreMN <- function(data, lambda, mixPrior, mixApprox, probK, dimTheta = 4,
                         S = 50, maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01, ...){
  if(!is.matrix(lambda)){
    lambda <- matrix(lambda, ncol = 1)
  }
  dimLambda <- nrow(lambda)
  N <- ncol(data)
  T <- nrow(data)
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- rep(0, dimLambda)
  V <- rep(0, dimLambda)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold | iter < 100){
    if(iter > maxIter){
      break
    }
    grad <- matrix(0, dimLambda, S)
    h <- matrix(0, dimLambda, S)
    eval <- numeric(S)
    z <- lambda[2*dimTheta*mixApprox + 1:mixApprox]
    approxWeights <- exp(z) / sum(exp(z))

    for(s in 1:S){
      component <- sample(1:mixApprox, 1, prob = approxWeights)
      mu <- lambda[(component-1)*dimTheta + 1:dimTheta]
      sd <- exp(lambda[dimTheta*mixApprox + (component-1)*dimTheta + 1:dimTheta])
      theta <- mu + sd * rnorm(dimTheta)
      
       derivs <- mixtureNormal(data = data,
                              lambda = lambda,
                              probK = probK,
                              theta = theta,
                              mix = mixApprox, 
                              ...)
      
       
       eval[s] <- derivs$val
       grad[,s] <- derivs$grad
       h[,s] <- derivs$grad / derivs$val
    }
    a <- vapply(1:dimLambda, function(x) cov(grad[x, ], h[x,]) / var(h[x,]), 1)
    a[is.na(a)] <- 0
    gradient <- rowMeans(grad - a * h, na.rm = TRUE) 
    gradientSq <- rowMeans((grad - a * h)^2, na.rm = TRUE)
    LB[iter] <- mean(eval, na.rm = TRUE)
    
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    if(any(is.na(alpha * Mst / sqrt(Vst + e)))){
      print('Break - NA gradient value')
      return(matrix(c(lambda, Mst, Vst), ncol = 3))
    }
    lambda <- lambda + alpha * Mst / sqrt(Vst + e)
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 100 == 0){
      print(paste0('iter: ', iter, ' ELBO: ', LB[iter]))
    }
    iter <- iter + 1
  }
  print(paste0('iter: ', min(iter-1, maxIter), ' ELBO: ', LB[min(iter-1, maxIter)]))
  return(list(lambda=lambda, LB = LB, iter = iter))
}
