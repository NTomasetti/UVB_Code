
# Basic linear interpolation
linearInterpolate <- function(x0, x1, y0, y1, x){
  y0 + (x - x0) * (y1 - y0) / (x1 - x0)
}

ISschools <- function(lambda, qScore, samples, dSamples, schoolsNew, batch = 1, maxIter = 5000,
                      alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01){
  if(!is.matrix(lambda)){
    lambda <- matrix(lambda, ncol = 1)
  }
  dimLambda <- nrow(lambda)
  if(!is.matrix(samples)){
    # samples may be passed in as a vector if theta is one-dimensional
    samples <- matrix(samples, ncol = 1)
  }
  diff <- threshold + 1
  iter <- 1
  LB <- ESS <- numeric(maxIter)
  M <- V <- numeric(dimLambda)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  S <- nrow(samples)
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    
    eval <- numeric(S)
    dLq <- matrix(0, S, dimLambda)
    lQ <- rep(0, S)
    
    # Simulate theta for the new schools
    newTheta <- matrix(0, S, batch)
    logPrior <- rep(0, S)
    for(s in 1:S){
      newTheta[s, ] <- rnorm(batch, lambda[dimLambda/2-((batch-1):0)], exp(lambda[dimLambda - ((batch-1):0)]))
      logPrior[s] <- sum(dnorm(newTheta[s, ], samples[s, 2], exp(samples[s, 1]), log = TRUE))
    }
  
    # Calculate the log density and score of the approximating distribution across all schools
    
    for(s in 1:S){
      model <- qScore(c(samples[s,], newTheta[s, ]), lambda)
      lQ[s] <- model$val
      dLq[s,] <- model$grad
    }
    
    # Calculate the importance sampling weights.
    
    w <- exp(lQ - dSamples - logPrior)
    
    # Calculate the log-likelihood and prior of the new samples
    
    loglik <- rep(0, S)
    for(s in 1:S){
      loglik[s] <- sum(dnorm(schoolsNew[, 1], newTheta[s, ], schoolsNew[, 2], log = TRUE))
    }
    logjoint <- loglik + logPrior + dSamples
    logjoint <- logjoint - max(logjoint)
    
    eval <- w * (logjoint - lQ)
    grad <- w * dLq * (logjoint - lQ) 
    f <- w * dLq
    
    a <-  diag(cov(grad, f) / var(f))
    a[is.na(a)] <- 0
    
    gradient <- colMeans(grad - a * f)
    gradientSq <- colMeans((grad - a * f)^2)
    LB[iter] <- mean(eval)
    
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mstar <- M / (1 - beta1^iter)
    Vstar <- V / (1 - beta2^iter)
    update <- alpha * Mstar / (sqrt(Vstar) + e)
    if(any(is.na(update))){
      print('Break')
      break
    }
    lambda <- lambda + update
    
    wNorm <- w / sum(w)
    ESS[iter] <- 1 / sum(wNorm^2)
    
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 100 == 0){
       print(paste0('Iteration: ', iter, ' ELBO: ', LB[iter]))
    }
    iter <- iter + 1
  }
  iter <- iter - 1
  print(paste0('iter: ', iter, ' ELBO: ', LB[iter]))
  return(list(lambda=lambda,
              LB = LB[1:iter],
              iter = iter,
              ESS = ESS[1:iter]))
}

# Calculate KL divergence from q to p for each posterior marginal
# Inputs VBdens as given by gaussianDensity.R, and MCMCdens in the same format
# ie. both require names, support and density columns
KLdiv <- function(VBdens, MCMCdens){
  dim <- length(unique(VBdens$var))
  names <- unique(VBdens$var)
  results <- rep(0, dim)
  
  for(i in 1:dim){
    sub <- VBdens %>% filter(var == names[i])
    subMC <- MCMCdens %>% filter(var == names[i])
    prob <- sub$density / sum(sub$density)
    theta <- data.frame(t = sample(sub$support, 1000, TRUE, prob))
    theta %>% 
      left_join(sub, by = c('t' = 'support')) %>%
      rename(q = density) %>%
      left_join(subMC, by = c('t' = 'support')) %>%
      rename(p = density) %>%
      summarise(KL = mean(log(q / p))) %>% 
      .$KL-> results[i]
  }
  data.frame(var = names, kl = results)
}

# A wrapper for density
calcDens <- function(df){
  density(df$draw, n = 1000)
}


