#' Score VB gradient ascent
#' 
#' This function runs a Variational Bayes gradient ascent algorithm to find the parameters of an approximating distribution that minimises the KL divergence to the true posterior, where the gradients are given by the Score gradient estimator of Ranganath et al 2014, with control variates. It applies a diagonal covariance multivariate normal approximation, and may apply a mixture of diagonal multivariate normals with a suitable specification of lambda.
#' Given a model function that calculates the derivative of the ELBO with resepct to lambda, it will return a list containing the optimal lambda, history of ELBO values, and iterations required to converge for a given model function and dataset.
#' @param lambda Parameter vector or matrix of auxiliary parameters. The first n rows correspond to means of an n dimensional theta, and the second n rows corresponds to the log of the standard deviations. If given as a matrix, it is assumed that a mixture approximation is to be applied, where each column corresponds to parameters for one component of the mixture. In this case an additional row is required at the end of the matrix, which corresponds to the unnormalised mixture weights 'z', where normalised mixture weights 'pi' = exp(z) / sum(exp(z))
#' @model A function with first two arguments lambda, theta, that retursn a list with two elements 'grad', the vector of the gradient of the ELBO with respect to lambda, and 'val', the value of the ELBO for this lambda.
#' @param S The number of Monte Carlo estimates per iteration. Defaults to 50.
#' @param maxIter Interger. Maximum number of gradient ascent iterations. Defaults to 5000.
#' @param alpha adam optimisation control parameter. Defaults to 0.01.
#' @param beta1 adam optimisation control parameter. Defaults to 0.9.
#' @param beta2 adam optimisation control parameter. Defaults to 0.99.
#' @param threshold Maximum difference in mean value of ELBO before convergence is achieved. Defaults to 0.01.
#' @param suppressProgress Boolean, if true the program will not print ELBO values to console during optimisation. Defaults to FALSE
#' @param ... Extra arguments passed into model (e.g., data, hyperparameters)
#' @export


scoreVB <- function(lambda, model, S = 50, maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01, suppressProgress = FALSE, ...){
  if(!is.matrix(lambda)){
    lambda <- matrix(lambda, ncol = 1)
  }
  dimLambda <- nrow(lambda)
  mix <- ncol(lambda)
  if(mix == 1){
    dimTheta = dimLambda / 2
  } else {
    dimTheta = (dimLambda - 1) / 2
  }
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- rep(0, dimLambda)
  V <- rep(0, dimLambda)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  var <- matrix(0, dimLambda, maxIter)
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    grad <- h <- matrix(0, dimLambda, S)
    eval <- numeric(S)
    if(mix != 1){
     z <- lambda[2 * dimTheta + 1, ]
    } else {
     z <- 0
    }
    
    pi <- exp(z) / sum(exp(z))
    s <- 0
  
    epsilon <- matrix(rnorm(dimTheta *S), nrow = S)
    
    for(s in 1:S){
      k <- sample(1:mix, 1, prob=pi)
      Qmean <- lambda[1:dimTheta, k]
      Qsd <- exp(lambda[dimTheta + 1:dimTheta, k])
      if(dimTheta > 1){
        theta <- c(Qmean + Qsd * epsilon[s,])
      } else {
        theta <- Qmean + Qsd * epsilon[s,]
      }
      derivs <- model(lambda, theta, mix = mix, ...)
      eval[s] <- derivs$val
      grad[,s] <- derivs$grad
      h[,s] <- grad[,s] / eval[s]
    }
    a <- vapply(1:dimLambda, function(x) cov(grad[x, ], h[x,], use = 'complete.obs') / var(h[x, ], na.rm = TRUE), runif(1))
    a[is.na(a)] <- 0
    var[,iter] <- apply(grad - a * h, 1, var)
    gradient <- rowMeans(grad - a * h, na.rm = TRUE)
    gradientSq <- rowMeans((grad - a * h)^2, na.rm = TRUE)
    LB[iter] <- mean(eval, na.rm = TRUE)

    
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    if(any(is.na(alpha * Mst / sqrt(Vst + e)))){
      print('Break - NA in gradient vector')
      break
    }
    lambda <- lambda + alpha * Mst / sqrt(Vst + e)
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 100 == 0 & !suppressProgress){
      print(paste0('iter: ', iter, ' ELBO: ', LB[iter]))
    }
    iter <- iter + 1
  }
  iter <- iter - 1
  if(!suppressProgress){
    print(paste0('iter: ', iter, ' ELBO: ', meanLB))
  }
  return(list(lambda=lambda,
              LB = LB[1:iter], 
              iter = iter))
}
