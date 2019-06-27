#' Calculates the marginal density of a transformed multivariate gaussian.
#'
#' This function returns a dataframe with columns for the support of each marginal distribution, the density at that value, and the name of that variable if provided.
#' @param mu the length n mean vector of the multivariate gaussian
#' @param sd the length n vector of standard deviations.
#' @param transform a length n vector with the marginal transform of each variable from c('identity', 'exp', 'sigmoid', stretchedSigmoid'). If null, defaults each variable to 'identity'
#' @param names a length n vector of the names of each variable. If null, defaults to V1, V2, .., Vn
#' @param supports a matrix with n columns that gives the support of each variable (per column) to calculate the densities over. If NULL, automatically find a suitable support range.
#' @export
gaussianDensity <- function(mu, sd, transform = NULL, names = NULL, supports = NULL){
  if(!is.vector(mu)){
    stop('mu must be a vector')
  }
  if(!is.vector(sd)){
    stop('sd must be a vector')
  }
  n <- length(mu)
  if(length(sd) != n){
    stop('mu and sd have different lengths')
  }
  if(is.null(transform)){
    transform <- rep('identity', n)
  } else {
    if(!is.vector(transform)){
      stop('transform must be a vector')
    }
    if(length(transform) != n){
      stop('mu and transform have different lengths')
    }
  }
  if(!all(transform %in% c('identity', 'exp', 'sigmoid', 'stretchedSigmoid'))){
    stop('invalid transformation')
  }
  if(is.null(names)){
    names <- paste0('V', 1:n)
  } else {
    if(!is.vector(names)){
      stop('names must be a vector')
    }
    if(length(names) != n){
      stop('mu and names have different lengths')
    }
    if(!is.character(names)){
      names <- as.character(names)
    }
  }
  if(!is.null(supports)){
    if(n == 1 & !is.matrix(supports) & !is.vector(supports)){
      stop('supports must be a matrix')
    } else if(!is.matrix(supports) & n != 1){
      stop('supports must be a matrix')
    }
    if(n > 1){
      if(ncol(supports) < n){
        stop('supports must have n columns')
      }
    }
  }

  dens <- data.frame()
  for(i in 1:n){
    if(transform[i] == 'identity') {
      if(is.null(supports)){
        support <- seq(mu[i] - 5*sd[i], mu[i] + 5*sd[i], length.out=1000)
      } else {
        support <- supports[,i]
      }
      density <- dnorm(support, mu[i], sd[i])
    } else if(transform[i] == 'exp'){
      if(is.null(supports)){
        mean <- exp(mu[i] + 0.5 * sd[i]^2)
        stdev <- sqrt((exp(sd[i]^2) - 1)*exp(2*mu[i]+sd[i]^2))
        support <- seq(max(1e-08, mean - 5*stdev), mean+5*stdev, length.out=1000)
      } else {
        support <- supports[,i]
      }
      density <- dlnorm(support, mu[i], sd[i])
    } else if (transform[i] == 'sigmoid') {
      if(is.null(supports)){
        sample <- 1 / (1 + exp(-rnorm(1000, mu[i], sd[i])))
        mean <- mean(sample)
        stdev <- sd(sample)
        support <- seq(max(0.001, mean-5*stdev), min(0.999, mean+5*stdev), length.out=1000)
      } else {
        support <- supports[,i]
      }
      density <- dnorm(log(support / (1-support)), mu[i], sd[i]) / (support - support^2)
    } else if(transform[i] == 'stretchedSigmoid'){
      if(is.null(supports)){
        sample <- 2 / (1 + exp(-rnorm(1000, mu[i], sd[i]))) - 1
        mean <- mean(sample)
        stdev <- sd(sample)
        support <- seq(max(-0.999, mean-5*stdev), min(0.999, mean+5*stdev), length.out=1000)
      } else {
        support <- supports[,i]
      }
      density <- dnorm(-log(2/(support+1)-1), mu[i], sd[i]) * 2 / (2*(support+1) - (support+1)^2)
    }
    df <- data.frame(support = support, density = density, var = names[i])
    dens <- rbind(dens, df)
  }
  dens
}
