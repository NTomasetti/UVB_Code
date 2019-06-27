#' Block Random Walk Metrolpolis Hastings Sampler
#'
#' Jointly samples a candidate from a Random Walk Metrolpolis Hastings algorithm with a given acceptance rate to be targeted.
#' @param start A vector of starting values for the parameters. If a named vector is inputted, the names will be attached to the output
#' @param iter The total number of iterations,
#' @param model A function with first argument theta that evaluates the log joint density at that value of theta
#' @param targetAcceptance The targetted MH acceptance rate, defaults to 0.234
#' @param burnIn The number of discarded iterations, should be lower than iter. Defaults to 0.
#' @param thin The thinning factor. Draws will only be saved when the iteration number equals zero, modulus thin. Defaults to 1.
#' @param stepsize The starting standard deviation for the Random Walk sampler. Will change dynamically through the algorithm
#' @param suppressProgress boolean. If true, supress estimated time remaining messages that appear each 1000 iterations. Defaults to false.
#' @param ... Extra arguments passed to model, eg. Data, Prior parameters
#' @export
metropolisHastings <- function(start, iter, model, targetAcceptance = 0.234, 
                               burnIn = 0, thin = 1, stepsize = 0.01, suppressProgress = FALSE,  ...){

  # Set up matrix of saved draws
  dim <- length(start)
  nSave <- floor((iter - burnIn) / thin)
  saveDraws <- matrix(0, nSave, dim)

  # changing MH acceptance rate
  alpha <- - qnorm(targetAcceptance/2)
  stepsizeCons <- (1 - 1/dim) * sqrt(2*pi) * exp(alpha^2/2) / (2 * alpha)  +  1 / (dim * 0.234 * (1 - 0.234))

  # Starting draw density
  draw <- start
  oldDens <- model(draw, ...)
  for(i in 1:iter){

    # Runtime timing
    if(i == 50){
      startTime <- Sys.time()
    } else if(i == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (iter - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    }
    # print progress
    if(!suppressProgress & i %% 1000 == 0){
      mins <-  (iter - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else if(mins > 1){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      } else {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', ceiling(mins * 60), ' seconds.'))
      }
    }

    # Metropolis Hastings Iteration
    candidate <- draw + stepsize * rnorm(dim)
    canDens <- model(candidate, ...)
    ratio <- exp(canDens - oldDens)
    c <- stepsize * stepsizeCons
    if(runif(1) < ratio){
      draw <- candidate
      oldDens <- canDens
      stepsize <- stepsize + c * (1 - targetAcceptance) / (28 + i)
    } else {
      stepsize <- stepsize - c * targetAcceptance / (28 + i)
    }
    if(i > burnIn & i %% thin == 0){
      saveDraws[(i - burnIn)/thin, ] <- draw
    }
  }

  draws <- data.frame(saveDraws)
  if(is.null(names(start))){
    colnames(draws) <- paste0('V', 1:dim)
  } else {
    colnames(draws) <- names(start)
  }
  draws$iter <- seq(burnIn + thin, iter, thin)
  draws
}






