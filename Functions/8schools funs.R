
# Basic linear interpolation
linearInterpolate <- function(x0, x1, y0, y1, x){
  y0 + (x - x0) * (y1 - y0) / (x1 - x0)
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


