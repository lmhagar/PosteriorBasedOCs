## functions for multinomial numerical study, where Algorithm 3 is used
## instead of Algorithm 1 (i.e., when the MLEs are not approximately normal)

## BEGIN SETUP ##

## load necessary packages
require(qrng)

## this helper function implements Algorithm 3 with the illustrative example to
## return the logit of the posterior probability (used with the function FindGamma())
targetPower <- function(u, params, delta, n_val, q, hyper){
  
  if (n_val <= 0){return(-1.5)}
  
  ## convert parameter values on p-scale to the Z-scale for each group
  p.1.1 <- params[1]; p.1.2 <- params[2]
  p.1.3 <- params[3]; p.1.4 <- params[4]
  z.1.2 <- p.1.2/(1-p.1.1)
  z.1.3 <- p.1.3/(1-p.1.1-p.1.2)
  z.1.4 <- p.1.4/(1-p.1.1-p.1.2-p.1.3)

  p.2.1 <- params[5]; p.2.2 <- params[6]
  p.2.3 <- params[7]; p.2.4 <- params[8]
  z.2.2 <- p.2.2/(1-p.2.1)
  z.2.3 <- p.2.3/(1-p.2.1-p.2.2)
  z.2.4 <- p.2.4/(1-p.2.1-p.2.2-p.2.3)
  
  dat1 <- rep(0, 5); dat2 <- rep(0, 5)
  
  q1 <- qbinom(u[1], n_val, p.1.1)
  c1 <- pbinom(q1 - 1, n_val, p.1.1); c2 <- pbinom(q1, n_val, p.1.1)
  
  ## the CDF of X^* is piecewise
  dat1[1] <- q1 - ifelse(q1 > 0, 0.5, 0) + ifelse(q1 > 0 & q1 < n_val, 1, 0.5)*(u[1] - c1)/(c2 - c1)
  
  ## f1 is the ceiling of the number of observations that we must still allocate
  f1 <- n_val - floor(dat1[1]); q1 <- qbinom(u[2], f1, z.1.2)
  c1 <- pbinom(q1 - 1, f1, z.1.2); c2 <- pbinom(q1, f1, z.1.2)
  
  ## the second step here applies a proportion decrease to account for there being a
  ## noninteger number of successes remaining
  dat1[2] <- q1 - ifelse(q1 > 0, 0.5, 0) + ifelse(q1 > 0 & q1 < f1, 1, 0.5)*(u[2] - c1)/(c2 - c1)
  dat1[2] <- dat1[2]*(n_val - dat1[1])/f1
  dat1[2] <- ifelse(is.na(dat1[2]) | !(is.finite(dat1[2])), 0, dat1[2])
  
  ## this process is repeated for each category
  f1 <- n_val - floor(sum(dat1[1:2])); q1 <- qbinom(u[3], f1, z.1.3)
  c1 <- pbinom(q1 - 1, f1, z.1.3); c2 <- pbinom(q1, f1, z.1.3)

  dat1[3] <- q1 - ifelse(q1 > 0, 0.5, 0) + ifelse(q1 > 0 & q1 < f1, 1, 0.5)*(u[3] - c1)/(c2 - c1)
  dat1[3] <- dat1[3]*(n_val - sum(dat1[1:2]))/f1
  dat1[3] <- ifelse(is.na(dat1[3]) | !(is.finite(dat1[3])), 0, dat1[3])
  
  f1 <- n_val - floor(sum(dat1[1:3])); q1 <- qbinom(u[4], f1, z.1.4)
  c1 <- pbinom(q1 - 1, f1, z.1.4); c2 <- pbinom(q1, f1, z.1.4)
  
  dat1[4] <- q1 - ifelse(q1 > 0, 0.5, 0) + ifelse(q1 > 0 & q1 < f1, 1, 0.5)*(u[4] - c1)/(c2 - c1)
  dat1[4] <- dat1[4]*(n_val - sum(dat1[1:3]))/f1
  dat1[4] <- ifelse(is.na(dat1[4]) | !(is.finite(dat1[4])), 0, dat1[4])
  
  dat1[5] <- n_val - sum(dat1[1:4])
  
  ## repeat for group 2
  N1 <- n_val
  N2 <- round(q*n_val)
  q1 <- qbinom(u[5], N2, p.2.1)
  c1 <- pbinom(q1 - 1, N2, p.2.1); c2 <- pbinom(q1, N2, p.2.1)
  
  dat2[1] <- q1 - ifelse(q1 > 0, 0.5, 0) + ifelse(q1 > 0 & q1 < N2, 1, 0.5)*(u[5] - c1)/(c2 - c1)
  
  f1 <- N2 - floor(dat2[1]); q1 <- qbinom(u[6], f1, z.2.2)
  c1 <- pbinom(q1 - 1, f1, z.2.2); c2 <- pbinom(q1, f1, z.2.2)
  
  dat2[2] <- q1 - ifelse(q1 > 0, 0.5, 0) + ifelse(q1 > 0 & q1 < f1, 1, 0.5)*(u[6] - c1)/(c2 - c1)
  dat2[2] <- dat2[2]*(N2 - dat2[1])/f1
  dat2[2] <- ifelse(is.na(dat2[2]) | !(is.finite(dat2[2])), 0, dat2[2])
  
  f1 <- N2 - floor(sum(dat2[1:2])); q1 <- qbinom(u[7], f1, z.2.3)
  c1 <- pbinom(q1 - 1, f1, z.2.3); c2 <- pbinom(q1, f1, z.2.3)
  
  dat2[3] <- q1 - ifelse(q1 > 0, 0.5, 0) + ifelse(q1 > 0 & q1 < f1, 1, 0.5)*(u[7] - c1)/(c2 - c1)
  dat2[3] <- dat2[3]*(N2 - sum(dat2[1:2]))/f1
  dat2[3] <- ifelse(is.na(dat2[3]) | !(is.finite(dat2[3])), 0, dat2[3])
  
  f1 <- N2 - floor(sum(dat2[1:3])); q1 <- qbinom(u[8], f1, z.2.4)
  c1 <- pbinom(q1 - 1, f1, z.2.4); c2 <- pbinom(q1, f1, z.2.4)
  
  dat2[4] <- q1 - ifelse(q1 > 0, 0.5, 0) + ifelse(q1 > 0 & q1 < f1, 1, 0.5)*(u[8] - c1)/(c2 - c1)
  dat2[4] <- dat2[4]*(N2 - sum(dat2[1:3]))/f1
  dat2[4] <- ifelse(is.na(dat2[4]) | !(is.finite(dat2[4])), 0, dat2[4])
  
  dat2[5] <- N2 - sum(dat2[1:4])
  
  ## extract hyperparameters
  alphas1 <- hyper[,1]
  betas1 <- hyper[,2]
  alphas2 <- hyper[,3]
  betas2 <- hyper[,4]
  
  ## obtain the posterior modes on the logit(Z) scale along with the variances
  ## corresponding to the Laplace approximation to the posterior (the variables
  ## on the (logit)Z-scale are independent)
  for (j in 1:4){
    assign(paste0("z1",j),
           (dat1[j] + alphas1[j])/(sum(dat1[j:5]) + alphas1[j] + betas1[j]))
    
    assign(paste0("sigma1",j), 1/sqrt((sum(dat1[j:5]) + alphas1[j] + betas1[j])*get(paste0("z1",j))*(1-get(paste0("z1",j)))))
  }
  
  for (j in 1:4){
    assign(paste0("z2",j),
           (dat2[j] + alphas2[j])/(sum(dat2[j:5]) + alphas2[j] + betas2[j]))
    
    assign(paste0("sigma2",j), 1/sqrt((sum(dat2[j:5]) + alphas2[j] + betas2[j])*get(paste0("z2",j))*(1-get(paste0("z2",j)))))
  }
  
  ## convert from logit(Z) to p-scale
  p12 <- z12*(1-z11); p13 <- z13*(1-z11)*(1-z12)
  p14 <- z14*(1-z11)*(1-z12)*(1-z13)
  
  p22 <- z22*(1-z21); p23 <- z23*(1-z21)*(1-z22)
  p24 <- z24*(1-z21)*(1-z22)*(1-z23)
  
  ## use the delta method to get the covariance matrix on the p-scale (not independent)
  deriv1_mat <- z11*sigma11*c(1 - z11, -p12, -p13, -p14)
  deriv1_mat <- rbind(deriv1_mat, z12*sigma12*c(0, (1-z11 -p12), -p13, -p14))
  deriv1_mat <- rbind(deriv1_mat, z13*sigma13*c(0, 0, (1-z11 -p12 - p13), -p14))
  deriv1_mat <- rbind(deriv1_mat, sigma14*c(0, 0, 0, p14*(1-z14)))
  
  deriv2_mat <- z21*sigma21*c(1 - z21, -p22, -p23, -p24)
  deriv2_mat <- rbind(deriv2_mat, z22*sigma22*c(0, (1-z21 -p22), -p23, -p24))
  deriv2_mat <- rbind(deriv2_mat, z23*sigma23*c(0, 0, (1-z21 -p22 - p23), -p24))
  deriv2_mat <- rbind(deriv2_mat, sigma24*c(0, 0, 0, p24*(1-z24)))
  
  ## covariance matrix for each group
  Fish1 <- t(deriv1_mat)%*%deriv1_mat
  Fish2 <- t(deriv2_mat)%*%deriv2_mat
  
  ## define theta metrics in terms of the transformed posterior modes
  theta1 <- z11 + 2*p12 + 3*p13 + 4*p14 + 5*(1 - z11 - p12 - p13 - p14)
  theta2 <- z21 + 2*p22 + 3*p23 + 4*p24 + 5*(1 - z21 - p22 - p23 - p24)
  
  ## compute partial derivatives for log(theta + 4) - log(theta - 4), which ensures
  ## that we have support over the entire real line
  star1 <- (4*(z21 - z11) + 3*(p22 - p12) + 2*(p23 - p13) + (p24 - p14) + 4)^(-1)
  star2 <- (4 - 4*(z21 - z11) - 3*(p22 - p12) - 2*(p23 - p13) - (p24 - p14))^(-1)
  
  ## partial derivatives for the first group for another application of the delta method
  ## derivatives for group 2 are just -1 times derivatives for group 1
  dp11 <- -4*(star1 + star2)
  dp12 <- -3*(star1 + star2)
  dp13 <- -2*(star1 + star2)
  dp14 <- -1*(star1 + star2)
  
  ## apply the delta method to get the limiting variance for each theta_j metric
  avar1 <- t(c(dp11, dp12, dp13, dp14))%*%Fish1%*%c(dp11, dp12, dp13, dp14)
  avar2 <- t(c(dp11, dp12, dp13, dp14))%*%Fish2%*%c(dp11, dp12, dp13, dp14)
  
  ## apply the delta method to get the limiting variance for theta = theta1 - theta2
  avar <- avar1 + avar2
  
  ## compute posterior probability given this normal approximation
  realP <- pnorm(log(delta[2] + 4) - log(4 - delta[2]), log(theta1 - theta2 + 4) - log(4 - theta1 + theta2), sqrt(avar)) - pnorm(log(delta[1] + 4) - log(4 - delta[1]), log(theta1 - theta2 + 4) - log(4 - theta1 + theta2), sqrt(avar))
  
  ## slight perturbation of the posterior probability if it is too close to 0 or 1
  ## this ensures that the logit of the posterior probability is finite.
  if (realP > 1 - 10^(-7)){
    realP <- 1 - 10^(-7)
  }
  else if (realP < .Machine$double.eps){
    realP <- .Machine$double.eps
  }
  
  ## return the logit of the posterior probability
  return(log(realP) - log(1 - realP))
}

## function where we return what we need to get the contour plots
## Code to implement Algorithm 2 for the illustrative example. We explore sample sizes
## and return the optimal design (i.e., the (n, gamma) combination)
findGamma <- function(green_par, red_par, pwr, typeI, deltas, q, 
                      alphas1, betas1, alphas2, betas2, m = 8192, m0 = 512,
                      seed_green = 1, seed_red = 2, upper_init = 1000, prng = FALSE, 
                      contour = FALSE){
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option available for comparison
  if (prng == FALSE){
    sob_green <- sobol(m, d = 9, randomize = "digital.shift", seed = seed_green)
    sob_red <- sobol(m, d = 9, randomize = "digital.shift", seed = seed_red)
  }
  else{
    set.seed(seed_green); sob_green <- matrix(runif(9*m), ncol = 9)
    
    set.seed(seed_red); sob_red <- matrix(runif(9*m), ncol = 9)
  }
  
  ## index the eta values from the design priors such that the anticipated value for
  ## theta increases with the index; we do this using the first dimension of the Sobol' sequence
  green_means1 <- 5 - 4*green_par[,1] - 3*green_par[,2] - 2*green_par[,3] - green_par[,4]
  green_means2 <- 5 - 4*green_par[,5] - 3*green_par[,6] - 2*green_par[,7] - green_par[,8]
  green_diff <- green_means1 - green_means2
  green_par <- green_par[order(green_diff),]
  green_par <- green_par[ceiling(length(green_diff)*sob_green[,1]),]
  
  ## repeat this process for the red region
  red_means1 <- 5 - 4*red_par[,1] - 3*red_par[,2] - 2*red_par[,3] - red_par[,4]
  red_means2 <- 5 - 4*red_par[,5] - 3*red_par[,6] - 2*red_par[,7] - red_par[,8]
  red_diff <- red_means1 - red_means2
  red_par <- red_par[order(red_diff),]
  red_par <- red_par[ceiling(length(red_diff)*sob_red[,1]),]
  
  if (m < m0){stop("m0 should be less than m")}
  lower <- 1
  n <- upper_init
  upper <- upper_init
  green_vec <- rep(0, m0)
  red_vec <- rep(1, m0)
  
  mid1 <- ceiling((1 - typeI)*m0)
  mid2 <- floor((1-pwr)*m0)
  
  ## find a larger upper bound if the initial upper bound provided
  ## is not large enough
  while (sort(red_vec)[mid1] > sort(green_vec)[mid2]){
    green_vec <- NULL
    red_vec <- NULL
    ## implement Algorithm 3 for each of the first m0 points
    for (i in 1:m0){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    ## increase upper bound if the criterion on the order statistics is not satisfied
    if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
      green_max <- green_vec; red_max <- red_vec
      upper <- n
      n_max <- n
    }
    else{
      lower <- n
      green_min <- green_vec; red_min <- red_vec
      n_min <- n
      n <- 2*n
    }
  }
  
  ## get a lower point to construct linear approximations to the logits
  ## of the posterior probabilities (this is only required if the initial
  ## upper bound for n was sufficiently large)
  if(lower == 1){
    n <- ceiling(0.5*(upper + lower))
    green_vec <- NULL
    red_vec <- NULL
    for (i in 1:m0){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    ## if the criterion is still satisfied then this smaller sample size is still
    ## an upper bound
    if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
      green_min <- green_vec; red_min <- red_vec
      upper <- n
      n_min <- n
    }
    else{
      lower <- n
      green_min <- green_vec; red_min <- red_vec
      n_min <- n
    }
  }
  
  ## construct the linear approximations for the green and red groups
  green_slope <- (green_max - green_min)/(n_max-n_min)
  red_slope <- (red_max - red_min)/(n_max-n_min)
  green_int <- green_min - green_slope*n_min
  red_int <- red_min - red_slope*n_min
  
  ## get an initial sample size using the linear approximations to the posterior 
  ## probabilities (no Algorithm 3)
  upper_temp <- upper
  lower_temp <- lower
  while ((upper_temp - lower_temp) > 1){
    n <- ceiling(0.5*(upper_temp + lower_temp))
    green_vec <- NULL
    red_vec <- NULL
    for (i in 1:m0){
      green_vec <- green_int + green_slope*n
      
      red_vec <- red_int + red_slope*n
    }
    if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
      upper_temp <- n
    }
    else{
      lower_temp <- n
    }
  }
  
  ii <- 0
  n <- max(10,upper_temp)
  n_last <- n
  while (ii < 1 | (abs(n_last - n)/n) >= 0.1 | abs(n_last - n) > 5){
    ii <- ii + 1
    n_last <- n
    green_vec <- NULL
    red_vec <- NULL
    ## we next run Algorithm 3 at that initialize sample size for the first m0 points
    for (i in 1:m0){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    ## these posterior probabilities will be used later to create better linear approximations
    ## on the logit scale (we need to determine whether n is a lower or upper bound here with
    ## respect to the first m0 points)
    if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
      upper <- n
      if (n < n_min){
        n_max <- n_min
        green_max <- green_min; red_max <- red_min
        n_min <- n
        green_min <- green_vec; red_min <- red_vec
      }
      else{
        n_max <- n
        green_max <- green_vec; red_max <- red_vec
      }
    }
    else{
      lower <- n
      if (n > n_min){
        n_min <- n
        green_min <- green_vec; red_min <- red_vec
      }
      else{
        n_max <- n_min
        green_max <- green_min; red_max <- red_min
        n_min <- n
        green_min <- green_vec; red_min <- red_vec
      }
    }
    
    ## now use the probabilities at n and the other bound for binary search to get
    ## better linear approximations on the logit scale. We then explore sample sizes
    ## again with the first m0 points using these linear approximations (no Algorithm 3)
    green_slope <- (green_max - green_min)/(n_max-n_min)
    red_slope <- (red_max - red_min)/(n_max-n_min)
    green_int <- green_min - green_slope*n_min
    red_int <- red_min - red_slope*n_min
    upper_temp <- upper
    lower_temp <- lower
    while ((upper_temp - lower_temp) > 1){
      n <- ceiling(0.5*(upper_temp + lower_temp))
      green_vec <- NULL
      red_vec <- NULL
      for (i in 1:m0){
        green_vec <- green_int + green_slope*n
        
        red_vec <- red_int + red_slope*n
      }
      if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
        upper_temp <- n
      }
      else{
        lower_temp <- n
      }
    }
    n <- upper_temp
    ## we check how much the sample size recommendation changed with respect to
    ## what we started this iteration of the floor loop with; if there is not
    ## much change, then we just take this new recommendation as n(0)
  }
  
  ## compute the posterior probabilities using Algorithm 3 for all m points for
  ## each Sobol' sequence at the initial sample size n(0)
  n_init <- n
  green_vec <- NULL
  red_vec <- NULL
  for (i in 1:m){
    green_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(green_par[i,]),
                                delta = deltas, u = sob_green[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    red_vec[i] <- targetPower(n_val = n, q = q,
                              params = as.numeric(red_par[i,]),
                              delta = deltas, u = sob_red[i,2:9],
                              hyper = cbind(alphas1, betas1, alphas2, betas2))
  }
  
  ## save the logits of the probabilities
  green_vec_n0 <- green_vec
  red_vec_n0 <- red_vec
  n0 <- n
  
  ## these flags denote whether alpha and beta are too small to explore
  ## symmetrically around the relevant order statistic
  flag1 <- ifelse(2*typeI*m < m0, 1, 0)
  flag2 <- ifelse(2*(1 - pwr)*m < m0, 1, 0)
  
  ## if beta is too small, we take the smallest m0 order statistics
  if (flag2){
    low2 <- 1; high2 <- m0
  }
  ## otherwise we explore symmetrically around the order statistic
  else{
    low2 <- floor((1 - pwr)*m) - 0.5*m0 + 1
    high2 <- low2 + m0 - 1
  }
  
  ## if alpha is too small, we take the largest m0 order statistics
  if (flag1){
    low1 <- m - m0 + 1; high1 <- m
  }
  else{
    low1 <- floor((1 - typeI)*m) - 0.5*m0 + 1
    high1 <- low1 + m0 - 1
  }
  
  ## we update the indices of the order statistics to reflect the larger
  ## set of points (m vs. m0)
  mid1 <- ceiling((1 - typeI)*m)
  mid2 <- floor((1-pwr)*m)
  
  ## if we do not satisfy the criterion on the order statistics, we create the linear
  ## approximation using a larger second sample size 
  if (sort(red_vec)[mid1] > sort(green_vec)[mid2]){
    green_low <- green_vec; red_low <- red_vec
    lower <- n
    n_low <- n
    ## choose a larger sample size n(1)
    n <- max(ceiling(1.1*n), n + 5)
    green_vec <- NULL
    red_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    green_vec_n1 <- green_vec
    red_vec_n1 <- red_vec
    n1 <- n
    
    n_high <- n
    green_high <- green_vec; red_high <- red_vec
    
    ## obtain the linear approximations for all m points
    green_slope <- (green_high - green_low)/(n_high-n_low)
    red_slope <- (red_high - red_low)/(n_high-n_low)
    green_int <- green_low - green_slope*n_low
    red_int <- red_low - red_slope*n_low
    
    ## if the criterion on the order statistics is still not satisfied for this
    ## larger sample size, we need to find a larger upper bound for the binary search
    ## we do not do this using Algorithm 3 -- just the linear approximations to the
    ## posterior probabilities on the logit scale.
    if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
      upper <- n
    }
    else {
      n <- ceiling(1.5*n_low)
      green_vec <- green_int + green_slope*n
      red_vec <- red_int + red_slope*n
      upper <- n
      
      while(sort(red_vec)[mid1] > sort(green_vec)[mid2]){
        n <- ceiling(1.5*n)
        green_vec <- green_int + green_slope*n
        red_vec <- red_int + red_slope*n
        if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
          upper <- n
        }
      }
    }
  }
  ## if we satisfy the criterion on the order statistics, we create the linear
  ## approximation using a smaller second sample size 
  else{
    green_high <- green_vec; red_high <- red_vec
    upper <- n
    n_high <- n
    ## choose a smaller sample size n(1)
    n <- min(floor(0.9*n), n - 5)
    green_vec <- NULL
    red_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    green_vec_n1 <- green_vec
    red_vec_n1 <- red_vec
    n1 <- n
    
    n_low <- n
    green_low <- green_vec; red_low <- red_vec
    
    ## obtain the linear approximations for all m points
    green_slope <- (green_high - green_low)/(n_high-n_low)
    red_slope <- (red_high - red_low)/(n_high-n_low)
    green_int <- green_low - green_slope*n_low
    red_int <- red_low - red_slope*n_low
    
    ## if the criterion on the order statistics is still satisfied for this
    ## larger sample size, we need to find a smaller lower bound for the binary search
    if (sort(red_vec)[mid1] > sort(green_vec)[mid2]){
      lower <- n
    }
    else{
      n <- floor(0.5*n_high)
      green_vec <- green_int + green_slope*n
      red_vec <- red_int + red_slope*n
      lower <- n
      
      while(sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
        n <- floor(0.5*n)
        green_vec <- green_int + green_slope*n
        red_vec <- red_int + red_slope*n
        if (sort(red_vec)[mid1] > sort(green_vec)[mid2]){
          lower <- n
        }
      }
    }
  }
  
  ## implement binary search given these bounds for the sample size
  while ((upper - lower) > 1){
    n <- ceiling(0.5*(upper + lower))
    
    ## use the linear approximations to select the points from the Sobol' 
    ## sequence in a targeted way; these points will be used to explore the 
    ## sampling distributions of posterior probabilities in a targeted way
    green_vec <- green_int + green_slope*n
    red_vec <- red_int + red_slope*n
    ## get the indices for these points
    sub_green <- which(rank(green_vec) >= low2 & rank(green_vec) <= high2)
    sub_red <- which(rank(red_vec) >= low1 & rank(red_vec) <= high1)
    
    green_vec <- NULL
    red_vec <- NULL
    for (i in 1:length(sub_green)){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[sub_green[i],]),
                                  delta = deltas, u = sob_green[sub_green[i],2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[sub_red[i],]),
                                delta = deltas, u = sob_red[sub_red[i],2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## for the unselected points, we just assume that their probabilities are less (greater)
    ## than the relevant order statistic if their anticipated probability based on the 
    ## linear approximations on the logit scale is less (greater) than the relevant order statistic
    red_aug <- c(rep(c(min(red_vec)-1), low1 - 1), red_vec, rep(c(max(red_vec)+1), m - high1))
    green_aug <- c(rep(c(min(green_vec)-1), low2 - 1), green_vec, rep(c(max(green_vec)+1), m - high2))
    if (sort(red_aug)[mid1] <= sort(green_aug)[mid2]){
      green_max <- green_vec; red_max <- red_vec
      sub_green_max <- sub_green; sub_red_max <- sub_red
      upper <- n
    }
    else{
      lower <- n
    }
  }
  
  n <- upper
  ## we don't need to reapproximate the same posterior probabilities if we have already explored
  ## this sample size
  if (n == n_high){
    
    ## do not need to do this; it is just so we have three sample sizes to construct each contour plot
    n <- ifelse(n0 == n_high, max(ceiling(1.1*n0), n0 + 5), max(ceiling(1.1*n1), n1 + 5))
    green_vec <- NULL
    red_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    green_vec_n2 <- green_vec
    red_vec_n2 <- red_vec
    n2 <- n
    
    n_final <- ifelse(n1 > n0, 1, 0)
    
    ## first component of list is (n, gamma, relevant order statistic of green probabilities, confirmatory type I error
    ## estimate, confirmatory power estimate, n^((0)), n^((1)), n^((2)))
    ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities from
    ## the red and green regions for the sample sizes n^((0)), n^((1)), and n^((2)))
    ## we sort the sample sizes and the posterior probabilities to help with making the plot
    results <- list(c(get(paste0("n", n_final)), 1/(1+exp(-as.numeric(sort(get(paste0("red_vec_n", n_final)))[mid1]))), 
                      1/(1+exp(-as.numeric(sort(get(paste0("green_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("red_vec_n", n_final)) > sort(get(paste0("red_vec_n", n_final)))[mid1])), 
                      as.numeric(mean(get(paste0("green_vec_n", n_final)) > sort(get(paste0("red_vec_n", n_final)))[mid1])),
                      sort(c(n0, n1, n2))),
                    get(paste0("red_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("green_vec_n", which.min(c(n0, n1, n2))-1)), 
                    get(paste0("red_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("green_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                    get(paste0("red_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("green_vec_n", which.max(c(n0, n1, n2))-1)))
    
    if (contour == TRUE){
      return(results)
    }
    else{
      return(results[[1]][1:5])
    }
  }
  
  if (n == n_low){
    
    ## do not need to do this; it is just so we have three sample sizes to construct each contour plot
    n <- ifelse(n0 == n_low, min(floor(0.9*n0), n0 - 5), min(floor(0.9*n1), n1 - 5))
    green_vec <- NULL
    red_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    green_vec_n2 <- green_vec
    red_vec_n2 <- red_vec
    n2 <- n
    
    n_final <- ifelse(n1 < n0, 1, 0)
    
    ## first component of list is (n, gamma, relevant order statistic of green probabilities, confirmatory type I error
    ## estimate, confirmatory power estimate, n^((0)), n^((1)), n^((2)))
    ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities from
    ## the red and green regions for the sample sizes n^((0)), n^((1)), and n^((2)))
    ## we sort the sample sizes and the posterior probabilities to help with making the plot
    results <- list(c(get(paste0("n", n_final)), 1/(1+exp(-as.numeric(sort(get(paste0("red_vec_n", n_final)))[mid1]))), 
                      1/(1+exp(-as.numeric(sort(get(paste0("green_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("red_vec_n", n_final)) > sort(get(paste0("red_vec_n", n_final)))[mid1])), 
                      as.numeric(mean(get(paste0("green_vec_n", n_final)) > sort(get(paste0("red_vec_n", n_final)))[mid1])),
                      sort(c(n0, n1, n2))),
                    get(paste0("red_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("green_vec_n", which.min(c(n0, n1, n2))-1)), 
                    get(paste0("red_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("green_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                    get(paste0("red_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("green_vec_n", which.max(c(n0, n1, n2))-1)))
    
    if (contour == TRUE){
      return(results)
    }
    else{
      return(results[[1]][1:5])
    }
  }
  
  ## otherwise, we approximate the remaining posterior probabilities that were not
  ## selected at the final sample size recommendation
  green_vec <- NULL
  red_vec <- NULL
  remain_green <- subset(1:m, ! 1:m %in% sub_green_max)
  remain_red <- subset(1:m, ! 1:m %in% sub_red_max)
  for (i in remain_green){
    green_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(green_par[i,]),
                                delta = deltas, u = sob_green[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
  }
  for (i in remain_red){
    red_vec[i] <- targetPower(n_val = n, q = q,
                              params = as.numeric(red_par[i,]),
                              delta = deltas, u = sob_red[i,2:9],
                              hyper = cbind(alphas1, betas1, alphas2, betas2))
  }
  green_vec[sub_green_max] <- green_max; red_vec[sub_red_max] <- red_max
  
  ## save the logits of the posterior probabilities
  green_vec_n2 <- green_vec
  red_vec_n2 <- red_vec
  n2 <- n
  
  ## first component of list is (n, gamma, relevant order statistic of green probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate, n^((0)), n^((1)), n^((2)))
  ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities from
  ## the red and green regions for the sample sizes n^((0)), n^((1)), and n^((2)))
  ## we sort the sample sizes and the posterior probabilities to help with making the plot
  
  results <- list(c(n, 1/(1 + exp(-as.numeric(sort(red_vec)[mid1]))),
                    1/(1 + exp(-as.numeric(sort(green_vec)[mid2]))),
                    as.numeric(mean(red_vec > sort(red_vec)[mid1])), as.numeric(mean(green_vec > sort(red_vec)[mid1])),
                    sort(c(n0, n1, n2))),
                  get(paste0("red_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("green_vec_n", which.min(c(n0, n1, n2))-1)), 
                  get(paste0("red_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("green_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                  get(paste0("red_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("green_vec_n", which.max(c(n0, n1, n2))-1)))
  
  if (contour == TRUE){
    return(results)
  }
  else{
    return(results[[1]][1:5])
  }
}

## alternative function to explore sample sizes using all points from the 
## hypercube for each sample size that we consider
findGammaFull <- function(green_par, red_par, pwr, typeI, deltas, q, 
                          alphas1, betas1, alphas2, betas2, m = 8192, m0 = 512,
                          seed_green = 1, seed_red = 2, upper_init = 1000){
  
  ## generate Sobol' sequence for each sampling distribution
  sob_green <- sobol(m, d = 9, randomize = "digital.shift", seed = seed_green)
  sob_red <- sobol(m, d = 9, randomize = "digital.shift", seed = seed_red)
  
  ## index the eta values from the design priors such that the anticipated value for
  ## theta increases with the index; we do this using the first dimension of the Sobol' sequence
  green_means1 <- 5 - 4*green_par[,1] - 3*green_par[,2] - 2*green_par[,3] - green_par[,4]
  green_means2 <- 5 - 4*green_par[,5] - 3*green_par[,6] - 2*green_par[,7] - green_par[,8]
  green_diff <- green_means1 - green_means2
  green_par <- green_par[order(green_diff),]
  green_par <- green_par[ceiling(length(green_diff)*sob_green[,1]),]
  
  ## repeat this process for the red region
  red_means1 <- 5 - 4*red_par[,1] - 3*red_par[,2] - 2*red_par[,3] - red_par[,4]
  red_means2 <- 5 - 4*red_par[,5] - 3*red_par[,6] - 2*red_par[,7] - red_par[,8]
  red_diff <- red_means1 - red_means2
  red_par <- red_par[order(red_diff),]
  red_par <- red_par[ceiling(length(red_diff)*sob_red[,1]),]
  
  if (m < m0){stop("m0 should be less than m")}
  lower <- 1
  n <- upper_init
  upper <- upper_init
  green_vec <- rep(0, m)
  red_vec <- rep(1, m)
  
  mid1 <- ceiling((1 - typeI)*m)
  mid2 <- floor((1-pwr)*m)
  
  ## find a larger upper bound if the initial upper bound provided
  ## is not large enough
  while (sort(red_vec)[mid1] > sort(green_vec)[mid2]){
    green_vec <- NULL
    red_vec <- NULL
    ## implement Algorithm 3 for each of the first m0 points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    ## increase upper bound if the criterion on the order statistics is not satisfied
    if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
      green_max <- green_vec; red_max <- red_vec
      upper <- n
      n_max <- n
    }
    else{
      lower <- n
      green_min <- green_vec; red_min <- red_vec
      n_min <- n
      n <- 2*n
    }
  }
  
  ## implement binary search given these bounds for the sample size
  while ((upper - lower) > 1){
    n <- ceiling(0.5*(upper + lower))
    
    green_vec <- NULL
    red_vec <- NULL
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:9],
                                  hyper = cbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:9],
                                hyper = cbind(alphas1, betas1, alphas2, betas2))
    }
    
    if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
      green_max <- green_vec; red_max <- red_vec
      upper <- n
    }
    else{
      lower <- n
    }
  }
  
  n <- upper
  
  ## list is (n, gamma, relevant order statistic of green probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate)
  return(list(c(n, 1/(1 + exp(-as.numeric(sort(red_vec)[mid1]))),
                1/(1 + exp(-as.numeric(sort(green_vec)[mid2]))),
                as.numeric(mean(red_vec > sort(red_vec)[mid1])), as.numeric(mean(green_vec > sort(red_vec)[mid1])))))
  
}