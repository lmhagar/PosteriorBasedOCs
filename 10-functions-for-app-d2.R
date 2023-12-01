## functions for Weibull numerical study, where Algorithm 4 is used
## instead of Algorithm 1

## BEGIN SETUP ##

## load necessary packages
require(qrng)
require(nleqslv)

## function to solve for the posterior mode of each group
## x and y are parameterized as for the previous function
fn_grad <- function(x, y, mu, tau, kappa, lambda) {
  
  res1 <- -y[3]*(x[1] - y[1])*exp(y[2])^2 + y[3]*(x[2] - y[2])*exp(y[2])*(1 + digamma(1)) - exp(x[1])*tau + mu
  
  res2 <- y[3]*(x[1] - y[1])*exp(y[2])*(1 + digamma(1)) - y[3]*(x[2] - y[2])*y[4] - exp(x[2])*lambda + kappa
  
  return(c(res1, res2))
}

calc_covar <- function(u, v, yy_star, hyper){
  l <- exp(u); k <- exp(v); n <- yy_star[3]
  tau <- hyper[1]; lambda <- hyper[2]
  d11 <- n*k^2 + l*tau
  d12 <- -n*k*(1 + digamma(1))
  d22 <- n*yy_star[4] + k*lambda
  mat <- rbind(c(d11, d12), c(d12, d22))
  return(solve(mat))
}

## this helper function implements Algorithm 4 with the Weibull example to
## return the logit of the posterior probability (used with the function FindGamma())
targetPower <- function(u, params, deltas, n_val, hyper, q, threshold){
  
  ## return negative power if sample size is not positive
  if (n_val <= 0){return(-1.5)}
  
  wei_lambda.1 <- params[1]
  wei_k.1 <- params[2]
  wei_lambda.2 <- params[3]
  wei_k.2 <- params[4]
  
  emc <- pi^2/6 + (digamma(1))^2 + 2*digamma(1)
  ## generate approximately normal MLEs for group 1 using delta method
  rho1 <- (1 + digamma(1))/sqrt(1 + emc)
  mat1 <- matrix((sqrt(6)/(pi*wei_k.1))^2*c(1 + emc, wei_k.1*(1 + digamma(1)), wei_k.1*(1 + digamma(1)), wei_k.1^2), nrow = 2)
  l1 <- qnorm(u[1], log(wei_lambda.1), sqrt(mat1[1,1]/n_val))
  k1 <- qnorm(u[2], log(wei_k.1) + rho1*(l1 - log(wei_lambda.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/n_val))
  
  ## generate approximately normal MLEs for group 2 using delta method
  rho2 <- (1 + digamma(1))/sqrt(1 + emc)
  mat2 <- matrix((sqrt(6)/(pi*wei_k.2))^2*c(1 + emc, wei_k.2*(1 + digamma(1)), wei_k.2*(1 + digamma(1)), wei_k.2^2), nrow = 2)
  l2 <- qnorm(u[3], log(wei_lambda.2), sqrt(mat2[1,1]/round(q*n_val)))
  k2 <- qnorm(u[4], log(wei_k.2) + rho2*(l2 - log(wei_lambda.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/round(q*n_val)))
  
  ## summarize information from first group of data (faster computation)
  yy_star1 <- c(l1, k1, n_val, 1 + emc)
  ## find posterior modes for the first group (logalpha and logbeta)
  modes <- nleqslv(c(l1, k1), fn_grad, y = yy_star1, mu = hyper[1,1], tau = hyper[1,2],
                   kappa = hyper[2,1], lambda = hyper[2,2] )$x
  
  mat1_new <- calc_covar(modes[1], modes[2], yy_star1, c(hyper[1,2], hyper[2,2]))
  ## exponentiate modes to return to standard scale
  modes1 <- exp(modes)
  
  ## repeat all steps for the second group
  yy_star2 <- c(l2, k2, round(q*n_val), 1 + emc)
  modes <- nleqslv(c(l2, k2), fn_grad, y = yy_star2, mu = hyper[3,1], tau = hyper[3,2],
                   kappa = hyper[4,1], lambda = hyper[4,2] )$x
  
  mat2_new <- calc_covar(modes[1], modes[2], yy_star2, c(hyper[3,2], hyper[4,2]))
  modes2 <- exp(modes)
  l1 <- modes1[1]; k1 <- modes1[2]; l2 <- modes2[1]; k2 <- modes2[2]
  
  ## ensure no modes are 0 due to underflow errors
  if(max(l1 <= 0, k1 <= 0, l2 <= 0, k2<= 0)){return(-1.5)}
  
  ## define theta metrics in terms of design values
  theta1 <- l1*log(1/(1 - threshold))^(1/k1)
  theta2 <- l2*log(1/(1 - threshold))^(1/k2)
  
  if(max(theta1 <= 0, !is.finite(theta1), theta2 <= 0, !is.finite(theta2))){return(-1.5)}
  
  ## compute partial derivatives of theta with respect to logalpha and logbeta for each group
  ## this is different from the previous function that computes the derivatives with respect
  ## to alpha and beta
  d_lambda1 <- l1*log(1/(1 - threshold))^(1/k1)
  
  d_k1 <- -(l1/k1)*(-log(1 - threshold))^(1/k1)*log(-log(1 - threshold))
  
  d_lambda2 <- l2*log(1/(1 - threshold))^(1/k2)
  
  d_k2 <- -(l2/k2)*(-log(1 - threshold))^(1/k2)*log(-log(1 - threshold))
  
  ## apply the delta method to get the limiting variance for each theta_j metric
  avar1 <- t(c(d_lambda1, d_k1))%*%mat1_new%*%c(d_lambda1, d_k1)
  
  avar2 <- t(c(d_lambda2, d_k2))%*%mat2_new%*%c(d_lambda2, d_k2)
  
  ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
  Fish_ratio_mu <- avar1/theta1^2 + avar2/theta2^2 
  
  ## return negative power if division causes NA/Inf values
  if (is.na(Fish_ratio_mu)){return(-1.5)}
  if (Fish_ratio_mu < -10000){return(-1.5)}
  
  ## return power based on normal approximation induced by Bernstein-von Mises
  realP <- pnorm(deltas[2], log(theta1/theta2), sqrt(Fish_ratio_mu)) - pnorm(deltas[1], log(theta1/theta2), sqrt(Fish_ratio_mu))
  
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

## Code to implement Algorithm 2 for the Weibull example. We explore sample sizes
## and return the optimal design (i.e., the (n, gamma) combination)
findGamma <- function(green_par, red_par, pwr, typeI, deltas, q,
                      alphas1, betas1, alphas2, betas2, m = 8192, m0 = 512, threshold = 0.9,
                      seed_green = 1, seed_red = 2, upper_init = 1000, prng = FALSE, 
                      contour = FALSE){
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option available for comparison
  if (prng == FALSE){
    sob_green <- sobol(m, d = 5, randomize = "digital.shift", seed = seed_green)
    sob_red <- sobol(m, d = 5, randomize = "digital.shift", seed = seed_red)
  }
  else{
    set.seed(seed_green); sob_green <- matrix(runif(5*m), ncol = 5)
    
    set.seed(seed_red); sob_red <- matrix(runif(5*m), ncol = 5)
  }
  
  ## index the eta values from the design priors such that the anticipated value for
  ## theta increases with the index; we do this using the first dimension of the Sobol' sequence
  ## we have a difference here on the logarithmic scale
  green_means1 <- log(green_par[,1]*log(1/(1 - threshold))^(1/green_par[,2]))
  green_means2 <- log(green_par[,3]*log(1/(1 - threshold))^(1/green_par[,4]))
  green_diff <- green_means1 - green_means2
  green_par <- green_par[order(green_diff),]
  green_par <- green_par[ceiling(length(green_diff)*sob_green[,1]),]
  
  ## repeat this process for the red region
  red_means1 <- log(red_par[,1]*log(1/(1 - threshold))^(1/red_par[,2]))
  red_means2 <- log(red_par[,3]*log(1/(1 - threshold))^(1/red_par[,4]))
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
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
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
    print(c(lower, upper, sort(red_vec)[mid1], sort(green_vec)[mid2]))
  }
  
  ## get a lower point to construct linear approximations to the logits
  ## of the posterior probabilities (this is only required if the initial
  ## upper bound for n was sufficiently large)
  if(lower == 1){
    n <- ceiling(0.5*(upper + lower))
    green_vec <- NULL
    red_vec <- NULL
    for (i in 1:m0){
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
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
    print(c(lower, upper, sort(red_vec)[mid1], sort(green_vec)[mid2]))
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
    print(c(lower_temp, upper_temp, sort(red_vec)[mid1], sort(green_vec)[mid2]))
  }
  
  print(c(lower_temp, upper_temp))
  
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
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
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
    print(c(lower, upper, sort(red_vec)[mid1], sort(green_vec)[mid2]))
    
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
      print(c(lower_temp, upper_temp, sort(red_vec)[mid1], sort(green_vec)[mid2]))
    }
    n <- upper_temp
    ## we check how much the sample size recommendation changed with respect to
    ## what we started this iteration of the floor loop with; if there is not
    ## much change, then we just take this new recommendation as n(0)
    print(c(n_last, n))
  }
  
  print(c(n, lower, upper, lower_temp, upper_temp, n_last == n_min, n_last == n_max))
  
  ## compute the posterior probabilities using Algorithm 3 for all m points for
  ## each Sobol' sequence at the initial sample size n(0)
  n_init <- n
  green_vec <- NULL
  red_vec <- NULL
  for (i in 1:m){
    green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(green_par[i,]),
                                delta = deltas, u = sob_green[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                              params = as.numeric(red_par[i,]),
                              delta = deltas, u = sob_red[i,2:5],
                              hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  
  print(c(n, quantile(red_vec, 1 - typeI), quantile(green_vec, 1 - pwr)))
  
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
    print(n)
    green_vec <- NULL
    red_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
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
        print(c(lower, upper, sort(red_vec)[mid1] , sort(green_vec)[mid2]))
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
    print(n)
    green_vec <- NULL
    red_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
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
        print(c(lower, upper, sort(red_vec)[mid1], sort(green_vec)[mid2]))
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
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[sub_green[i],]),
                                  delta = deltas, u = sob_green[sub_green[i],2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[sub_red[i],]),
                                delta = deltas, u = sob_red[sub_red[i],2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
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
    print(c(lower, upper, sort(red_aug)[mid1], sort(green_aug)[mid2]))
  }
  
  n <- upper
  print(n)
  ## we don't need to reapproximate the same posterior probabilities if we have already explored
  ## this sample size
  if (n == n_high){
    
    ## do not need to do this; it is just so we have three sample sizes to construct each contour plot
    n <- ifelse(n0 == n_high, max(ceiling(1.1*n0), n0 + 5), max(ceiling(1.1*n1), n1 + 5))
    print(n)
    green_vec <- NULL
    red_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
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
    print(n)
    green_vec <- NULL
    red_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    green_vec_n2 <- green_vec
    red_vec_n2 <- red_vec
    n2 <- n
    
    n_final <- ifelse(n1 < n0, 1, 0)
    
    print(c("End:", n0, n1, n2))
    
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
    green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(green_par[i,]),
                                delta = deltas, u = sob_green[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  for (i in remain_red){
    red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                              params = as.numeric(red_par[i,]),
                              delta = deltas, u = sob_red[i,2:5],
                              hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  green_vec[sub_green_max] <- green_max; red_vec[sub_red_max] <- red_max
  
  ## save the logits of the posterior probabilities
  green_vec_n2 <- green_vec
  red_vec_n2 <- red_vec
  n2 <- n
  
  print(c(n, sort(red_vec)[mid1], sort(green_vec)[mid2]))
  
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
## hypercube for each sample size that we consider (slightly modified for use
## with the Weibull example)
findGammaFull <- function(green_par, red_par, pwr, typeI, deltas, q,
                          alphas1, betas1, alphas2, betas2, m = 8192, m0 = 512, threshold = 0.9,
                          seed_green = 1, seed_red = 2, upper_init = 1000){
  
  ## generate Sobol' sequence for each sampling distribution
  sob_green <- sobol(m, d = 5, randomize = "digital.shift", seed = seed_green)
  sob_red <- sobol(m, d = 5, randomize = "digital.shift", seed = seed_red)
  
  ## index the eta values from the design priors such that the anticipated value for
  ## theta increases with the index; we do this using the first dimension of the Sobol' sequence
  green_means1 <- log(green_par[,1]*log(1/(1 - threshold))^(1/green_par[,2]))
  green_means2 <- log(green_par[,3]*log(1/(1 - threshold))^(1/green_par[,4]))
  green_diff <- green_means1 - green_means2
  green_par <- green_par[order(green_diff),]
  green_par <- green_par[ceiling(length(green_diff)*sob_green[,1]),]
  
  ## repeat this process for the red region
  red_means1 <- log(red_par[,1]*log(1/(1 - threshold))^(1/red_par[,2]))
  red_means2 <- log(red_par[,3]*log(1/(1 - threshold))^(1/red_par[,4]))
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
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
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
    print(c(lower, upper, sort(red_vec)[mid1], sort(green_vec)[mid2]))
  }
  
  ## implement binary search given these bounds for the sample size
  while ((upper - lower) > 1){
    n <- ceiling(0.5*(upper + lower))
    
    green_vec <- NULL
    red_vec <- NULL
    for (i in 1:m){
      green_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                  params = as.numeric(green_par[i,]),
                                  delta = deltas, u = sob_green[i,2:5],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      red_vec[i] <- targetPower(n_val = n, q = q, threshold = threshold,
                                params = as.numeric(red_par[i,]),
                                delta = deltas, u = sob_red[i,2:5],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    if (sort(red_vec)[mid1] <= sort(green_vec)[mid2]){
      green_max <- green_vec; red_max <- red_vec
      upper <- n
    }
    else{
      lower <- n
    }
    print(c(lower, upper, sort(red_vec)[mid1], sort(green_vec)[mid2]))
  }
  
  n <- upper
  print(n)
  
  ## list is (n, gamma, relevant order statistic of green probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate)
  return(list(c(n, 1/(1 + exp(-as.numeric(sort(red_vec)[mid1]))),
                1/(1 + exp(-as.numeric(sort(green_vec)[mid2]))),
                as.numeric(mean(red_vec > sort(red_vec)[mid1])), as.numeric(mean(green_vec > sort(red_vec)[mid1])))))
  
}