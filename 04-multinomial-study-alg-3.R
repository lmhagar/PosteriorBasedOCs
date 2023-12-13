## code to produce figures in multinomial numerical study, where Algorithm 3 is used
## instead of Algorithm 1 (Figure 3 and the figure in Appendix C.2)

## run this file AFTER running file 03-functions-for-sec-4-5.R

## BEGIN SETUP ##

## load necessary packages
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(ggplot2)
require(ggpubr)
require(cowplot)

## we first generate the eta values from the design priors
## extract design prior parameters for groups 1 and 2
## from the Shiny app
alphas1 <- read.csv("priors_group1.csv")$alpha
betas1 <- read.csv("priors_group1.csv")$beta

alphas2 <- read.csv("priors_group2.csv")$alpha
betas2 <- read.csv("priors_group2.csv")$beta

## generate observations from the design priors
z11 = rbeta(5000000, alphas1[1], betas1[1]); z12 = rbeta(5000000, alphas1[2], betas1[2])
z13 = rbeta(5000000, alphas1[3], betas1[3]); z14 = rbeta(5000000, alphas1[4], betas1[4])
z21 = rbeta(5000000, alphas2[1], betas2[1]); z22 = rbeta(5000000, alphas2[2], betas2[2])
z23 = rbeta(5000000, alphas2[3], betas2[3]); z24 = rbeta(5000000, alphas2[4], betas2[4])

## convert from the Z-scale to the p-scale for both groups
p11 <- z11
p12 <- z12*(1-z11)
p13 <- z13*(1-z12)*(1-z11)
p14 <- z14*(1-z13)*(1-z12)*(1-z11)

p21 <- z21
p22 <- z22*(1-z21)
p23 <- z23*(1-z22)*(1-z21)
p24 <- z24*(1-z23)*(1-z22)*(1-z21)


## tweak any multinomial probabilities that round to 0 or 1
## (to four decimal places); the "problem_inds" are the indices
## that we need to take care of
p11 <- pmax(z11, 0.001)
p12 <- pmax(z12*(1-z11), 0.001)
p13 <- pmax(z13*(1-z12)*(1-z11), 0.001)
p14 <- pmax(z14*(1-z13)*(1-z12)*(1-z11), 0.001)

## slightly rescale if necessary for group 1
problem_inds <- which(p11 + p12 + p13 + p14 > 1 - 0.001)
sums1 <- p11 + p12 + p13 + p14

if (length(problem_inds) > 0){
  for (i in 1:length(problem_inds)){
    p11[i] <- (1 - 0.001)*p11[i]/sums1[problem_inds[i]]
    p12[i] <- (1 - 0.001)*p12[i]/sums1[problem_inds[i]]
    p13[i] <- (1 - 0.001)*p13[i]/sums1[problem_inds[i]]
    p14[i] <- (1 - 0.001)*p14[i]/sums1[problem_inds[i]]
  }
}

## repeat this process for the second group
p21 <- pmax(z21, 0.001)
p22 <- pmax(z22*(1-z21), 0.001)
p23 <- pmax(z23*(1-z22)*(1-z21), 0.001)
p24 <- pmax(z24*(1-z23)*(1-z22)*(1-z21), 0.001)

## slightly rescale if necessary for group 2
problem_inds <- which(p21 + p22 + p23 + p24 > 1 - 0.001)
sums2 <- p21 + p22 + p23 + p24 

if (length(problem_inds) > 0){
  for (i in 1:length(problem_inds)){
    p21[problem_inds[i]] <- (1 - 0.001)*p21[problem_inds[i]]/sums2[problem_inds[i]]
    p22[problem_inds[i]] <- (1 - 0.001)*p22[problem_inds[i]]/sums2[problem_inds[i]]
    p23[problem_inds[i]] <- (1 - 0.001)*p23[problem_inds[i]]/sums2[problem_inds[i]]
    p24[problem_inds[i]] <- (1 - 0.001)*p24[problem_inds[i]]/sums2[problem_inds[i]]
  }
}

## obtain theta values for each draw to segment the priors
theta1 <- p11 + 2*p12 + 3*p13 + 4*p14 + 5*(1 - p11 - p12 - p13 - p14)
theta2 <- p21 + 2*p22 + 3*p23 + 4*p24 + 5*(1 - p21 - p22 - p23 - p24)

## define the green and red regions
delta_L <- -0.5; green_mid <- -0.2; green_hw <- 0.1
green_low <- green_mid - green_hw; green_high <- green_mid + green_hw
red_high <- delta_L; red_low <- red_high - 0.05

## approximate the induced prior density function on theta for sampling-resampling
diff_dens <- density(theta1 - theta2)
diff_pdf <- approxfun(x = diff_dens$x,
                      y = diff_dens$y)

## conduct sampling-resampling to uniformly sample from the green region
diff_num <- ifelse(theta1 - theta2 <= green_high, theta1 - theta2 >= green_low, 0)
diff_denom <- diff_pdf(theta1 - theta2)

weights <- diff_num/diff_denom
weights <- weights/sum(weights)

## resample to get 8192 parameter values from green region
inds_region <- sample(1:5000000, size = 8192, replace = TRUE, prob = weights)

## output parameter values to a .csv file
write.csv(cbind(z11, z12, z13, z14)[inds_region,], "zs1_green.csv", row.names = FALSE)
write.csv(cbind(p11, p12, p13, p14)[inds_region,], "ps1_green.csv", row.names = FALSE)
write.csv(theta1[inds_region], "theta1_green.csv", row.names = FALSE)

write.csv(cbind(z21, z22, z23, z24)[inds_region,], "zs2_green.csv", row.names = FALSE)
write.csv(cbind(p21, p22, p23, p24)[inds_region,], "ps2_green.csv", row.names = FALSE)
write.csv(theta2[inds_region], "theta2_green.csv", row.names = FALSE)

## conduct sampling-resampling to uniformly sample from the red region
diff_num <- ifelse(theta1 - theta2 <= red_high, theta1 - theta2 >= red_low, 0)
diff_denom <- diff_pdf(theta1 - theta2)

weights <- diff_num/diff_denom
weights <- weights/sum(weights)

## resample to get 8192 parameter values from red region
inds_region2 <- sample(1:5000000, size = 8192, replace = TRUE, prob = weights)

## output parameter values to a .csv file
write.csv(cbind(z11, z12, z13, z14)[inds_region2,], "zs1_red.csv", row.names = FALSE)
write.csv(cbind(p11, p12, p13, p14)[inds_region2,], "ps1_red.csv", row.names = FALSE)
write.csv(theta1[inds_region2], "theta1_red.csv", row.names = FALSE)

write.csv(cbind(z21, z22, z23, z24)[inds_region2,], "zs2_red.csv", row.names = FALSE)
write.csv(cbind(p21, p22, p23, p24)[inds_region2,], "ps2_red.csv", row.names = FALSE)
write.csv(theta2[inds_region2], "theta2_red.csv", row.names = FALSE)

## save these parameter values to one matrix for each region
greens <- read.csv("ps1_green.csv")
greens <- cbind(greens, read.csv("ps2_green.csv"))
reds <- read.csv("ps1_red.csv")
reds <- cbind(reds, read.csv("ps2_red.csv"))
greens <- matrix(unlist(greens), ncol = 8)
reds <- matrix(unlist(reds), ncol = 8)

## define hyperparameters for analysis priors (given in Section 2)
alphas1 <- 0.8*rep(1,4)
betas1 <- 0.8*c(4,3,2,1)
alphas2 <- 0.8*rep(1,4)
betas2 <- 0.8*c(4,3,2,1)

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## repeat the sample size calculation for the illustrative example 1000 times
## with the different Sobol' sequences
## delta_U is 4 instead of Inf since 4 is the maximum value for theta
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGamma(green_par = greens, red_par = reds, 
                                        pwr = 0.8, typeI = 0.05, deltas = c(-0.5, 4), q = 1.25, 
                                        alphas1, betas1, alphas2, betas2, m = 8192,
                                        seed_green = i, seed_red = i + 1000, contour = TRUE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp[,1:9], "sec4p3_sobol_8192_summary.csv", row.names = FALSE)

## output the red and green posterior probabilities for each of the three sample
## sizes n^((0)) = small, n^((1)) = mid, and n^((2)) = large
write.csv(temp[,10:8201], "sec4p3_sobol_8192_red_small.csv", row.names = FALSE)
write.csv(temp[,8202:16393], "sec4p3_sobol_8192_green_small.csv", row.names = FALSE)
write.csv(temp[,16394:24585], "sec4p3_sobol_8192_red_mid.csv", row.names = FALSE)
write.csv(temp[,24586:32777], "sec4p3_sobol_8192_green_mid.csv", row.names = FALSE)
write.csv(temp[,32778:40969], "sec4p3_sobol_8192_red_large.csv", row.names = FALSE)
write.csv(temp[,40970:49161], "sec4p3_sobol_8192_green_large.csv", row.names = FALSE)

## repeat the 1000 sample size calculations with the SAME Sobol' sequences but
## explore the entire hypercube and each sample size considered
tempFull <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGammaFull(green_par = greens, red_par = reds, 
                                            pwr = 0.8, typeI = 0.05, deltas = c(-0.5, 4), q = 1.25, 
                                            alphas1, betas1, alphas2, betas2, m = 8192,
                                            seed_green = i, seed_red = i + 1000)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

write.csv(tempFull, "sec4p3_sobol_8192_full_summary.csv", row.names = FALSE)

## repeat the sample size calculation for the illustrative example 1000 times
## using pseudorandom sequences of length 8192
## delta_U is 4 instead of Inf since 4 is the maximum value for theta
tempPRNG <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGamma(green_par = greens, red_par = reds, 
                                        pwr = 0.8, typeI = 0.05, deltas = c(-0.5, 4), q = 1.25, 
                                        alphas1, betas1, alphas2, betas2, m = 8192, prng = TRUE,
                                        seed_green = i + 2000, seed_red = i + 3000, contour = FALSE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(tempPRNG, "sec4p3_prng_8192_summary.csv", row.names = FALSE)

## repeat the sample size calculation for the illustrative example 1000 times
## using pseudorandom sequences of length 24000
## delta_U is 4 instead of Inf since 4 is the maximum value for theta
tempPRNG24k <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                    .options.snow=opts, .errorhandling = "remove") %dopar% {
                      temp_res <- findGamma(green_par = greens, red_par = reds, 
                                            pwr = 0.8, typeI = 0.05, deltas = c(-0.5, 4), q = 1.25, 
                                            alphas1, betas1, alphas2, betas2, m = 24000, prng = TRUE,
                                            seed_green = i + 4000, seed_red = i + 5000, contour = FALSE)
                      
                      c(i, as.numeric(unlist(temp_res)))
                    }

## output the optimal design and summary for each simulation repetition
write.csv(tempPRNG24k, "sec4p3_prng_24k_summary.csv", row.names = FALSE)

## approximate the sampling distributions of posterior probabilities by
## generating data directly (for the contour plots)
resamps <- 81920
pb <- txtProgressBar(max = resamps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get parameter values from design priors to generate
## the multinomial data for the green region
samps <- seq(100, 140, 1)
paramss_in <- read.csv("ps1_green.csv")
paramss_in <- cbind(paramss_in, read.csv("ps2_green.csv"))
paramss1 <- cbind(paramss_in[1:4], 1 - rowSums(paramss_in[1:4]))
paramss2 <- cbind(paramss_in[5:8], 1 - rowSums(paramss_in[5:8]))

paramss1 <- rbind(paramss1, paramss1, paramss1, paramss1, paramss1, 
                  paramss1, paramss1, paramss1, paramss1, paramss1)

paramss2 <- rbind(paramss2, paramss2, paramss2, paramss2, paramss2, 
                  paramss2, paramss2, paramss2, paramss2, paramss2)
delta_L <- -0.5
q <- 1.25

for (k in 1:length(samps)){
  probs_in <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        ## generate multinomial data and generate from the beta posteriors
                        dat1 <- rmultinom(1, samps[k], prob = as.numeric(paramss1[j,]))
                        z1 <- rbeta(100000, dat1[1] + alphas1[1], sum(dat1[2:5]) + betas1[1])
                        z2 <- rbeta(100000, dat1[2] + alphas1[2], sum(dat1[3:5]) + betas1[2])
                        z3 <- rbeta(100000, dat1[3] + alphas1[3], sum(dat1[4:5]) + betas1[3])
                        z4 <- rbeta(100000, dat1[4] + alphas1[4], dat1[5]  + betas1[4])
                        
                        ## convert from the Z-scale to the p-scale
                        p2 <- z2*(1-z1); p3 <- z3*(1-z1)*(1-z2); p4 <- z4*(1-z3)*(1-z1)*(1-z2)
                        p5 <- 1 - z1 - p2 - p3 - p4
                        
                        ## get posterior draws for theta1
                        theta.1 <- z1 + 2*p2 + 3*p3 + 4*p4 + 5*p5
                        
                        ## repeat for model 2
                        dat2 <- rmultinom(1, round(q*samps[k]), prob = as.numeric(paramss2[j,]))
                        z1 <- rbeta(100000, dat2[1] + alphas2[1], sum(dat2[2:5]) + betas2[1])
                        z2 <- rbeta(100000, dat2[2] + alphas2[2], sum(dat2[3:5]) + betas2[2])
                        z3 <- rbeta(100000, dat2[3] + alphas2[3], sum(dat2[4:5]) + betas2[3])
                        z4 <- rbeta(100000, dat2[4] + alphas2[4], dat2[5]  + betas2[4])
                        
                        p2 <- z2*(1-z1); p3 <- z3*(1-z1)*(1-z2); p4 <- z4*(1-z3)*(1-z1)*(1-z2)
                        p5 <- 1 - z1 - p2 - p3 - p4
                        
                        theta.2 <- z1 + 2*p2 + 3*p3 + 4*p4 + 5*p5
                        
                        ## output posterior probability, where delta_U is Inf
                        mean(ifelse(theta.1 - theta.2 >= delta_L, 1, 0))
                      }
  
  ## output results to .csv file
  write.csv(probs_in, paste0("probs_green_sec5_", samps[k], ".csv"), row.names = FALSE)
}

## get parameter values from design priors to generate
## the multinomial data for the red region
paramss_out <- read.csv("ps1_red.csv")
paramss_out <- cbind(paramss_out, read.csv("ps2_red.csv"))
paramss1 <- cbind(paramss_out[1:4], 1 - rowSums(paramss_out[1:4]))
paramss2 <- cbind(paramss_out[5:8], 1 - rowSums(paramss_out[5:8]))

paramss1 <- rbind(paramss1, paramss1, paramss1, paramss1, paramss1, 
                  paramss1, paramss1, paramss1, paramss1, paramss1)

paramss2 <- rbind(paramss2, paramss2, paramss2, paramss2, paramss2, 
                  paramss2, paramss2, paramss2, paramss2, paramss2)

for (k in 1:length(samps)){

  probs_out <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts, .errorhandling = "remove") %dopar% {
                         
                         ## generate multinomial data and generate from the beta posteriors
                         dat1 <- rmultinom(1, samps[k], prob = as.numeric(paramss1[j,]))
                         z1 <- rbeta(100000, dat1[1] + alphas1[1], sum(dat1[2:5]) + betas1[1])
                         z2 <- rbeta(100000, dat1[2] + alphas1[2], sum(dat1[3:5]) + betas1[2])
                         z3 <- rbeta(100000, dat1[3] + alphas1[3], sum(dat1[4:5]) + betas1[3])
                         z4 <- rbeta(100000, dat1[4] + alphas1[4], dat1[5]  + betas1[4])
                         
                         ## convert from the Z-scale to the p-scale
                         p2 <- z2*(1-z1); p3 <- z3*(1-z1)*(1-z2); p4 <- z4*(1-z3)*(1-z1)*(1-z2)
                         p5 <- 1 - z1 - p2 - p3 - p4
                         
                         ## get posterior draws for theta1
                         theta.1 <- z1 + 2*p2 + 3*p3 + 4*p4 + 5*p5
                         
                         ## repeat for model 2
                         dat2 <- rmultinom(1, round(q*samps[k]), prob = as.numeric(paramss2[j,]))
                         z1 <- rbeta(100000, dat2[1] + alphas2[1], sum(dat2[2:5]) + betas2[1])
                         z2 <- rbeta(100000, dat2[2] + alphas2[2], sum(dat2[3:5]) + betas2[2])
                         z3 <- rbeta(100000, dat2[3] + alphas2[3], sum(dat2[4:5]) + betas2[3])
                         z4 <- rbeta(100000, dat2[4] + alphas2[4], dat2[5]  + betas2[4])
                         
                         p2 <- z2*(1-z1); p3 <- z3*(1-z1)*(1-z2); p4 <- z4*(1-z3)*(1-z1)*(1-z2)
                         p5 <- 1 - z1 - p2 - p3 - p4
                         
                         theta.2 <- z1 + 2*p2 + 3*p3 + 4*p4 + 5*p5
                         
                         ## output posterior probability, where delta_U is Inf
                         mean(ifelse(theta.1 - theta.2 >= delta_L, 1, 0))
                       }
  
  ## output results to .csv file
  write.csv(probs_out, paste0("probs_red_sec5_", samps[k], ".csv"), row.names = FALSE)
}

## create matrices for left contour plot in Figure 3 (based on the first
## sample size calculation)
first_rep <- as.numeric(read.csv("sec4p3_sobol_8192_summary.csv")[1,])
n_low <- first_rep[7]; n_mid <- first_rep[8];
n_high <- first_rep[9]

## read in the posterior probabilities corresponding to n^((0)) and n^((1))
green_mid <- unlist(read.csv("sec4p3_sobol_8192_green_mid.csv")[1,])
red_mid <- unlist(read.csv("sec4p3_sobol_8192_red_mid.csv")[1,])
green_low <- unlist(read.csv("sec4p3_sobol_8192_green_small.csv")[1,])
red_low <- unlist(read.csv("sec4p3_sobol_8192_red_small.csv")[1,])

## get the slopes and intercepts for the linear approximations
green_slope <- (green_mid - green_low)/(n_mid-n_low)
red_slope <- (red_mid - red_low)/(n_mid-n_low)
green_int <- green_low - green_slope*n_low
red_int <- red_low - red_slope*n_low

## approximate the sampling distributions of posterior probabilities
## on the logit scale using these approximations
for (i in seq(100, as.numeric(n_mid))){
  assign(paste0("green_vec_", i), green_int + green_slope*i)
  assign(paste0("red_vec_", i), red_int + red_slope*i)
}

## read in the posterior probabilities corresponding to n^((2))
green_high <- unlist(read.csv("sec4p3_sobol_8192_green_large.csv")[1,])
red_high <- unlist(read.csv("sec4p3_sobol_8192_red_large.csv")[1,])

## get the slopes and intercepts for the linear approximations
## use n^((2)) instead of n^((0)) this time
green_slope <- (green_high - green_mid)/(n_high-n_mid)
red_slope <- (red_high - red_mid)/(n_high-n_mid)
green_int <- green_mid - green_slope*n_mid
red_int <- red_mid - red_slope*n_mid

## approximate the sampling distributions of posterior probabilities
## on the logit scale using these approximations
for (i in seq(as.numeric(n_mid) + 1, 125)){
  assign(paste0("green_vec_", i), green_int + green_slope*i)
  assign(paste0("red_vec_", i), red_int + red_slope*i)
}

## create a vector of gamma values on the logit scale to compute power
## and type I error rate estimates
opt_gamma <- as.numeric(first_rep[3])
opt_gamma <- log(opt_gamma) - log(1 - opt_gamma)

gammas <- seq(log(0.91) - log(0.09), log(0.96) - log(0.04), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))


x <- seq(100, 120, 1)
y <- 1/(1 + exp(-gammas))

## z matrix is for power
z_mat <- NULL
for (i in 1:length(x)){
  z_mat <- rbind(z_mat, get(paste0("green_vec_", x[i])))
}
z <- NULL
for (j in 1:length(y)){
  z <- cbind(z, rowMeans(z_mat > gammas[j]))
}

## w matrix is for type I error
w_mat <- NULL
for (i in 1:length(x)){
  w_mat <- rbind(w_mat, get(paste0("red_vec_", x[i])))
}
w <- NULL
for (j in 1:length(y)){
  w <- cbind(w, rowMeans(w_mat > gammas[j]))
}

write.csv(w, "w_mat1.csv", row.names = FALSE)
write.csv(z, "z_mat1.csv", row.names = FALSE)

## create matrices for centre contour plot in Figure 3 (based on the average
## of 1000 sample size calculations)

## read in the sample sizes and posterior probabilities (logit scale) for 
## all 1000 repetitions
n_lows <- as.numeric(read.csv("sec4p3_sobol_8192_summary.csv")[,7])
n_mids <- as.numeric(read.csv("sec4p3_sobol_8192_summary.csv")[,8])
n_highs <- as.numeric(read.csv("sec4p3_sobol_8192_summary.csv")[,9])

green_mids <- read.csv("sec4p3_sobol_8192_green_mid.csv")
red_mids <- read.csv("sec4p3_sobol_8192_red_mid.csv")
green_lows <- read.csv("sec4p3_sobol_8192_green_small.csv")
red_lows <- read.csv("sec4p3_sobol_8192_red_small.csv")
green_highs <- read.csv("sec4p3_sobol_8192_green_large.csv")
red_highs <- read.csv("sec4p3_sobol_8192_red_large.csv")

z_full <- matrix(0, nrow = 21, ncol = 50)
w_full <- matrix(0, nrow = 21, ncol = 50)
for (k in 1:1000){
  ## the process below repeats the process detailed for sample size calculation
  ## 1 above for all repetitions k = {1, 2, ..., 1000}
  green_low <- unlist(green_lows[k,]); green_mid <- unlist(green_mids[k,]); 
  green_high <- unlist(green_highs[k,])
  
  red_low <- unlist(red_lows[k,]); red_mid <- unlist(red_mids[k,]); 
  red_high <- unlist(red_highs[k,])
  
  n_low <- n_lows[k]; n_mid <- n_mids[k]; n_high <- n_highs[k]
  
  green_slope <- (green_mid - green_low)/(n_mid-n_low)
  red_slope <- (red_mid - red_low)/(n_mid-n_low)
  green_int <- green_low - green_slope*n_low
  red_int <- red_low - red_slope*n_low
  
  for (i in seq(100, as.numeric(n_mid))){
    assign(paste0("green_vec_", i), green_int + green_slope*i)
    assign(paste0("red_vec_", i), red_int + red_slope*i)
  }
  
  green_slope <- (green_high - green_mid)/(n_high-n_mid)
  red_slope <- (red_high - red_mid)/(n_high-n_mid)
  green_int <- green_mid - green_slope*n_mid
  red_int <- red_mid - red_slope*n_mid
  
  for (i in seq(as.numeric(n_mid) + 1, 125)){
    assign(paste0("green_vec_", i), green_int + green_slope*i)
    assign(paste0("red_vec_", i), red_int + red_slope*i)
  }
  
  gammas <- seq(log(0.91) - log(0.09), log(0.96) - log(0.04), length.out = 50)
  
  x <- seq(100, 120, 1)
  y <- 1/(1 + exp(-gammas))
  
  z_mat <- NULL
  for (i in 1:length(x)){
    z_mat <- rbind(z_mat, get(paste0("green_vec_", x[i])))
  }
  z <- NULL
  for (j in 1:length(y)){
    z <- cbind(z, rowMeans(z_mat > gammas[j]))
  }
  
  ## now we multiply the matrix for each sample size calculation by 1/1000
  ## to get the matrices corresponding to the average
  z_full <- 0.001*z + z_full
  
  w_mat <- NULL
  for (i in 1:length(x)){
    w_mat <- rbind(w_mat, get(paste0("red_vec_", x[i])))
  }
  w <- NULL
  for (j in 1:length(y)){
    w <- cbind(w, rowMeans(w_mat > gammas[j]))
  }
  
  w_full <- 0.001*w + w_full
  
}
write.csv(z_full, "z_full_mat.csv", row.names = FALSE)
write.csv(w_full, "w_full_mat.csv", row.names = FALSE)

## create matrices for right contour plot in Figure 3 (based on
## simulating multinomial data)
z_full2 <- matrix(0, nrow = 21, ncol = 50)
w_full2 <- matrix(0, nrow = 21, ncol = 50)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in seq(100, 120)){
    assign(paste0("green_vec_", i), 
           pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_green_sec5_", i,".csv")))),
                     .Machine$double.eps), 1 - 10^(-7)))
    assign(paste0("red_vec_", i), 
           pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_red_sec5_", i,".csv")))),
                     .Machine$double.eps), 1 - 10^(-7)))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.91) - log(0.09), log(0.96) - log(0.04), length.out = 50)
  
x <- seq(100, 120, 1)
y <- 1/(1 + exp(-gammas))
  
z_mat <- NULL
for (i in 1:length(x)){
  z_mat <- rbind(z_mat, get(paste0("green_vec_", x[i])))
}
z_mat <- log(z_mat) - log(1-z_mat)
z <- NULL
for (j in 1:length(y)){
  z <- cbind(z, rowMeans(z_mat > gammas[j]))
}
  
z_full2 <- z
  
w_mat <- NULL
for (i in 1:length(x)){
  w_mat <- rbind(w_mat, get(paste0("red_vec_", x[i])))
}
w_mat <- log(w_mat) - log(1-w_mat)
w <- NULL
for (j in 1:length(y)){
  w <- cbind(w, rowMeans(w_mat > gammas[j]))
}
  
w_full2 <- w

## write output to a .csv file  
write.csv(z_full2, "z_full_mat2.csv", row.names = FALSE)
write.csv(w_full2, "w_full_mat2.csv", row.names = FALSE)

## create the three contour plots and output as .pdf file for the article
pdf(file = "Figure3OC.pdf",   # The directory you want to save the file in
    width = 7.5, 
    height = 5) 

par(mfrow=c(2,3), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.25,1,0))

## read in matrices for left plot
z <- matrix(unlist(read.csv("z_mat1.csv")), nrow = 21, ncol = 51)
w <- matrix(unlist(read.csv("w_mat1.csv")), nrow = 21, ncol = 51)
gammas <- seq(log(0.91) - log(0.09), log(0.96) - log(0.04), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))
y <- 1/(1 + exp(-gammas))

contour(x, y, w, levels = c(seq(0.02, 0.035, 0.005), 0.045, seq(0.065, 0.08, 0.005)), 
        xlab = expression(italic("n")), xlim = c(100,120),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
points(x = first_rep[2], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))
contour(x, y, w, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, w, levels = c(0.055), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w, levels = c(0.06), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w, levels = c(0.04), col = "black", add = TRUE, labcex = 0.8)
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
box()

gammas <- seq(log(0.91) - log(0.09), log(0.96) - log(0.04), length.out = 50)
y <- 1/(1 + exp(-gammas))

z_full <- matrix(unlist(read.csv("z_full_mat.csv")), nrow = 21, ncol = 50)
w_full <- matrix(unlist(read.csv("w_full_mat.csv")), nrow = 21, ncol = 50)

contour(x, y, w_full, levels = seq(0.02, 0.08, 0.005), xlab = expression(italic("n")), 
        xlim = c(100,120), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
box()

z_full2 <- matrix(unlist(read.csv("z_full_mat2.csv")), nrow = 21, ncol = 50)
w_full2 <- matrix(unlist(read.csv("w_full_mat2.csv")), nrow = 21, ncol = 50)

contour(x, y, w_full2, levels = c(seq(0.02, 0.045, 0.005),seq(0.055, 0.08, 0.005)), 
        xlab = expression(italic("n")), xlim = c(100,120), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
box()

gammas <- seq(log(0.91) - log(0.09), log(0.96) - log(0.04), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))
y <- 1/(1 + exp(-gammas))

contour(x, y, z, levels = c(seq(0.68, 0.77, 0.015), 0.80, seq(0.83, 0.86, 0.015)), 
        xlab = expression(italic("n")), xlim = c(100,120),axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
points(x = first_rep[2], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))
contour(x, y, w, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8)
contour(x, y, z, levels = c(0.785), col = "black", add = TRUE, method = "edge", labcex = 0.8)
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
box()

gammas <- seq(log(0.91) - log(0.09), log(0.96) - log(0.04), length.out = 50)
y <- 1/(1 + exp(-gammas))

contour(x, y, z_full, levels = c(seq(0.68, 0.8, 0.015), seq(0.83, 0.86, 0.015)), 
        xlab = expression(italic("n")), xlim = c(100,120), axes = FALSE,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge", cex.lab = 1.25)
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
box()
contour(x, y, w_full, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8)
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
box()

contour(x, y, z_full2, levels = c(seq(0.68, 0.86, 0.015)), 
        xlab = expression(italic("n")), xlim = c(100,120), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8)
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()

## Generate the density plots for Figure C.1 to compare using Sobol' 
## sequences with pseudorandom ones

## first create density plot for the sample size recommendations
dSOB <- density(temp[,2])
dPRNG <- density(tempPRNG[,2])
dPRNG24k <- density(tempPRNG24k[,2])
d_frame1 <- data.frame(x = dSOB$x, y = dSOB$y, type = "ASOB")
d_frame1B <- data.frame(x = dPRNG$x, y = dPRNG$y, type = "BPRNG")
d_frame1C <- data.frame(x = dPRNG24k$x, y = dPRNG24k$y, type = "CPRNG24k")

plotC1a <- ggplot(data=d_frame1, aes(x=x, y= y)) + theme_bw() +
  geom_polygon(aes(y=y), col="gray10", fill="snow1", size=0.75, alpha=0) +
  geom_polygon(aes(x = d_frame1B$x, y= d_frame1B$y), col="#077DAA", fill="snow1", size=0.75, alpha=0) +
  geom_polygon(aes(x = d_frame1C$x, y= d_frame1C$y), col="#f49b3f", fill="snow1", size=0.75, alpha=0) +
  geom_polygon(aes(y=y), col="gray10", fill="snow1", size=0.75, alpha=0) +
  labs(x=bquote(italic('n')), y=" ") + 
  theme(plot.title = element_text(margin = unit(c(0, 0, 5, 0), "mm"),
                                  hjust = 0.5,size=18,face="bold", colour = "#FFFFFF")) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=16,face="bold")) + 
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))

## create density plot for the conviction thresholds
ddSOB <- density(temp[,3])
ddPRNG <- density(tempPRNG[,3])
ddPRNG24k <- density(tempPRNG24k[,3])
d_frame2 <- data.frame(x = ddSOB$x, y = ddSOB$y, type = "ASOB")
d_frame2B <- data.frame(x = ddPRNG$x, y = ddPRNG$y, type = "BPRNG")
d_frame2C <- data.frame(x = ddPRNG24k$x, y = ddPRNG24k$y, type = "CPRNG24k")

## create density plot for induced prior on theta1
plotC1b <- ggplot(data=d_frame2, aes(x=x, y= y)) + theme_bw() +
  geom_polygon(aes(y=y), col="gray10", fill="snow1", size=0.75, alpha=0) +
  geom_polygon(aes(x = d_frame2B$x, y= d_frame2B$y), col="#077DAA", fill="snow1", size=0.75, alpha=0) +
  geom_polygon(aes(x = d_frame2C$x, y= d_frame2C$y), col="#f49b3f", fill="snow1", size=0.75, alpha=0) +
  geom_polygon(aes(y=y), col="gray10", fill="snow1", size=0.75, alpha=0) +
  labs(x=bquote(gamma), y=" ") + 
  theme(plot.title = element_text(margin = unit(c(0, 0, 5, 0), "mm"),
                                  hjust = 0.5,size=18,face="bold", colour = "#FFFFFF")) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=16,face="bold")) + 
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))

## create common legend for the plots
dcomb <- rbind(d_frame2, d_frame2B, d_frame2C)
plotC1c <- ggplot(data=dcomb, aes(x=x, y= y)) + theme_bw() +
  geom_line(aes(y = y, color=as.factor(type)), size = 1) +
  theme(legend.position="bottom") +
  scale_color_manual(name = " ", 
                     c(expression("Sobol' ("*italic(m)*" = 8192)  "), 
                       expression("PRNG ("*italic(m)*" = 8192)  "), 
                       expression("PRNG ("*italic(m)*" = 24000)")),
                     values = c("black", "#077DAA", "#f49b3f")) +
  theme(legend.text=element_text(size=16)) 

mylegend <- get_legend(plotC1c)

## combine two plots
fig_pre <- plot_grid(plotC1a + theme(plot.margin=unit(c(0.125,0.25,0.125,0.25),"cm")), 
                     plotC1b + theme(plot.margin=unit(c(0.125,0.25,0.125,0.25),"cm")), 
                     ncol = 2, rel_widths = c(1,1))

fig_final <- plot_grid(fig_pre, mylegend, ncol = 1, rel_heights = c(2, .3))

# output as .pdf file for the article
pdf(file = "FigureC1OC.pdf",   # The directory you want to save the file in
    width = 7.75, # The width of the plot in inches
    height = 4.25) # The height of the plot in inches

fig_final

dev.off()