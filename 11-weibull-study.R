## code to produce figures in multinomial numerical study, where Algorithm 3 is used
## instead of Algorithm 1 (Figure 3 and the figure in Appendix C.2)

## run this file AFTER running file 09-prior-elicitation.R and 
## 10-functions-for-app-d2.R

## BEGIN SETUP ##

## load necessary packages
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(colorspace)

## save these parameter values to one matrix for each region
greens <- read.csv("weibull_green.csv")
reds <- read.csv("weibull_red.csv")
greens <- matrix(unlist(greens), ncol = 4)
reds <- matrix(unlist(reds), ncol = 4)

## define hyperparameters for analysis priors
alphas1 <- c(2, 1)
betas1 <- c(2, 1)
alphas2 <- c(2, 1)
betas2 <- c(2, 1)

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## repeat the sample size calculation for the Weibull example 1000 times
## with the different Sobol' sequences
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <-  findGamma(green_par = greens, red_par = reds, 
                                         pwr = 0.7, typeI = 0.1, deltas = c(log(1/1.2), log(1.2)), q = 1, 
                                         alphas1, betas1, alphas2, betas2, m = 8192,
                                         seed_green = i + 8000, seed_red = i + 9000, contour = TRUE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp[,1:9], "wei_sobol_8192_summary.csv", row.names = FALSE)

## output the red and green posterior probabilities for each of the three sample
## sizes n^((0)) = small, n^((1)) = mid, and n^((2)) = large
write.csv(temp[,10:8201], "wei_sobol_8192_red_small.csv", row.names = FALSE)
write.csv(temp[,8202:16393], "wei_sobol_8192_green_small.csv", row.names = FALSE)
write.csv(temp[,16394:24585], "wei_sobol_8192_red_mid.csv", row.names = FALSE)
write.csv(temp[,24586:32777], "wei_sobol_8192_green_mid.csv", row.names = FALSE)
write.csv(temp[,32778:40969], "wei_sobol_8192_red_large.csv", row.names = FALSE)
write.csv(temp[,40970:49161], "wei_sobol_8192_green_large.csv", row.names = FALSE)

## repeat the 1000 sample size calculations with the SAME Sobol' sequences but
## explore the entire hypercube and each sample size considered
tempFull <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                    .options.snow=opts, .errorhandling = "remove") %dopar% {
                      temp_res <- findGammaFull(green_par = greens, red_par = reds, 
                                                pwr = 0.7, typeI = 0.1, deltas = c(log(1/1.2), log(1.2)), q = 1, 
                                                alphas1, betas1, alphas2, betas2, m = 8192,
                                                seed_green = i + 8000, seed_red = i + 9000)
                      
                      c(i, as.numeric(unlist(temp_res)))
                    }

write.csv(tempFull, "wei_sobol_8192_full_summary.csv", row.names = FALSE)

## create matrices for left contour plot in Figure D.3 (based on the average
## of 1000 sample size calculations); we do not have a contour plot for the
## single sample size calculation in this figure

## read in the sample sizes and posterior probabilities (logit scale) for 
## all 1000 repetitions
n_lows <- as.numeric(read.csv("wei_sobol_8192_summary.csv")[,7])
n_mids <- as.numeric(read.csv("wei_sobol_8192_summary.csv")[,8])
n_highs <- as.numeric(read.csv("wei_sobol_8192_summary.csv")[,9])

green_mids <- read.csv("wei_sobol_8192_green_mid.csv")
red_mids <- read.csv("wei_sobol_8192_red_mid.csv")
green_lows <- read.csv("wei_sobol_8192_green_small.csv")
red_lows <- read.csv("wei_sobol_8192_red_small.csv")
green_highs <- read.csv("wei_sobol_8192_green_large.csv")
red_highs <- read.csv("wei_sobol_8192_red_large.csv")

z_full <- matrix(0, nrow = 26, ncol = 60)
w_full <- matrix(0, nrow = 26, ncol = 60)
for (k in 1:1000){
  if (k %% 25 == 0){print(k)}
  ## repeat the process from 04-multinomial-study-alg3.R with Weibull model
  green_low <- unlist(green_lows[k,]); green_mid <- unlist(green_mids[k,]); 
  green_high <- unlist(green_highs[k,])
  
  red_low <- unlist(red_lows[k,]); red_mid <- unlist(red_mids[k,]); 
  red_high <- unlist(red_highs[k,])
  
  n_low <- n_lows[k]; n_mid <- n_mids[k]; n_high <- n_highs[k]
  
  green_slope <- (green_mid - green_low)/(n_mid-n_low)
  red_slope <- (red_mid - red_low)/(n_mid-n_low)
  green_int <- green_low - green_slope*n_low
  red_int <- red_low - red_slope*n_low
  
  for (i in seq(145, as.numeric(n_mid))){
    assign(paste0("green_vec_", i), green_int + green_slope*i)
    assign(paste0("red_vec_", i), red_int + red_slope*i)
  }
  
  green_slope <- (green_high - green_mid)/(n_high-n_mid)
  red_slope <- (red_high - red_mid)/(n_high-n_mid)
  green_int <- green_mid - green_slope*n_mid
  red_int <- red_mid - red_slope*n_mid
  
  for (i in seq(as.numeric(n_mid) + 1, 170)){
    assign(paste0("green_vec_", i), green_int + green_slope*i)
    assign(paste0("red_vec_", i), red_int + red_slope*i)
  }
  
  gammas <- seq(log(0.84) - log(0.16), log(0.9) - log(0.1), length.out = 60)
  
  x <- seq(145, 170, 1)
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
write.csv(z_full[6:26,], "wei_z_full_mat.csv", row.names = FALSE)
write.csv(w_full[6:26,], "wei_w_full_mat.csv", row.names = FALSE)

## approximate the sampling distributions of posterior probabilities by
## generating data directly (for the contour plots)
## function to approximate posteriors using MCMC
WeibullPost <- function(y1, y2, mu1, tau1, kappa1, lambda1,
                        mu2, tau2, kappa2, lambda2, tau = 4.29, burnin = 1000, nchains = 1,
                        nthin = 1, ndraws = 5000){
  
  ## y1 and y2 are the food expenditure observations in each group
  ## lambda_j has a Gamma(mu_j, tau_j) prior, where tau_j is a rate
  ## k_j has a Gamma(kappa_j, lambda_j) prior, where lambda_j is a rate (nu_j in paper)
  ## tau is the threshold for the tail probability
  ## burnin is the number of MCMC iterations to discard at the start of each chain
  ## nchains is the number of chains to generate
  ## nthin is the thinning parameter for the MCMC process
  ## ndraws is the number of draws to generate (excluding burnin but including thinned draws)
  
  n1 <- length(y1)
  model1.fit <- jags.model(file="JAGS_weibull.txt",
                           data=list(n=n1, y = y1, 
                                     tau0 = tau1, mu0 = mu1,
                                     kappa0 = kappa1, lambda0 = lambda1), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model1.fit, burnin, progress.bar = "none")
  model1.samples <- coda.samples(model1.fit, c("nu", "l"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  nu.1 <- unlist(model1.samples[,2])
  lambda.1 <- unlist(model1.samples[,1])
  
  n2 <- length(y2)
  model2.fit <- jags.model(file="JAGS_weibull.txt",
                           data=list(n=n2, y = y2, 
                                     tau0 = tau2, mu0 = mu2,
                                     kappa0 = kappa2, lambda0 = lambda2), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model2.fit, burnin, progress.bar = "none")
  model2.samples <- coda.samples(model2.fit, c("nu", "l"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  nu.2 <- unlist(model2.samples[,2])
  lambda.2 <- unlist(model2.samples[,1])
  
  theta1 <- lambda.1*log(1/(1 - tau))^(1/nu.1) 
  theta2 <- lambda.2*log(1/(1 - tau))^(1/nu.2)
  theta <- theta1/theta2
  theta <- ifelse(is.na(theta), Inf, theta)
  return(theta)
}

resamps <- 40960
pb <- txtProgressBar(max = resamps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get parameter values from design priors to generate
## the Weibull data for the green region
samps <- seq(150, 170, 2)
params <- rbind(greens, greens, greens, greens, greens)

delta_L <- 1/1.2
delta_U <- 1.2

for (k in 1:length(samps)){
  probs_in <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        y_star1 <- rweibull(samps[k], params[j,2], params[j,1])
                        y_star2 <- rweibull(samps[k], params[j,4], params[j,3])
                        
                        theta.diff <- WeibullPost(y_star1, y_star2, 2, 1,
                                                  2, 1, 2, 1, 2, 1, 0.9)
                        
                        mean(ifelse(theta.diff > delta_L, theta.diff < delta_U,0))
                      }
  
  ## output results to .csv file
  write.csv(probs_in, paste0("probs_green_appd2_", samps[k], ".csv"), row.names = FALSE)
}

## get parameter values from design priors to generate
## the Weibull data for the red region
params <- rbind(reds, reds, reds, reds, reds)

for (k in 1:length(samps)){
  
  probs_out <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts, .errorhandling = "remove") %dopar% {
                         
                         y_star1 <- rweibull(samps[k], params[j,2], params[j,1])
                         y_star2 <- rweibull(samps[k], params[j,4], params[j,3])
                         
                         theta.diff <- WeibullPost(y_star1, y_star2, 2, 1,
                                                   2, 1, 2, 1, 2, 1, 0.9)
                         
                         mean(ifelse(theta.diff > delta_L, theta.diff < delta_U,0))
                       }
  
  ## output results to .csv file
  write.csv(probs_out, paste0("probs_red_appd2_", samps[k], ".csv"), row.names = FALSE)
}

## create matrices for right contour plot in Figure D.3 (based on
## Weibull data simulated above)
z_full2 <- matrix(0, nrow = 11, ncol = 60)
w_full2 <- matrix(0, nrow = 11, ncol = 60)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in seq(150,170,2)){
  assign(paste0("green_vec_", i), 
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_green_appd2_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
  assign(paste0("red_vec_", i),
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_red_appd2_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.84) - log(0.16), log(0.9) - log(0.1), length.out = 60)

x <- seq(150, 170, 2)
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
write.csv(z_full2, "wei_z_full_mat2.csv", row.names = FALSE)
write.csv(w_full2, "wei_w_full_mat2.csv", row.names = FALSE)

## create the three contour plots and output as .pdf file for the article
pdf(file = "FigureD3OC2.pdf",   # The directory you want to save the file in
    width = 6.5, 
    height = 6.5) 

par(mfrow=c(2,2), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.25,1,0))

x <- seq(150, 170, 1)
gammas <- seq(log(0.84) - log(0.16), log(0.9) - log(0.1), length.out = 60)
y <- 1/(1 + exp(-gammas))

z_full <- matrix(unlist(read.csv("wei_z_full_mat.csv")), nrow = 21, ncol = 60)
w_full <- matrix(unlist(read.csv("wei_w_full_mat.csv")), nrow = 21, ncol = 60)

contour(x, y, w_full, levels = c(seq(0.06, 0.08, 0.01), 0.1, 0.12, 0.13), xlab = expression(italic("n")), 
        xlim = c(150,170), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.1), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.7), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(150, 170, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.84, 0.9, 0.02), cex.axis = 1.15)
box()
contour(x, y, w_full, levels = c(0.09), col = "black", add = TRUE, labcex = 0.8)
contour(x, y, w_full, levels = c(0.11), col = "black", add = TRUE, labcex = 0.8)

z_full2 <- matrix(unlist(read.csv("wei_z_full_mat2.csv")), nrow = 11, ncol = 60)
w_full2 <- matrix(unlist(read.csv("wei_w_full_mat2.csv")), nrow = 11, ncol = 60)

x <- seq(150, 170, 2)

contour(x, y, w_full2, levels = c(seq(0.06, 0.08, 0.01), 0.1, 0.12, 0.13),
        xlab = expression(italic("n")), xlim = c(150,170), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.1), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.7), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(150, 170, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.84, 0.9, 0.02), cex.axis = 1.15)
box()
contour(x, y, w_full2, levels = c(0.09), col = "black", add = TRUE, labcex = 0.8)
contour(x, y, w_full2, levels = c(0.11), col = "black", add = TRUE, labcex = 0.8)

contour(x, y, z_full, levels = c(seq(0.56, 0.8, 0.02)), 
        xlab = expression(italic("n")), xlim = c(150,170), axes = FALSE,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge", cex.lab = 1.25)
axis(side = 1, at = seq(150, 170, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.84, 0.9, 0.02), cex.axis = 1.15)
box()
contour(x, y, w_full, levels = c(0.1), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full, levels = c(0.7), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")

x <- seq(150, 170, 1)
contour(x, y, z_full, levels = c(seq(0.56, 0.8, 0.02)), 
        xlab = expression(italic("n")), xlim = c(150,170), axes = FALSE,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge", cex.lab = 1.25)
axis(side = 1, at = seq(150, 170, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.84, 0.9, 0.02), cex.axis = 1.15)
box()
contour(x, y, w_full, levels = c(0.1), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full, levels = c(0.7), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")

x <- seq(150, 170, 2)
contour(x, y, z_full2, levels = c(seq(0.56, 0.8, 0.02)),
        xlab = expression(italic("n")), xlim = c(150,170), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.1), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge")
contour(x, y, z_full2, levels = c(0.7), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(150, 170, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.84, 0.9, 0.02), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()