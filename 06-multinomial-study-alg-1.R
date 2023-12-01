## code to produce figures in multinomial numerical study, where Algorithm 3 is used
## instead of Algorithm 1 (Figure 3 and the figure in Appendix C.2)

## run this file AFTER running file 05-functions-for-app-c3.R

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

## the commented out code here only needs to be run if the file
## 04-multinomial-study-alg-3.R has not already been run
## we first generate the eta values from the design priors
## extract design prior parameters for groups 1 and 2
## from the Shiny app
# alphas1 <- read.csv("priors_group1.csv")$alpha
# betas1 <- read.csv("priors_group1.csv")$beta
# 
# alphas2 <- read.csv("priors_group2.csv")$alpha
# betas2 <- read.csv("priors_group2.csv")$beta
# 
# ## generate observations from the design priors
# z11 = rbeta(5000000, alphas1[1], betas1[1]); z12 = rbeta(5000000, alphas1[2], betas1[2])
# z13 = rbeta(5000000, alphas1[3], betas1[3]); z14 = rbeta(5000000, alphas1[4], betas1[4])
# z21 = rbeta(5000000, alphas2[1], betas2[1]); z22 = rbeta(5000000, alphas2[2], betas2[2])
# z23 = rbeta(5000000, alphas2[3], betas2[3]); z24 = rbeta(5000000, alphas2[4], betas2[4])
# 
# ## convert from the Z-scale to the p-scale for both groups
# p11 <- z11
# p12 <- z12*(1-z11)
# p13 <- z13*(1-z12)*(1-z11)
# p14 <- z14*(1-z13)*(1-z12)*(1-z11)
# 
# p21 <- z21
# p22 <- z22*(1-z21)
# p23 <- z23*(1-z22)*(1-z21)
# p24 <- z24*(1-z23)*(1-z22)*(1-z21)
# 
# 
# ## tweak any multinomial probabilities that round to 0 or 1
# ## (to four decimal places); the "problem_inds" are the indices
# ## that we need to take care of
# p11 <- pmax(z11, 0.001)
# p12 <- pmax(z12*(1-z11), 0.001)
# p13 <- pmax(z13*(1-z12)*(1-z11), 0.001)
# p14 <- pmax(z14*(1-z13)*(1-z12)*(1-z11), 0.001)
# 
# ## slightly rescale if necessary for group 1
# problem_inds <- which(p11 + p12 + p13 + p14 > 1 - 0.001)
# sums1 <- p11 + p12 + p13 + p14
# 
# if (length(problem_inds) > 0){
#   for (i in 1:length(problem_inds)){
#     p11[i] <- (1 - 0.001)*p11[i]/sums1[problem_inds[i]]
#     p12[i] <- (1 - 0.001)*p12[i]/sums1[problem_inds[i]]
#     p13[i] <- (1 - 0.001)*p13[i]/sums1[problem_inds[i]]
#     p14[i] <- (1 - 0.001)*p14[i]/sums1[problem_inds[i]]
#   }
# }
# 
# ## repeat this process for the second group
# p21 <- pmax(z21, 0.001)
# p22 <- pmax(z22*(1-z21), 0.001)
# p23 <- pmax(z23*(1-z22)*(1-z21), 0.001)
# p24 <- pmax(z24*(1-z23)*(1-z22)*(1-z21), 0.001)
# 
# ## slightly rescale if necessary for group 2
# problem_inds <- which(p21 + p22 + p23 + p24 > 1 - 0.001)
# sums2 <- p21 + p22 + p23 + p24 
# 
# if (length(problem_inds) > 0){
#   for (i in 1:length(problem_inds)){
#     p21[problem_inds[i]] <- (1 - 0.001)*p21[problem_inds[i]]/sums2[problem_inds[i]]
#     p22[problem_inds[i]] <- (1 - 0.001)*p22[problem_inds[i]]/sums2[problem_inds[i]]
#     p23[problem_inds[i]] <- (1 - 0.001)*p23[problem_inds[i]]/sums2[problem_inds[i]]
#     p24[problem_inds[i]] <- (1 - 0.001)*p24[problem_inds[i]]/sums2[problem_inds[i]]
#   }
# }
# 
# ## obtain theta values for each draw to segment the priors
# theta1 <- p11 + 2*p12 + 3*p13 + 4*p14 + 5*(1 - p11 - p12 - p13 - p14)
# theta2 <- p21 + 2*p22 + 3*p23 + 4*p24 + 5*(1 - p21 - p22 - p23 - p24)
# 
# ## define the green and red regions
# delta_L <- -0.5; green_mid <- -0.2; green_hw <- 0.1
# green_low <- green_mid - green_hw; green_high <- green_mid + green_hw
# red_high <- delta_L; red_low <- red_high - 0.05
# 
# ## approximate the induced prior density function on theta for sampling-resampling
# diff_dens <- density(theta1 - theta2)
# diff_pdf <- approxfun(x = diff_dens$x,
#                       y = diff_dens$y)
# 
# ## conduct sampling-resampling to uniformly sample from the green region
# diff_num <- ifelse(theta1 - theta2 <= green_high, theta1 - theta2 >= green_low, 0)
# diff_denom <- diff_pdf(theta1 - theta2)
# 
# weights <- diff_num/diff_denom
# weights <- weights/sum(weights)
# 
# ## resample to get 8192 parameter values from green region
# inds_region <- sample(1:5000000, size = 8192, replace = TRUE, prob = weights)
# 
# ## output parameter values to a .csv file
# write.csv(cbind(z11, z12, z13, z14)[inds_region,], "zs1_green.csv", row.names = FALSE)
# write.csv(cbind(p11, p12, p13, p14)[inds_region,], "ps1_green.csv", row.names = FALSE)
# write.csv(theta1[inds_region], "theta1_green.csv", row.names = FALSE)
# 
# write.csv(cbind(z21, z22, z23, z24)[inds_region,], "zs2_green.csv", row.names = FALSE)
# write.csv(cbind(p21, p22, p23, p24)[inds_region,], "ps2_green.csv", row.names = FALSE)
# write.csv(theta2[inds_region], "theta2_green.csv", row.names = FALSE)
# 
# ## conduct sampling-resampling to uniformly sample from the red region
# diff_num <- ifelse(theta1 - theta2 <= red_high, theta1 - theta2 >= red_low, 0)
# diff_denom <- diff_pdf(theta1 - theta2)
# 
# weights <- diff_num/diff_denom
# weights <- weights/sum(weights)
# 
# ## resample to get 8192 parameter values from red region
# inds_region2 <- sample(1:5000000, size = 8192, replace = TRUE, prob = weights)
# 
# ## output parameter values to a .csv file
# write.csv(cbind(z11, z12, z13, z14)[inds_region2,], "zs1_red.csv", row.names = FALSE)
# write.csv(cbind(p11, p12, p13, p14)[inds_region2,], "ps1_red.csv", row.names = FALSE)
# write.csv(theta1[inds_region2], "theta1_red.csv", row.names = FALSE)
# 
# write.csv(cbind(z21, z22, z23, z24)[inds_region2,], "zs2_red.csv", row.names = FALSE)
# write.csv(cbind(p21, p22, p23, p24)[inds_region2,], "ps2_red.csv", row.names = FALSE)
# write.csv(theta2[inds_region2], "theta2_red.csv", row.names = FALSE)

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
                                        seed_green = i + 6000, seed_red = i + 7000, contour = TRUE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp[,1:9], "alg1_sobol_8192_summary.csv", row.names = FALSE)

## output the red and green posterior probabilities for each of the three sample
## sizes n^((0)) = small, n^((1)) = mid, and n^((2)) = large
write.csv(temp[,10:8201], "alg1_sobol_8192_red_small.csv", row.names = FALSE)
write.csv(temp[,8202:16393], "alg1_sobol_8192_green_small.csv", row.names = FALSE)
write.csv(temp[,16394:24585], "alg1_sobol_8192_red_mid.csv", row.names = FALSE)
write.csv(temp[,24586:32777], "alg1_sobol_8192_green_mid.csv", row.names = FALSE)
write.csv(temp[,32778:40969], "alg1_sobol_8192_red_large.csv", row.names = FALSE)
write.csv(temp[,40970:49161], "alg1_sobol_8192_green_large.csv", row.names = FALSE)

## repeat the 1000 sample size calculations with the SAME Sobol' sequences but
## explore the entire hypercube and each sample size considered
tempFull <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGammaFull(green_par = greens, red_par = reds, 
                                            pwr = 0.8, typeI = 0.05, deltas = c(-0.5, 4), q = 1.25, 
                                            alphas1, betas1, alphas2, betas2, m = 8192,
                                            seed_green = i + 6000, seed_red = i + 7000)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

write.csv(tempFull, "alg1_sobol_8192_full_summary.csv", row.names = FALSE)

## create matrices for left contour plot in Figure C.2 (based on the average
## of 1000 sample size calculations); we do not have a contour plot for the
## single sample size calculation in this figure

## read in the sample sizes and posterior probabilities (logit scale) for 
## all 1000 repetitions
n_lows <- as.numeric(read.csv("alg1_sobol_8192_summary.csv")[,7])
n_mids <- as.numeric(read.csv("alg1_sobol_8192_summary.csv")[,8])
n_highs <- as.numeric(read.csv("alg1_sobol_8192_summary.csv")[,9])

green_mids <- read.csv("alg1_sobol_8192_green_mid.csv")
red_mids <- read.csv("alg1_sobol_8192_red_mid.csv")
green_lows <- read.csv("alg1_sobol_8192_green_small.csv")
red_lows <- read.csv("alg1_sobol_8192_red_small.csv")
green_highs <- read.csv("alg1_sobol_8192_green_large.csv")
red_highs <- read.csv("alg1_sobol_8192_red_large.csv")

z_full <- matrix(0, nrow = 47, ncol = 60)
w_full <- matrix(0, nrow = 47, ncol = 60)
for (k in 1:1000){
  if (k %% 25 == 0){print(k)}
  ## repeat the process from 04-multinomial-study-alg3.R 
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
  
  for (i in seq(as.numeric(n_mid) + 1, 146)){
    assign(paste0("green_vec_", i), green_int + green_slope*i)
    assign(paste0("red_vec_", i), red_int + red_slope*i)
  }
  
  gammas <- seq(log(0.91) - log(0.09), log(0.97) - log(0.03), length.out = 60)
  
  x <- seq(100, 146, 1)
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
write.csv(z_full[1:41,], "alg1_z_full_mat.csv", row.names = FALSE)
write.csv(w_full[1:41,], "alg1_w_full_mat.csv", row.names = FALSE)

## create matrices for right contour plot in Figure C.2 (based on
## the multinomial data simulated in file 04-multinomial-study-alg3.R)
z_full2 <- matrix(0, nrow = 41, ncol = 60)
w_full2 <- matrix(0, nrow = 41, ncol = 60)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in seq(100, 140)){
    assign(paste0("green_vec_", i), 
           pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_green_sec5_", i,".csv")))),
                     .Machine$double.eps), 1 - 10^(-7)))
    assign(paste0("red_vec_", i), 
           pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_red_sec5_", i,".csv")))),
                     .Machine$double.eps), 1 - 10^(-7)))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.91) - log(0.09), log(0.97) - log(0.03), length.out = 60)
  
x <- seq(100, 140, 1)
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
write.csv(z_full2, "alg1_z_full_mat2.csv", row.names = FALSE)
write.csv(w_full2, "alg1_w_full_mat2.csv", row.names = FALSE)

## create the three contour plots and output as .pdf file for the article
pdf(file = "FigureC2OC.pdf",   # The directory you want to save the file in
    width = 6.5, 
    height = 6.5) 

par(mfrow=c(2,2), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.25,1,0))

gammas <- seq(log(0.91) - log(0.09), log(0.97) - log(0.03), length.out = 60)
y <- 1/(1 + exp(-gammas))

z_full <- matrix(unlist(read.csv("alg1_z_full_mat.csv")), nrow = 41, ncol = 60)
w_full <- matrix(unlist(read.csv("alg1_w_full_mat.csv")), nrow = 41, ncol = 60)

contour(x, y, w_full, levels = seq(0.02, 0.08, 0.005), xlab = expression(italic("n")), 
        xlim = c(100,140), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(100, 140, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.97, 0.015), cex.axis = 1.15)
box()

z_full2 <- matrix(unlist(read.csv("alg1_z_full_mat2.csv")), nrow = 41, ncol = 60)
w_full2 <- matrix(unlist(read.csv("alg1_w_full_mat2.csv")), nrow = 41, ncol = 60)

contour(x, y, w_full2, levels = c(seq(0.02, 0.045, 0.005),seq(0.055, 0.08, 0.005)), 
        xlab = expression(italic("n")), xlim = c(100,140), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(100, 140, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.97, 0.015), cex.axis = 1.15)
box()

contour(x, y, z_full, levels = c(seq(0.605, 0.8, 0.015), seq(0.83, 0.89, 0.015)), 
        xlab = expression(italic("n")), xlim = c(100,140), axes = FALSE,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge", cex.lab = 1.25)
axis(side = 1, at = seq(100, 140, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.97, 0.015), cex.axis = 1.15)
box()
contour(x, y, w_full, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8)

contour(x, y, z_full2, levels = c(seq(0.605, 0.89, 0.015)), 
        xlab = expression(italic("n")), xlim = c(100,140), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8)
axis(side = 1, at = seq(100, 140, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.97, 0.015), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()

## create Figure C3 to illustrate the differences between the sampling distributions 
## of the MLEs and their approximations

## generate MLEs for the first binomial distribution (i.e., whether or not
## an observation is assigned to category 1); use observed proportion of 
## observations in category 1 of comparison group for illustration
set.seed(1)
MLEs_exact <- rbinom(100000, 109, 3/108)
MLEs_round <- pmin(pmax(0.0001, MLEs_exact/109), 0.9999)

df_MLEs <- data.frame(x = log(MLEs_round) - log(1 - MLEs_round))
binwidth2 <- 0.35
bins2 <- seq(-9.65, -1.6, 0.35)

## generate MLEs according to the approximation in Algorithm 3 by inverting
## the CDF associated with the approximation X^*
fn_x_star <- approxfun(y = c(0,seq(1,109,1)-0.5,109), 
                       x = c(0,pbinom(0:109, 109, 3/108)))
set.seed(2)
x_stars <- fn_x_star(runif(10000000))
dens_xstar <- density(log(pmin(pmax(0.0001, x_stars/100), 0.9999)) - 
                log(1-pmin(pmax(0.0001, x_stars/100), 0.9999)))
d_xstar <- data.frame(x = dens_xstar$x, y = dens_xstar$y, type = "Algorithm 3")

## generate MLEs corresponding to Algorithm 1, where m_1 and var_1
## are the mean and variance on the logit scale
m_1 <- log(3/108) - log(1 - 3/108)
var_1 <- 1/sqrt(109*(3/108)*(1-3/108))
x_fit1 <- seq(floor(qnorm(0.001,m_1,sqrt(var_1))),
              ceiling(qnorm(1- 0.001,m_1,sqrt(var_1))), by = 0.01)
fit1 <- data.frame(x = x_fit1, y = dnorm(x_fit1,m_1,sqrt(var_1)), type = "Algorithm 1")

## make the plot for the female provider group data (top left corner)
plotC3a <-
  ggplot(data=df_MLEs, aes(df_MLEs$x)) + theme_bw() +
  geom_histogram(aes(y=..density..), position="identity", breaks = bins2, alpha=0.2, col = "grey25", size =1) +
  coord_cartesian(xlim = c(-9.65, -1)) +
  labs(title=bquote('Sampling Distribution of the MLE for logit('*italic(Z)[11]*')')) +
  labs(x=bquote('logit('*italic(Z)[11]*')'), y="Density\n") + 
  theme(plot.title = element_text(hjust = 0.5,size=17,face="bold")) + 
  geom_area(data = fit1, aes(x=x, y=y), fill="white", col="#077DAA", alpha=0, size = 1) +
  geom_area(data = d_xstar, aes(x=x, y=y), fill="white", col="#f49b3f", alpha=0, size = 1) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=16)) +
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))

d_alg <- rbind(d_xstar, fit1)
plotC3b <- ggplot(data=d_alg, aes(x=x, y= y)) + theme_bw() +
  geom_line(aes(y = y, color=as.factor(type)), size = 1) +
  theme(legend.position="bottom") +
  scale_color_manual(name = "Approximation  ", 
                     c("Algorithm 1   ", 
                       "Algorithm 3"),
                     values = c("#077DAA", "#f49b3f")) +
  theme(legend.text=element_text(size=16),
        legend.title=element_text(size=16)) 

mylegend <- get_legend(plotC3b)

## combine two plots
fig_final <- plot_grid(plotC3a, mylegend, ncol = 1, rel_heights = c(2, .3))

# output as .pdf file for the article
pdf(file = "FigureC3OC.pdf",   # The directory you want to save the file in
    width = 5.75, # The width of the plot in inches
    height = 4.25) # The height of the plot in inches

fig_final

dev.off()

## create Figure C4 to discuss the linear approximations on the logit scale

## we'll use three points from the green region for illustration
## while considering sample sizes between 90 and 150
delta_L <- -0.5; delta_U <- 4
params = greens
n_vals <- seq(90, 150, 1)

sob_pwr <- sobol(8192, d = 9, randomize = "digital.shift", seed = 1)

## first check the linear approximations using Algorithm 3
for (k in 1:length(n_vals)){
  n_val <- n_vals[k]
  assign(paste0("logit_probs_", n_val), NULL)
  for (i in c(1,2,12)){
    ## use point i from the Sobol' sequence
    assign(paste0("logit_probs_", n_val), c(get(paste0("logit_probs_", n_val)),
                                            targetPowerAlg3(n_val = n_val, q = 1.25,
                                                        params = as.numeric(params[i,]), 
                                                        hyper = cbind(alphas1, betas1, alphas2, betas2),
                                                        delta = c(delta_L, delta_U), u = sob_pwr[i,2:9])))
  }
}

## combine the results from different sample sizes for each point in 
## the Sobol' sequence
create_linear_approx <- function(start, stop, ind){
  temp <- NULL
  for (j in seq(start, stop, 1)){
    ## append value corresponding to the next sample size
    temp <- c(temp, get(paste0("logit_probs_", j))[ind])
  }
  temp
}

alg3_1 <- create_linear_approx(90, 150, 1)
alg3_2 <- create_linear_approx(90, 150, 2)
alg3_3 <- create_linear_approx(90, 150, 3)

## next check the linear approximations using the discrete process
for (k in 1:length(n_vals)){
  n_val <- n_vals[k]
  assign(paste0("logit_probs_", n_val), NULL)
  for (i in c(1,2,12)){
    ## use point i from the Sobol' sequence
    assign(paste0("logit_probs_", n_val), c(get(paste0("logit_probs_", n_val)),
                                            targetPowerDiscrete(n_val = n_val, q = 1.25,
                                                            params = as.numeric(params[i,]), 
                                                            hyper = cbind(alphas1, betas1, alphas2, betas2),
                                                            delta = c(delta_L, delta_U), u = sob_pwr[i,2:9])))
  }
}

disc_1 <- create_linear_approx(90, 150, 1)
disc_2 <- create_linear_approx(90, 150, 2)
disc_3 <- create_linear_approx(90, 150, 3)

df_lin <- data.frame(n = rep(90:150), 
                     y = c(disc_1, disc_2, disc_3, alg3_1, alg3_2, alg3_3),
                     type = rep(1:6, each = 61))

plotC4 <-
  ggplot(data=df_lin, aes(x=n, y = y)) + theme_bw() +
  geom_line(aes(y = y, color=as.factor(type), linetype=as.factor(type)), size = 1) +
  scale_color_manual(name = "", 
                     values = c("grey40", "#077DAA", "#f49b3f",
                                darken("grey40", 0.5), darken("#077DAA", 0.25), darken("#f49b3f", 0.25))) +
  scale_linetype_manual(name = "", 
                     values = rep(c(1,2), each = 3)) +
  labs(title=bquote('Linear Approximations on the Logit Scale')) +
  labs(x=bquote(italic('n')), y=bquote('logit( '*italic(p)[italic(n)*','*italic(q)*','*bolditalic(u)[r]]^{delta['U'] - delta['L']}*' )')) + 
  theme(plot.title = element_text(hjust = 0.5,size=17,face="bold")) + 
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=16)) +
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) +
  theme(axis.title.x = element_text(margin = unit(c(0, 10, 0, 0), "mm"))) +
  theme(legend.position="none") + ylim(0,11.25) +
  scale_x_continuous(breaks=seq(90,150,20))


pdf(file = "FigureC4OC.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4.25) # The height of the plot in inches

plotC4

dev.off()