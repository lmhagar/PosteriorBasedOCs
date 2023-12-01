## code to determine design values and informative priors
## this procedure involves processing data from the ENIGH 2018 survey

## BEGIN SETUP ##

## Process to obtain data file from ENIGH website
## Step 1: Go to this website: https://www.inegi.org.mx/programas/enigh/nc/2018/#Datos_abiertos
## Step 2: Click the "Download CSV" button to download a .zip file
## Step 3: Open .zip file
## Step 4: Open the "conjunto_de_datos_concentradohogar_enigh_2018_ns" subfolder
## Step 5: Open the "conjunto_de_datos" subfolder
## Step 6: Extract the lone .csv file: "conjunto_de_datos_concentradohogar_enigh_2018_ns.csv" and move it to the working directory

require(mvtnorm)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(rjags)
library(coda)
## read data set
enigh <- read.csv("conjunto_de_datos_concentradohogar_enigh_2018_ns.csv")

## create a numeric character for region
enigh$region <- as.numeric(substr(as.character(enigh$ubica_geo),1,nchar(as.character(enigh$ubica_geo))-3))

## keep only the households from Aguascalientes (Region = 1):
aguas <- enigh[enigh$region %in% c(1),]

## keep only columns of interest
aguas <- aguas[c("region","factor","est_socio","educa_jefe", "tot_integ", "alimentos", "sexo_jefe")]

## convert column titles to English
names(aguas)[4] <- "education"
names(aguas)[5] <- "total_people"
names(aguas)[6] <- "total_food"
names(aguas)[7] <- "sex"

## remove households with 0 quarterly expenditure on food
aguas_full <- aguas ## save all data for later
aguas <- aguas[aguas$total_food > 0,]

## keep only individuals with estimated socioeconomic class 2
aguas <- aguas[aguas$est_socio ==2 ,]

## create simplified weighting factor 
aguas$factor2 <- round(aguas$factor/75)

## repeat the observations according to new weighting factor
aguas_long <- aguas[1,]
for (i in 1:nrow(aguas)){
  if (i %% 100 == 0){
    print(i)
  }
  for (j in 1:aguas$factor2[i]){
    aguas_long <- rbind(aguas_long, aguas[i,])
  }
}
aguas_long <- aguas_long[-1,]

## calculate food expense per person in thousands of pesos (adjusted for 2 percent inflation)
aguas_long$food <- 0.001*aguas_long$total_food/aguas_long$total_people*(1.02)^2

## split based on sex of main household provider
aguas_F <- aguas_long[aguas_long$sex == 2,]
aguas_M <- aguas_long[aguas_long$sex ==1,]

## remove households with more than 20000 pesos per person per month;
## this is three households with female providers, and 2 with male providers
aguas_F <- subset(aguas_F, aguas_F$food <= 20)
aguas_M <- subset(aguas_M, aguas_M$food <= 20)

## save food expenditure data for both groups
write.csv(aguas_F[c("food")], "aguas2018_food_F.csv", row.names = FALSE)
write.csv(aguas_M[c("food")], "aguas2018_food_M.csv", row.names = FALSE)

## we now obtain design values for alpha1, beta1, alpha2, and beta2
## first load the .csv files just created
## save the food data in vector form
y1 <- read.csv("aguas2018_food_F.csv")$food 
y2 <- read.csv("aguas2018_food_M.csv")$food

## find design priors for the Weibull distribution
set.seed(3)
burnin <- 1000
nchains <- 1
nthin <- 2
ndraws <- nthin*100000
n1 <- length(y1)
mu0 <- 2; tau0 <- 1; kappa0 <- 2; lambda0 <- 1
model1.fit <- jags.model(file="JAGS_weibull.txt",
                         data=list(n=n1, y = y1, 
                                   tau0 = tau0, mu0 = mu0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model1.fit, burnin)
model1.samples <- coda.samples(model1.fit, c("nu", "l"), n.iter=ndraws, thin=nthin)

nu.1 <- unlist(model1.samples[,2])
lambda.1 <- unlist(model1.samples[,1])

set.seed(4)
n2 <- length(y2)
model2.fit <- jags.model(file="JAGS_weibull.txt",
                         data=list(n=n2, y = y2, 
                                   tau0 = tau0, mu0 = mu0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model2.fit, burnin)
model2.samples <- coda.samples(model2.fit, c("nu", "l"), n.iter=ndraws, thin=nthin)

nu.2 <- unlist(model2.samples[,2])
lambda.2 <- unlist(model2.samples[,1])

## find posterior means for design values
mean(nu.1) ## should be 1.41
mean(lambda.1) ## should be 3.39
mean(nu.2) ## should be 1.49
mean(lambda.2) ## should be 3.42

## save posterior draws to .csv files for reference
write.csv(nu.1, "nu1s_2018.csv", row.names = FALSE)
write.csv(lambda.1, "lambda1s_2018.csv", row.names = FALSE)
write.csv(nu.2, "nu2s_2018.csv", row.names = FALSE)
write.csv(lambda.2, "lambda2s_2018.csv", row.names = FALSE)

## now we inflate the variances of these approximately gamma
## marginal posteriors to find informative priors for the numerical study

## define factor by which to inflate the variance
var_inf <- 30

## consider nu1
## find the posterior mode and variance
a1dens <- density(nu.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(nu.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star

## should be GAMMA(38.07,  26.23)
informs <- c(a1_star, b1_star)

## consider lambda1
## find the posterior mode and variance
a1dens <- density(lambda.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(lambda.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star

## should be GAMMA(34.92,  10.02)
informs <- rbind(informs, c(a1_star, b1_star))

var_inf <- 100
## consider nu2
## find the posterior mode and variance
a1dens <- density(nu.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(nu.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star

## should be GAMMA(38.35, 25.09)
informs <- rbind(informs, c(a1_star, b1_star))

## consider lambda2
## find the posterior mode and variance
a1dens <- density(lambda.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(lambda.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star

## should be GAMMA(37.51, 10.70)
informs <- rbind(informs, c(a1_star, b1_star))

## output informative weibull distribution hyperparameters to .csv file
write.csv(informs, "informs_weibull.csv", row.names = FALSE)

## obtain parameter values from the priors
set.seed(10)
## generate observations from the design priors
k1s <- rgamma(5000000, informs[1,1], informs[1,2])
l1s <- rgamma(5000000, informs[2,1], informs[2,2])
k2s <- rgamma(5000000, informs[3,1], informs[3,2])
l2s <- rgamma(5000000, informs[4,1], informs[4,2])

## obtain theta values for each draw to segment the priors
theta1 <- l1s*log(1/(1 - threshold))^(1/k1s)
theta2 <- l2s*log(1/(1 - threshold))^(1/k2s)

## do the sampling on the log-scale
## define the green and red regions
delta_L <- 1/1.2; delta_U <- 1.2
green_high <- 1.05; green_low <- 1/green_high
red_low <- delta_U; red_high <- red_low + 0.025

## approximate the induced prior density function on theta for sampling-resampling
diff_dens <- density(log(theta1) - log(theta2))
diff_pdf <- approxfun(x = diff_dens$x,
                      y = diff_dens$y)

## conduct sampling-resampling to uniformly sample from the green region
diff_num <- ifelse(theta1/theta2 <= green_high, theta1/theta2 >= green_low, 0)
diff_denom <- diff_pdf(log(theta1) - log(theta2))

weights <- diff_num/diff_denom
weights <- weights/sum(weights)

set.seed(1)
## resample to get 8192 parameter values from green region
inds_region <- sample(1:5000000, size = 8192, replace = TRUE, prob = weights)

## output parameter values to a .csv file
write.csv(cbind(l1s, k1s, l2s, k2s)[inds_region,], "weibull_green.csv", row.names = FALSE)
write.csv(theta1[inds_region], "wei_theta1_green.csv", row.names = FALSE)
write.csv(theta2[inds_region], "wei_theta2_green.csv", row.names = FALSE)

## conduct sampling-resampling to uniformly sample from the red region
diff_num <- ifelse((theta1/theta2 <= red_high & theta1/theta2 >= red_low) | 
                   (theta1/theta2 <= 1/red_low & theta1/theta2 >= 1/red_high), 1, 0)
diff_denom <- diff_pdf(log(theta1) - log(theta2))

weights <- diff_num/diff_denom
weights <- weights/sum(weights)

## resample to get 8192 parameter values from red region
set.seed(2)
inds_region2 <- sample(1:5000000, size = 8192, replace = TRUE, prob = weights)

## output parameter values to a .csv file
write.csv(cbind(l1s, k1s, l2s, k2s)[inds_region2,], "weibull_red.csv", row.names = FALSE)
write.csv(theta1[inds_region2], "wei_theta1_red.csv", row.names = FALSE)
write.csv(theta2[inds_region2], "wei_theta2_red.csv", row.names = FALSE)

## create the plots for figure D.2
## approximate the posterior using nonparametric density estimation
d1 <- density(theta1)
d_frame1 <- data.frame(x = d1$x, y = d1$y)

## create density plot for induced prior on theta1
plot2a <- ggplot(data=d_frame1, aes(x=x)) + theme_bw() +
  geom_polygon(aes(y=y), col="gray10", fill="snow1", size=0.75, alpha=0.9) +
  labs(x=bquote(bold('\n Prior for'~ theta[1])), y=" ") + 
  theme(plot.title = element_text(margin = unit(c(0, 0, 5, 0), "mm"),
                                  hjust = 0.5,size=18,face="bold", colour = "#FFFFFF")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) + 
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 13)) + ylim(0, 0.4) 

## create density plot for induced prior on theta2
d2 <- density(theta2)
d_frame2 <- data.frame(x = d2$x, y = d2$y)

plot2b <- ggplot(data=d_frame2, aes(x=x)) + theme_bw() +
  geom_polygon(aes(y=y), col="gray10", fill="snow1", size=0.75, alpha=0.9) +
  labs(x=bquote(bold('\n Prior for'~ theta[2])), y=" ") + 
  theme(plot.title = element_text(margin = unit(c(0, 0, 5, 0), "mm"),
                                  hjust = 0.5,size=18,face="bold", colour = "#FFFFFF")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) + 
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) + ylim(0, 0.4) 

## create density plot for induced prior on theta (with red and green shading)
d3 <- density(log(theta1) - log(theta2))
d_frame3 <- data.frame(x = d3$x, y = d3$y)

plot2c <- ggplot(data=d_frame3, aes(x=x)) + theme_bw() +
  geom_polygon(aes(y=y), col="gray10", fill="snow1", size=0.75, alpha=0.9) +
  labs(x=bquote(bold('\n Prior for log('*theta*')')), y=" ") + 
  theme(plot.title = element_text(margin = unit(c(0, 0, 5, 0), "mm"),
                                  hjust = 0.5,size=18,face="bold", colour = "#FFFFFF")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) + 
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.5), limits = c(-1.1, 1.1))

## create additional data frames to add the green and red
## shading on the plot
G.L <- -1*log(1.05); G.U <- log(1.05)
CPM_subset = subset(d_frame3, x > G.L & x < G.U)
x_new = 
  c(max(G.L,min(CPM_subset$x)), unlist(CPM_subset$x), 
    min(G.U,max(CPM_subset$x)))
y_new = c(0.005,unlist(CPM_subset$y),0.005)
CPM_frame = data.frame(x = x_new,y = y_new)

R.L <- -1*log(1.225); R.U <- -1*log(1.2)
CPM_subset2 = subset(d_frame3, x > R.L & x < R.U)
x_new2 = 
  c(max(R.L,min(CPM_subset2$x)), unlist(CPM_subset2$x), 
    min(R.U,max(CPM_subset2$x)))
y_new2 = c(0.005,unlist(CPM_subset2$y),0.005)
CPM_frame2 = data.frame(x = x_new2,y = y_new2)

R.U <- log(1.225); R.L <- log(1.2)
CPM_subset3 = subset(d_frame3, x > R.L & x < R.U)
x_new3 = 
  c(max(R.L,min(CPM_subset3$x)), unlist(CPM_subset3$x), 
    min(R.U,max(CPM_subset3$x)))
y_new3 = c(0.005,unlist(CPM_subset3$y),0.005)
CPM_frame3 = data.frame(x = x_new3,y = y_new3)

## make plot for one-sided superiority test (top right corner)
plot2d <- plot2c +
  geom_polygon(data = CPM_frame, aes(x=x, y=y), fill="seagreen", col="gray10", size = 0, alpha = 0.75) +
  geom_polygon(data = CPM_frame2, aes(x=x, y=y), fill="firebrick", col="gray10", size = 0, alpha = 0.75) +
  geom_polygon(data = CPM_frame3, aes(x=x, y=y), fill="firebrick", col="gray10", size = 0, alpha = 0.75)

fig1 <- plot_grid(plot2a + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")), 
                  plot2b + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")),
                  plot2d + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")),
                  rel_widths = c(1, 1, 1), ncol = 3)

## output as .pdf file for the article
pdf(file = "FigureD2OC.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches (12.41)
    height = 3.65) # The height of the plot in inches (10.7)

fig1

dev.off()