## code to reproduce Figure 1 in the main text
## this file should be run after "01-food-data-2020.R"

## BEGIN SETUP ##

## load necessary packages
require(ggplot2)
require(cowplot)
require(rjags)
require(coda)

## load the .csv files created using "01-food-data-2020.R" 
## save the food data in vector form
y1 <- read.csv("aguas_food_F.csv")$food 
y2 <- read.csv("aguas_food_M.csv")$food

## save the data as data frame (for plotting later)
foodF <- read.csv("aguas_food_F.csv")
foodM <- read.csv("aguas_food_M.csv")

## use JAGS to obtain posterior draws via MCMC
set.seed(1)
burnin <- 1000
nchains <- 1
nthin <- 2
ndraws <- nthin*100000
n1 <- length(y1)
tau1 <- 2; mu1 <- 1; kappa1 <- 2; lambda1 <- 1
model1.fit <- jags.model(file="JAGS_weibull.txt",
                         data=list(n=n1, y = y1, 
                                   tau0 = tau1, mu0 = mu1,
                                   kappa0 = kappa1, lambda0 = lambda1), 
                         n.chains = nchains, quiet = TRUE)

update(model1.fit, burnin, progress.bar = "none")
model1.samples <- coda.samples(model1.fit, c("nu", "l"), n.iter=ndraws, thin=nthin)

nu.1 <- unlist(model1.samples[,2])
lambda.1 <- unlist(model1.samples[,1])

write.csv(data.frame(nu.1, lambda.1), "jags_output_group1.csv", row.names = FALSE)

set.seed(2)
n2 <- length(y2)
tau2 <- 2; mu2 <- 1; kappa2 <- 2; lambda2 <- 1
model2.fit <- jags.model(file="JAGS_weibull.txt",
                         data=list(n=n2, y = y2, 
                                   tau0 = tau2, mu0 = mu2,
                                   kappa0 = kappa2, lambda0 = lambda2), 
                         n.chains = nchains, quiet = TRUE)

update(model2.fit, burnin, progress.bar = "none")
model2.samples <- coda.samples(model2.fit, c("nu", "l"), n.iter=ndraws, thin=nthin)

nu.2 <- unlist(model2.samples[,2])
lambda.2 <- unlist(model2.samples[,1])

write.csv(data.frame(nu.2, lambda.2), "jags_output_group2.csv", row.names = FALSE)

## find good histogram breaks based on the data
breaks1 <- hist(y1, plot=FALSE)$breaks
breaks2 <- hist(y2, plot=FALSE)$breaks
binwidth2 <-0.64
bins2 <- seq(0, 19.2, binwidth2)

## create data frames for plotting
x_fit1 <- seq(floor(qweibull(0.001,mean(nu.1),mean(lambda.1))),
              ceiling(qweibull(0.999,mean(nu.1),mean(lambda.1))), by = 0.01)
fit1 <- data.frame(xfit = x_fit1, yfit = dweibull(x_fit1,mean(nu.1),mean(lambda.1)))

x_fit2 <- seq(floor(qweibull(0.001,mean(nu.2),mean(lambda.2))),
              ceiling(qweibull(0.999,mean(nu.2),mean(lambda.2))), by = 0.01)
fit2 <- data.frame(xfit = x_fit2, yfit = dweibull(x_fit2,mean(nu.2),mean(lambda.2)))

hist_lower2 <- 0
hist_upper2 <- 20

group1_name <- "Female Household Provider"
n1 <- length(y1)

use_pctl <- 1
metric_plot <- paste("0.9-Quantile", sep="")
metric_value1 <- round(quantile(y1, 0.9),3)
percentile <- 0.9

## make the plot for the female provider group data (top left corner)
plot1a <-
  ggplot(data=foodF, aes(foodF$food)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))/binwidth2), breaks = bins2,
                 col="#1B86B4", 
                 fill="#077DAA", 
                 alpha = 0, size = 1) + 
  coord_cartesian(xlim = c(hist_lower2, hist_upper2)) +
  labs(title=paste(group1_name, "\n")) +
  labs(x=bquote(atop(' ', atop(textstyle('Food Expenditure per Person (MXN $1000)'),
                               textstyle('n = '*.(n1)*", "*.(metric_plot)*' ('*hat(theta)[1]*  
                                           ') = '*.(metric_value1))))), y="Density\n") + 
  theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")) + 
  theme(plot.subtitle = element_text(hjust = 0.5,size=14)) +
  geom_area(data = fit1, aes(x=xfit, y=yfit), fill="#94D3FF", col="#077DAA", alpha=0.45, size = 1) +
  geom_segment(aes(x = quantile(y1, 0.9), y = 0, 
                   xend=quantile(y1, 0.9), yend = Inf), color="grey16", size=1.5) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) + ylim(0, 0.3)

group2_name <- "Male Household Provider"
n2 <- length(y2)
metric_value2 <- round(quantile(y2, 0.9),3)

## make the plot for the male provider group (bottom left corner)
plot1c <- ggplot(data=foodM, aes(foodM$food)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))/binwidth2), breaks = bins2,
                 col="#9E8203", 
                 fill="#FED758", 
                 alpha = 0, size = 1) + 
  coord_cartesian(xlim = c(hist_lower2, hist_upper2)) +
  labs(title=paste(group2_name,"\n")) +
  labs(x=bquote(atop(' ', atop(textstyle('Food Expenditure per Person (MXN $1000)'),
                               textstyle('n = '*.(n2)*", "*.(metric_plot)*' ('*hat(theta)[2]*  
                                           ') = '*.(metric_value2))))), y="Density\n") + 
  theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_area(data = fit2, aes(x=xfit, y=yfit), fill="#FED758", col="#9E8203", alpha=0.45, size = 1) + 
  geom_segment(aes(x = quantile(y2, 0.9), y = 0, 
                   xend=quantile(y2, 0.9), yend = Inf), color="grey16", size=1.5) + ylim(0, 0.3)

## combine into grid
figD1 <- plot_grid(plot1a + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")), 
                       plot1c + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")),
                       rel_widths = c(1, 1))

## output as .pdf file for the article
pdf(file = "FigureD1OC.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches (12.41)
    height = 4.5) # The height of the plot in inches (10.7)

figD1

dev.off()