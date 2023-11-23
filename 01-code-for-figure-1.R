## code to prepare figure for Section 2

## BEGIN SETUP ##

## load necessary packages
library(ggplot2)
library(cowplot)
library(ggpubr)
library(rjags)
library(coda)

## create data frame for observed Likert data
df.maize <- data.frame(Variety = c(rep("One", 5), rep("Two", 5)),
                       Prop = c(c(3,7,14,28,56)/108, c(2,2,14,43,76)/137),
                       Score = c(seq(1,5,1), seq(1,5,1)))

## create barplot
dat_plot <- ggplot(df.maize,aes(Score,Prop,fill=Variety)) + theme_bw() +
  geom_bar(stat="identity",position='dodge') + 
  theme(legend.position="top") +
  theme(legend.text=element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  scale_fill_manual(name = " ", 
                     labels = c(expression("MH44 A: "*hat(theta)[1]*" = 4.18  "), 
                                expression("MH43 A: "*hat(theta)[2]*" = 4.38")),
                     values = c(
                       'One' = '#077DAA',
                       'Two' = '#FED758')) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16)) +
  labs(y = "Proportion", xlab = "Likert Score") + ylim(0, 0.6)

## get the posterior of theta_1 - theta_2
## define priors
scl <- 0.8
alphas1 <- scl*rep(1,4)
betas1 <- scl*seq(4,1)

alphas2 <- scl*rep(1,4)
betas2 <- scl*seq(4,1)

nchains <- 1
burnin <- 1000
nthin <- 2
ndraws <- 100000*nthin

## use JAGS to obtain the posterior approximations with the model in 
## the JAGS_ordinal.txt file
set.seed(1)
dat1 <- c(3,7,14,28,56)
model1.fit <- jags.model(file="JAGS_ordinal.txt",
                         data=list(N=sum(dat1), x=as.numeric(dat1), a=alphas1, b=betas1), 
                         n.chains = nchains)
## burnin
update(model1.fit, burnin)

m1.samp <- coda.samples(model1.fit, c("p"), n.iter=ndraws, thin=nthin)[[1]]

## get posterior draws for ordinal mean in group 1
theta.1 <- m1.samp[,1] + 2*m1.samp[,2] + 3*m1.samp[,3] + 4*m1.samp[,4] + 5*m1.samp[,5]

## repeat for model 2
set.seed(2)
dat2 <- c(2,2,14,43,76)
model2.fit <- jags.model(file="JAGS_ordinal.txt",
                         data=list(N=sum(dat2), x=as.numeric(dat2), a=alphas2, b=betas2), 
                         n.chains = nchains)
update(model2.fit, burnin)

m2.samp <- coda.samples(model2.fit, c("p"), n.iter=ndraws, thin=nthin)[[1]]

theta.2 <- m2.samp[,1] + 2*m2.samp[,2] + 3*m2.samp[,3] + 4*m2.samp[,4] + 5*m2.samp[,5]

## create data frame to plot density function of posterior for theta = theta1 - theta2
d <- density(theta.1 - theta.2)
d_frame <- data.frame(x = d$x, y = d$y)

## make plot with no grey shading first (additional layers to be added later)
plot1e <- ggplot(data=d_frame, aes(x=d_frame$x)) + theme_bw() +
  geom_polygon(aes(y=d_frame$y), col="gray10", fill="snow1", size=0.75, alpha=0.9) +
  labs(title = bquote(bold("Posterior pdf of "~ theta[1] ~ '-' ~ theta[2]))) +
  labs(x=bquote(bold('\n Posterior pdf of'~ theta[1] ~ '-' ~ theta[2])), y=" ") + 
  theme(plot.title = element_text(margin = unit(c(0, 0, 5, 0), "mm"),
                                  hjust = 0.5,size=18,face="bold", colour = "#FFFFFF")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) + 
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) + 
  scale_x_continuous(breaks = seq(-0.6, 0.3, by = 0.3))

## compute posterior probability numerically for the plot title
delta.L <- -0.5
delta.U <- Inf
theta <- theta.1 - theta.2
CPM_4 <- round(mean(ifelse(theta < delta.U, theta > delta.L,0)),4)
CPM_title <- bquote(bold(italic(Pr)*'( '* theta >-0.5*' | '*italic(data)*') =' ~ .(CPM_4)))

x_title <- bquote(bold('\n Posterior pdf of'~ theta[1] ~ '-' ~ theta[2]))

## define second data frame for the shaded portion of the posterior
CPM_subset = subset(d_frame, x > delta.L & x < delta.U)
x_new = 
  c(max(delta.L,min(CPM_subset$x)), unlist(CPM_subset$x), 
    min(delta.U,max(CPM_subset$x)))
y_new = c(0,unlist(CPM_subset$y),0)
CPM_frame = data.frame(x = x_new,y = y_new)

## make plot for one-sided noninferiority test
plot1b <- 
  plot1e + labs(title = NULL) + labs(x = NULL) +
  labs(title=CPM_title) +
  labs(x= x_title) +
  geom_polygon(data = CPM_frame, aes(x=x, y=y), fill="gray", col="gray10", size = 0.75, alpha = 0.75) +
  theme(plot.title = element_text(hjust = 0.5,size=18,face="bold", colour = "#000000"))

fig1 <- plot_grid(dat_plot + theme(plot.margin=unit(c(0.25,0.25,0.45,0.25),"cm")), 
                       plot1b + theme(plot.margin=unit(c(0.45,0.25,0.25,0.25),"cm")),
                       rel_widths = c(1.45, 1.15))

## output as .pdf file for the article
pdf(file = "Figure1OC.pdf",   # The directory you want to save the file in
    width = 9.7895, # The width of the plot in inches (12.41)
    height = 4.15) # The height of the plot in inches (10.7)

fig1

dev.off()