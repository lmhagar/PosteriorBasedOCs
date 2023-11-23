## code to prepare Figure 2

## BEGIN SETUP ##

## load necessary packages
library(ggplot2)
library(cowplot)
library(ggpubr)
library(rjags)
library(coda)

## extract design prior parameters for groups 1 and 2
## from the Shiny app
alphas1 <- round(read.csv("beta_priors_new.csv")$alpha,2)
betas1 <- round(read.csv("beta_priors_new.csv")$beta,2)

alphas2 <- round(read.csv("beta_priors_existing.csv")$alpha,2)
betas2 <- round(read.csv("beta_priors_existing.csv")$beta,2)

## generate observations from the design priors
set.seed(3)
z11 = rbeta(1000000, alphas1[1], betas1[1]); z12 = rbeta(1000000, alphas1[2], betas1[2])
z13 = rbeta(1000000, alphas1[3], betas1[3]); z14 = rbeta(1000000, alphas1[4], betas1[4])
z21 = rbeta(1000000, alphas2[1], betas2[1]); z22 = rbeta(1000000, alphas2[2], betas2[2])
z23 = rbeta(1000000, alphas2[3], betas2[3]); z24 = rbeta(1000000, alphas2[4], betas2[4])

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

## get draws from the induced design priors on theta1 and theta2
theta1 <- p11 + 2*p12 + 3*p13 + 4*p14 + 5*(1 - p11 - p12 - p13 - p14)
theta2 <- p21 + 2*p22 + 3*p23 + 4*p24 + 5*(1 - p21 - p22 - p23 - p24)

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
  scale_x_continuous(breaks = seq(3.6, 4.8, by = 0.3), limits = c(3.5, 4.9)) + ylim(0, 3.05) 

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
  scale_x_continuous(breaks = seq(3.6, 4.8, by = 0.3), limits = c(3.5, 4.9)) + ylim(0, 3.05)

## create density plot for induced prior on theta (with red and green shading)
d3 <- density(theta1 - theta2)
d_frame3 <- data.frame(x = d3$x, y = d3$y)

plot2c <- ggplot(data=d_frame3, aes(x=x)) + theme_bw() +
  geom_polygon(aes(y=y), col="gray10", fill="snow1", size=0.75, alpha=0.9) +
  labs(x=bquote(bold('\n Prior for'~ theta)), y=" ") + 
  theme(plot.title = element_text(margin = unit(c(0, 0, 5, 0), "mm"),
                                  hjust = 0.5,size=18,face="bold", colour = "#FFFFFF")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) + 
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) +
  scale_x_continuous(breaks = seq(-0.8, 0.4, by = 0.3), limits = c(-0.95, 0.45))

## create additional data frames to add the green and red
## shading on the plot
G.L <- -0.3; G.U <- -0.1
CPM_subset = subset(d_frame3, x > G.L & x < G.U)
x_new = 
  c(max(G.L,min(CPM_subset$x)), unlist(CPM_subset$x), 
    min(G.U,max(CPM_subset$x)))
y_new = c(0.005,unlist(CPM_subset$y),0.005)
CPM_frame = data.frame(x = x_new,y = y_new)

R.L <- -0.55; R.U <- -0.5
CPM_subset2 = subset(d_frame3, x > R.L & x < R.U)
x_new2 = 
  c(max(R.L,min(CPM_subset2$x)), unlist(CPM_subset2$x), 
    min(R.U,max(CPM_subset2$x)))
y_new2 = c(0.005,unlist(CPM_subset2$y),0.005)
CPM_frame2 = data.frame(x = x_new2,y = y_new2)

## make plot for one-sided superiority test (top right corner)
plot2d <- plot2c +
  geom_polygon(data = CPM_frame, aes(x=x, y=y), fill="seagreen", col="gray10", size = 0.75, alpha = 0.75) +
  geom_polygon(data = CPM_frame2, aes(x=x, y=y), fill="firebrick", col="gray10", size = 0.75, alpha = 0.75)

fig1 <- plot_grid(plot2a + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")), 
                       plot2b + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")),
                        plot2d + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")),
                       rel_widths = c(1, 1, 1), ncol = 3)

## output as .pdf file for the article
pdf(file = "Figure2OC.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches (12.41)
    height = 3.65) # The height of the plot in inches (10.7)

fig1

dev.off()