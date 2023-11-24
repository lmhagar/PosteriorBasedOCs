# PosteriorBasedOCs
Supplementary code for Scalable Design with Posterior-Based Operating Characteristics manuscript

The code files for this manuscript are divided into several groups.

Group 1: create Figure 1 for illustrative example
- 01-code-for-figure-1: processes the Likert data and approximates posteriors for Figure 1
- JAGS_ordinal.txt: model to approximate the posteriors for the ordinal model using MCMC

Group 2: create Figure 2 to visualize the design priors
- 02-code-for-figure-2: takes parameters for the specified design priors and visualizes the induced priors

Group 3: conduct numerical studies in Sections 4, 5, and Appendix C.2
- 03-functions-for-sec-4-5: defines functions used to conduct the numerical studies with the multinomial example
- 04-multinomial-study-alg-3: conducts the numerical studies and generates the figures for these sections
- priors-group-1.csv: parameters for the beta design priors (group 1)
- priors-group-2.csv: parameters for the beta design priors (group 1)

Group 4: conduct numerical studies in Appendix C.3