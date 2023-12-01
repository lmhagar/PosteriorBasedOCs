# PosteriorBasedOCs
Supplementary code for Scalable Design with Posterior-Based Operating Characteristics manuscript

The code files for this manuscript are divided into five groups.

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
- 05-functions-for-app-c3: defines functions for the supplementary numerical studies with the multinomial example
- 06-multinomial-study-alg-1: conducts the numerical studies and generates the figures for these sections

Group 5: conduct numerical studies in Appendix D
- 07-weibull-data-2020: processes the ENIGH 2020 food data used for illustration
- 08-code-for-figure-D1: reproduces Figure D.1 in the online supplement
- 09-prior-elicitation: processes the ENIGH 2018 food data and uses the relevant posteriors
                        to return design and analysis priors for this example
- 10-functions-for-app-d2: defines functions used to conduct the numerical studies with the Weibull example
- 11-weibull-study: conducts the numerical studies and generates the remaining figures for this section