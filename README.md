# BaySemiCompeting
A Bayesian Nonparametric Approach for Evaluating the Causal Effect of Treatment in Randomized Trials with Semi-Competing Risks

This repository contains a simulated dataset that has similar design as the real data in the paper "A Bayesian Nonparametric Approach for Evaluating the Causal Effect of Treatment in Randomized Trials with Semi-Competing Risks" and the code implementing the proposed BNP method and the algorithm to compute the causal estimand, as well as the code to reproduce the results in the supplemtary material.  

1. "data.Rdata" contains the simulated dataset. 

2. The function “main.R” includes the functions to implement the proposed BNP model, compute the marginal survival distributions under treatments, calculate the causal estimand, and reproduce the results in the supplemtary material. 

3. "mcmc.R" is the function to obtain the MCMC results from the proposed BNP model. 

4. “Marginal_Surv.R” is the function to compute the mariginal survival curves. It needs te posterior samples from the MCMC results. 

5. "Estimand.R" is the function to calculate the proposed causal estimand. It needs the posterior samples from the MCMC results. 

6. The function “update_pos.R” includes all the supporting functions for "mcmc.R", “Marginal_Surv.R”, and "Estimand.R". 

7. "saved_mcmc.RData" is the saved MCMC results. Since the "mcmc.R" needs time to run, we saved the MCMC results so that the readers can run "“Marginal_Surv.R" and "Estimand.R" without running "mcmc.R". 

8. "Figure1.RData" saves the results from computing the mariginal survival probabilities so that the readers can easily reproduce Figure S1 in the supplemtary material.

9. "hu_rho1.RData" and "hu_rho1.RData" save the results from computing the proposed causal estimand for rho=0.2 and rho=0.8 respectively, so that the readers can easily reproduce Figure S2 in the supplemtary material. 

10. "supp.pdf" is the supplemtary material for the paper "A Bayesian Nonparametric Approach for Evaluating the Causal Effect of Treatment in Randomized Trials with Semi-Competing Risks". 
