Authors: Denis Tverskoi and Grzegorz A. Rempala, 2025. This repository contains matlab program files that reproduce the results of numerical simulations and data fit presented in the research work 
"Model Fit vs. Predictive Reliability: A Case Study of the 1978 Influenza Outbreak".

1. The simulations.m file simulates the posterior distributions of the stochastic SIR model parameters using the Vanilla ABC method. It uses the Gil3.m file to obtain realizations of the stochastic process employing the modified Gillespie algorithm for delayed reactions.
2. The Analysis.m file visualizes the modeling result and computes the goodness-of-fit. It uses the empirical_hdp_4d.m file to compute credible intervals for the parameter estimates.
