Instructions for using the code

MC.m is the main file with the code for the simulations.
loglik.m is the corresponding log-likelihood function.
All the other files are utilities called by the code for the computation of the log-determinant, for computing the standard errors, etcÉ (references and description can be found therein).

To use the code for empirical applications, load the data in the MC file modifying the code including new rows for the number of X and Z variables needed. At the moment, the MC file just considers one X variable and one Z variable. Then, add them also in the equations in the log-likelihood function. 


