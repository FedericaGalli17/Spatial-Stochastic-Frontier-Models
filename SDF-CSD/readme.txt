Instructions for using the code

MC.m is the main file with the code of the simulations for 5 time periods.
Loglik2_T5.m is the corresponding log-likelihood function.
All the other files are utilities called by the code for the computation of the log-determinant, for computing the standard errors, etcÉ (references and usage are described therein).

To use the code for empirical applications, load the data in the MC file modifying the code including new rows for the number of X and Z variables needed as well as for the number of time periods. At the moment, the MC file just considers one X variable and one Z variable. Then, add them also in the equations in the log-likelihood function. 


