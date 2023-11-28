This folder contains the codes used the estimate the model developed in the paper "Accounting for Unobserved Individual Heterogeneity in Spatial Stochastic Frontier Models: the case of Italian Innovative Start-ups". Specifically,
-"MC.m" allows to test the model on simulated data and "loglik.m" is the related log-likelihood function.
-"codeapplication.m" allows to estimate the empirical model is section 5 (see below for a description of the data) and "loglikstartupWP" is the related log-likelihood function. 
-"scores.m" is used for the computation of the efficiency scores.
-all others codes have been borrowed for other authors and allow to obtain robust standard errors, work with spatial weight matrices, estimate the log-determinant, etc. The detailed references and usage can be found therein.


----------------------------------------
Information on the data fileData are structured as follows:
-each observation is reported by row;
-columns 1 to 5 contain data on value added (logarithm of real value added) from 2016 to 2020; 
-columns 6 to 10 contain data on labour (logarithm of total wages) from 2016 to 2020;
-columns 11 to 15 contain data on capital (logarithm of gross fixed capital) from 2016 to 2020;
-columns 16 to 20 contain data on size (logarithm of number of employees) from 2016 to 2020;
-columns 21 to 25 contain data on patents (logarithm of patent rights) from 2016 to 2020;-columns 26 to 30 contain data on financial leverage (logarithm of long debts over net equity) from 2016 to 2020;
-columns 31 to 35 contain data on financial flexibility (logarithm of liquid assets over total assets) from 2016 to 2020;
-columns 36 to 40 contain data on intellectual capital (logarithm of intangible capital over total capital) from 2016 to 2020;
-columns 41 to 45 contain data on R&D (logarithm of R&D expenditure) from 2016 to 2020;

Informations are retrieved from the AIDA Bureau Van Dijk database (at this link https://login.bvdinfo.com/R0/AidaNeo). Data refer to all Italian innovative start-ups in the time period 2016-2020. Information with zero or negative values for the value added have been deleted. 

The spatial weight matrix used for the main results is build as an inverse distance matrix truncated at 50km. The matrix imported in Matlab refers to the whole panel (NxN) and is not normalised since normalisation is performed in the code after rows and columns deletion for missing observations in each panel. Starting from the complete matrix based on all observations, for each panel, rows and columns referring to missing observations are removed.





 