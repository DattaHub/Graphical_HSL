Functions
BayesGLasso_Columnwise.m
BayesGLassoGDP.m
glasso_adaptive_cv.m
glasso_cv.m
glasso_FTH.m
glasso_SCAD_cv.m
rand_ig.m
are necessary for running BGLASSO and GSCAD for comparisons. 

1. The R file 'GHS_generate_data.R' generates data for the simulations. 
This file generates one data set from all precision matrix structures when p=100 and 
n=120, and saves them to a .RData file (for GLASSO estimates using R) and  .csv files
 (for HSL_ECM, HSL_MCMC, GHS, BGLASSO, and GSCAD estimates using Matlab). Code for other dimensions
 can be found in the comments in this file.

2. The Matlab file 'GHS_sim_GHS.m' reads .csv data files and estimates the precision 
matrix by graphical horseshoe. Matlab function 'GHS.m' is called for estimation. 
This file prints the mean and standard deviation of Stein's loss, Frobenius norm, 
true positive rate (TPR), false positive rate (FPR), MCC and time of the estimates. 
For demonstration, 500 burn-in samples and 1000 MCMC samples are used in estimation.

3. The Matlab file 'HSL_MCMC_sim_HSL.m' reads .csv data files and estimates the precision 
matrix by graphical horseshoe like penalty (from MC samples). 
Matlab function 'HSL_MCMC.m' is called for estimation. 
This file prints the mean and standard deviation of Stein's loss, Frobenius norm, 
true positive rate (TPR), false positive rate (FPR), MCC and time of the estimates. 
For demonstration, 500 burn-in samples and 1000 MCMC samples are used in estimation.

4. The Matlab file 'HSL_ECM_50.m' reads .csv data files and estimates the precision 
matrix by graphical horseshoe like penalty (using Expectation Conditional Maximization). 
Matlab function 'Multi_start_point_Fixed_b_EM_HS_like.m' is called for estimation. 
This file prints the mean and standard deviation of Stein's loss, Frobenius norm, 
true positive rate (TPR), false positive rate (FPR), MCC and time of the estimates. 

5. The R file 'GHS_sim_GLASSO.r' reads .RData data file and estimates the precision
 matrix by graphical lasso, prints the mean and standard deviation of Stein's loss, Frobenius norm, 
true positive rate (TPR), false positive rate (FPR), MCC and time of the estimates. 

6. The Matlab file 'GHS_sim_GLASSO.m' reads .csv data files and estimates the precision 
matrix by Bayesian graphical lasso. Matlab function 'BayesGLasso_Columnwise.m' is called for estimation. 
This file prints the mean and standard deviation of Stein's loss, Frobenius norm, 
true positive rate (TPR), false positive rate (FPR), MCC and time of the estimates. 

7. The Matlab file 'GHS_sim_GSCAD.m' reads .csv data files and estimates the precision 
matrix by graphical SCAD. Matlab function 'glasso_SCAD_cv.m' is called for estimation. 
This file prints the mean and standard deviation of Stein's loss, Frobenius norm, 
true positive rate (TPR), false positive rate (FPR), MCC and time of the estimates. 

8. The R file scale_factor_computation.R computes the global scale (shrinkage)
parameter for given values of n, p. If you are trying out for dimensions which are not in the paper,
please run this R code and make corresponding changes in 'HSL_MCMC.m' and 
'Multi_start_point_Fixed_b_EM_HS_like.m'. 