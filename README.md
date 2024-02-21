This repository contains the scripts used for the simulation study in Vigié (2022): "A kernel-based first stage in linear instrumental variable models" (pdf included in the repository).
Description of the scripts (more about each function inside the scripts):
- auxiliary functions.R: contains small functions that are called in the simulation study
- WMD.R: Implements the Weighted Minimum Distance (WMD) estimator of Antoine and Lavergne (2014) 
- kernel estimators.R: contains "lcll", the functions that computes the kernel regression estimates of f(x) in y = f(x) + u. Contains the same ridging process as the npreg function of the np package of Hayfield and Racine, if numerical issues arise
- kernel 2SLS.R: contains the function that estimates the coefficients beta in y = [X, W] * beta + u, where X is an endogenous variable, W are exogenous variables, and Z are the instruments used for the two stage least squares estimates.
                 The first stage is estimated via kernel regression, using the npreg function of the np package of Hayfield and Racine
- DGP.R: contains the data generating process for the simulation study. Multiple options are available, see details in the script
- Inference.R: contains various inference-related functions (needs to be polished)
- bandwidth choice criteria.R: contains various objective functions for the choice of the optimal bandwidth in linear instrumental variable models where the first stage regressions uses kernel estimators (needs to be polished)
- leave one out.R: contains various functions for the computation of cross validation criterion functions for the optimal bandwidth (needs to be polished)
