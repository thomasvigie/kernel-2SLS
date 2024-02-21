# This code contains the parameters of the simulation study ---------------------

# Sample size
N = [100, 1000]
# Parameter of interest
beta_0 = 1
# Number of replications
nsim = 2
# Different first stage functional forms
fs_fun = ["log", "linear", "A&L", "binary_x", "normal density"]
#fs_fun = ["A&L"]

# Bandwidths
bandrange = [0.01, 0.1, 1, 10, 100]
 mu = [0.0, 0]
# Error terms Covariances
sigma_low = [1.0, 0.1, 0.1, 1.0]
sigma_mid = [1.0, 0.5, 0.5, 1.0]
sigma_high = [1.0, 0.9, 0.9, 1.0]

sigma_low = reshape(sigma_low, 2, 2)
sigma_mid = reshape(sigma_mid, 2, 2)
sigma_high = reshape(sigma_high, 2, 2)

SIGMA = [sigma_low, sigma_mid, sigma_high]
