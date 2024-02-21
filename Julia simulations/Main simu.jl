# This script is a test on some data, to compare if my functions in Julia and R do the same thing
# Pkg.add("Nemo")
# Pkg.add("Polynomials")
# Pkg.add("SharedArrays")
# Pkg.add("Distributed")
# Pkg.add("StatsModels")

using Polynomials
#using StatsModels
using DataFrames
using Distributions
#using KernelEstimator
#using BlackBoxOptim
using LinearAlgebra
using Optim
using CSV
#using NLsolve
using SharedArrays      # Used to create the matrix the simulation will fill
#using Distributed

include("auxiliary functions.jl")
include("parameters.jl")
include("Estimators.jl")

# Different covariances for the
covs = [0.1, 0.5, 0.9]
# Number of estimators used in the simulation (x2 since we include the standard errors. WMDF doesn't include the standard error... yet)
n_estimators = 14 + 1
# how many bandwidth columns we report
band = 5  # I want to report all the bandwidths: the one we loopover, 1 from using np, 1 from my loo criterion, and 2 coming from the AIC criteria
mat = SharedArray(zeros(nsim*length(bandrange), (4 + n_estimators + band), length(N), length(SIGMA), length(fs_fun)) )
start = time()

# Threads.@threads
#@time begin
# @distributed
 for i = 1:nsim
        for j = 1:length(SIGMA)
        for n in 1:length(N)
		for f in 1:length(fs_fun)
    data = dgp(N[n], beta_0, mu, SIGMA[j], design = fs_fun[f])
    # Got loop over the bandwidth here to keep the same data set throughout

    OLS = ols(data.x, data.y, intercept = "false")
    TSLS = tsls(data.x, data.y, data.z, intercept = "false")
#    bw = myest_bw(data.x, data.y, data.z; regtype = "lc")

    bw_np = myest_bw(data.x, data.y, data.z; regtype = "ll").bandwidth[1]
    np_2SLS = myest(data.x, data.y, data.z, bw_np; regtype = "ll")

# AIC version for the IV case
aic_bw = min_AIC_iv( data.y, data.x, data.z; lower = 0.01, upper = 1000, type = "ll").opt_h
aic =  myest(data.x, data.y, data.z, aic_bw; regtype = "ll")

# AIC version on the first stage
aic_bw_stage1 = min_AIC( data.x, data.z; lower = 0.01, upper = 1000, type = "ll").opt_h
aic_stage1 =  myest(data.x, data.y, data.z, aic_bw_stage1; regtype = "ll")

# WMDF
wmd = wmdf(data.x, data.y, data.z; intercept = "false").WMDF


# Add that to the simulation
kernel_TV = loo_opti(data.x, data.y, data.z, regtype = "ll")
bw_TV = kernel_TV.Bandwidth[1]
# Maybe for later, the whole formatting thing is still tricky
#    kernel_TV_h = map(h -> myest( data.x, data.y, data.z, h, regtype = "lc"), bandrange )
#    gna = groupby(kernel_TV_h, :bandwidth)
#    combine(gna), :estimate => mean)

     # NP = np_tsls(data.x, data.y, data.z, option = locallinear)
     # LIML = liml(data.x, data.y, data.z)
     # WMDF = wmdf(data.x, data.y, data.z; intercept = "false").WMDF
     # res_lc = TV_kernel_iv(data.y, data.x, data.z, option = "lc")
     # res_ll = TV_kernel_iv(data.y, data.x, data.z, option = "ll")
     # TV_lc = res_lc.bt_hat
     # TV_ll = res_ll.bt_hat
     # opt_h_lc = res_lc.opt_h
     # opt_h_ll = res_ll.opt_h
     # TV_lc_h = map(bw -> kernel_tv( data.x, data.y, data.z, "lc", bw), bandrange )
     # TV_ll_h = map(bw -> kernel_tv( data.x, data.y, data.z, "ll", bw), bandrange )
     bhat_ols = OLS.Estimates[1]
     bhat_tsls = TSLS.Estimates[1]
     bhat_np = np_2SLS.Estimates[1]
	 bhat_TV = kernel_TV.Estimates[1]
     bhat_aic = aic.Estimates[1]
     bhat_aic_stage1 = aic_stage1.Estimates[1]
     bhat_wmdf = wmd

     bhat_ols_se = OLS."Standard errors"[1]
     bhat_tsls_se = TSLS."Standard errors"[1]
     bhat_np_se = np_2SLS."Standard errors"[1]
     bhat_TV_se = kernel_TV."Standard errors"[1]
     bhat_aic_se = aic."Standard errors"[1]
     bhat_aic_stage1_se = aic_stage1."Standard errors"[1]



# Loop over the bandwidths
     for h in 1:length(bandrange)
         np_2SLS_h = myest(data.x, data.y, data.z, bandrange[h]; regtype = "ll")
         bhat_np_h = np_2SLS_h.Estimates[1]
         bhat_np_h_se = np_2SLS_h."Standard errors"[1]

    mat[(i-1)*length(bandrange) + h, 1, n, j, f] = i
	mat[(i-1)*length(bandrange) + h, 2, n, j, f] = N[n]
    mat[(i-1)*length(bandrange) + h, 3, n, j, f] = covs[j]
	mat[(i-1)*length(bandrange) + h, 4, n, j, f] = f
	mat[(i-1)*length(bandrange) + h, 5, n, j, f] = bhat_ols
    mat[(i-1)*length(bandrange) + h, 6, n, j, f] = bhat_ols_se

    mat[(i-1)*length(bandrange) + h, 7, n, j, f] = bhat_tsls
    mat[(i-1)*length(bandrange) + h, 8, n, j, f] = bhat_tsls_se

    mat[(i-1)*length(bandrange) + h, 9, n, j, f] = bandrange[h]

    mat[(i-1)*length(bandrange) + h, 10, n, j, f] = bhat_np_h
    mat[(i-1)*length(bandrange) + h, 11, n, j, f] = bhat_np_h_se


    mat[(i-1)*length(bandrange) + h, 12, n, j, f] = bhat_np
    mat[(i-1)*length(bandrange) + h, 13, n, j, f] = bhat_np_se

	mat[(i-1)*length(bandrange) + h, 14, n, j, f]= bhat_TV
    mat[(i-1)*length(bandrange) + h, 15, n, j, f]= bhat_TV_se

    mat[(i-1)*length(bandrange) + h, 16, n, j, f]= bhat_aic
    mat[(i-1)*length(bandrange) + h, 17, n, j, f]= bhat_aic_se

    mat[(i-1)*length(bandrange) + h, 18, n, j, f]= bhat_wmdf

    mat[(i-1)*length(bandrange) + h, 19, n, j, f]= bhat_aic_stage1
    mat[(i-1)*length(bandrange) + h, 20, n, j, f]= bhat_aic_stage1_se

    mat[(i-1)*length(bandrange) + h, 21, n, j, f]= aic_bw
    mat[(i-1)*length(bandrange) + h, 22, n, j, f]= aic_bw_stage1
    mat[(i-1)*length(bandrange) + h, 23, n, j, f]= bw_np
    mat[(i-1)*length(bandrange) + h, 24, n, j, f]= bw_TV


end
end
end
end
end

#end

elapsed = time() - start
# On nsim = 2
# no parallel: 3.4 sec
# using @distributed:

# Results formatting time
ncol = size(mat, 2)

mata = Array{Float64}(undef, 1, ncol)   # Make a row to concatenate the results vertically
for i in 1:length(fs_fun), j in 1:length(SIGMA), k in 1:length(N)
global mata = vcat(mata, mat[:,:, k, j, i])   # The "global" is needed on recent versions of Julia since otherwise, mata is reassigned every iteration and Julia does not recognize it
end
mat = mata[2:end,:]
#mat = vcat(mat[:,:,:, 1], mat[:,:,:, 2], mat[:,:,:, 3])
#mat = vcat(mat[:, :, 1], mat[:, :, 2])

results = DataFrame( mat , :auto)
results[!, :x4] = Int64.(results[:, :x4])
results[!, :x4] = fs_fun[results[:, 4]]
#ll_h_names = map(k -> name_k(k, "ll_h"), bandrange)
#lc_h_names = map(k -> name_k(k, "lc_h"), bandrange)

# To do:
# - Change names to match the R simulations
# - Add WMDF
# - Add the variances
# - Double check the cross validation for betahat formula
# - Rename current results
col_names = ["s", "sample", "covuv", "fs_fun",
             "ols",
             "ols_se",
             "tsls",
             "tsls_se",
             "bandwidth",
			 "ll",
             "ll_se",
            # "lc",
            # "lc_se"
             "np_tsls",
             "np_tsls_se",
             "kernel_TV",
             "kernel_TV_se",
             "bhat_aic",
             "bhat_aic_se",
             "WMDF",
             "bhat_aic_stage1",
             "bhat_aic_stage1_se",
             "bw_aic",
             "bw_aic_stage1",
             "bw_np",
             "bw_TV"
              ]
#col_names = vcat(col_names, ll_h_names, lc_h_names)
#names!(results, [Symbol("$i") for i in col_names ])
results = rename!(results, col_names)
results

#writedlm("caca.csv", bet, delim = ',')
CSV.write("F:/Kernel 2SLS/Julia simulations/first results (unfinished).csv", results)
