# Bayesian Hidden Markov Models in Julia -- A Quick Overview

### Author: Patrik Mirzai

This project can be used to perform Bayesian parameter estimation of hidden Markov models. A brief introduction to the different function are given below. Some applications include:

- Estimation of the mean and variance of the different states
- Use posterior draws from the latent state sequence to estimate the time point for a regime change in data


Initialize:

```julia
using CSV
using Distributions
using Statistics
using Plots
using Random
using LinearAlgebra
using DataFrames
using Pycall

#Import local functions
include("MC_simulation.jl")
include("state_counts.jl")
include("grouped_obs.jl")
include("sigma_post_par.jl")
include("posterior_variance_draw.jl")
include("mu_post_par.jl")
include("posterior_mean_draws.jl")
include("transition_counter.jl")
include("initial_dist_post.jl")
include("state_dep_diag.jl")
include("backward_function.jl")
include("main_function.jl")

#Optional plotting options:
a=font(20,"Computer Modern")
b=font(18,"Computer Modern")
c=font(16,"Computer Modern")
d=font(14,"Computer Modern")
pyplot(dpi=300)  #Set backend for plotting -- gives high resolution when saving the plot

```

We can now import and plot data


```julia
dat = CSV.read("My_data.csv", header=false);  #Read the data
plot(dat[:,1])  #Plot the data

#=
Optional for saving figure to directory:
savefig("name_of_figure")
=#
```

![grouped](https://github.com/mirzaipatrik/Bayesian_HMM/blob/master/Functions/Traceplot.png)

# Parameter estimation

Let us now perform a Bayesian parameter estimation of the HMM. We need to specify the number of states, number of iterations and the data, the initial transition probaiblity matrix, the raw data and the hyperparameters:

```julia
#Inputs:
n_states = 3  #Number of states
n_iter = 2000  #Number of iterations
Γ = [0.1 0.7 0.2; 0.3 0.4 0.3; 0.2 0.2 0.6];  #Transition prob. matrix
#These are our hyperparameters:
υ_hyper = [1,1,1]
σ2_hyper = [1,1,1]
μ_hyper = [0.8, 1, 1.2]
κ_hyper = [1,1,1]

MCMC_sim = main_function(n_states, n_iter, Γ, dat, υ_hyper, σ2_hyper, μ_hyper, κ_hyper)
```

We can now plot the MCMC chain. A plot of the posterior draws of the different means is given, but posterior draws of the variances can be obtained analogously.

```julia
#Retrieve the matrix of posterior draws of μ and σ2:
μ = MCMC_sim[1]
σ2 = MCMC_sim[2]

#We use a "burn-in" period of 200 when plotting the posterior draws from the Gibbs sampler:
plot(μ[200:length(μ[:,1]),1], xlabel="MCMC iteration", ylabel="μ [nl]", label="μ1", color="blue", linewidth=1.5, guidefont=b, titlefont=b, tickfont=b, legendfont=b, title="", ylim=[5, 8])
plot!(μ[200:length(μ[:,1]),2], label="μ2", color="red", linewidth=1.5)
plot!(μ[200:length(μ[:,1]),3], label="μ3", color="grey", linewidth=1.5)
```
![grouped](https://github.com/mirzaipatrik/Bayesian_HMM/blob/master/Functions/posterior_mean_draws.png)


# References
Following references have been used in the implementation of the Gibbs sampler:

- Cappé, O., Moulines, E. and Rydén, T. (2005), Inference in Hidden Markov Models, New York: Springer-Verlag.

- Rydén, T. (2008), ‘EM versus Markov chain Monte Carlo for estimation of hidden Markov models: a computational perspective’, Bayesian Analysis, Vol. 3.

- Zucchini, W. and MacDonald, I. (2009), Hidden Markov Models for Time Series:
An Introduction Using R, Chapman & Hall/CRC.
