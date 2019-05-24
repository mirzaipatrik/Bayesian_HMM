# Bayesian Hidden Markov Models in Julia -- A Quick Overview

### Author: Patrik Mirzai

This project can be used to perform Bayesian parameter estimation of hidden Markov models. A brief introduction to the different functions are given below. Some applications include:

- Estimation of the mean and the variance of the different states
- Use posterior draws from the latent state sequence to estimate the time point for a regime change in data


## Bayesian parameter estimation:

Let's start with initializing the functions that are used to run the Gibbs sampler:

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

## Parameter estimation

Let us now perform a Bayesian parameter estimation of the HMM. We need to specify the number of states, the number of iterations, the burn-in period, the initial transition probaiblity matrix, the raw data, and the hyperparameters:

```julia
#Inputs:
n_states = 3  #Number of states
n_iter = 2000  #Number of iterations
Γ = [0.1 0.7 0.2; 0.3 0.4 0.3; 0.2 0.2 0.6];  #Transition prob. matrix
#These are our hyperparameters:
υ_hyper = [1,1,1]
σ2_hyper = [1,1,1]
μ_hyper = [4, 5, 6]
κ_hyper = [1,1,1]

#The burn-in values are used for computing the AIC and BIC
burn_in = 1000

MCMC_sim = main_function(n_states, n_iter, burn_in, Γ, dat, υ_hyper, σ2_hyper, μ_hyper, κ_hyper)
```

We can now plot the MCMC chain. A plot of the posterior draws of the different means is given. The posterior draws of the variances can be obtained analogously.

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

Retrieve AIC and BIC scores:

```julia
#AIC and BIC scores
AIC_1 = MCMC_sim[7]
BIC_1 = MCMC_sim[8]
```


## Estimation of the time point for a regime change
Let us now estimate the time point for a regime change in the data. Let's run the Gibbs sampler again:

```julia
#Read the new data set:
regime_data = CSV.read("regime_data.txt")

#Plot the data set:
plot(regime_data[:,1], color="blue", label="", linewidth=2.5, ylabel="Volume [nl]", xlabel="Observation", guidefont=a, titlefont=a, tickfont=a, legendfont=a)
```
![grouped](https://github.com/mirzaipatrik/Bayesian_HMM/blob/master/Functions/BGA_plot.png)



```julia
#Set initial transition probability matrix:
Γ_2 = [0.5 0.5; 0.5 0.5]

#Setting hyperparameters:
υ_hyper = [1,1]
σ2_hyper = [1,1]
μ_hyper = [0.5, 1.8]
κ_hyper = [1,1]

#Run the Gibbs sampler:
new_test = main_function(2, 11000, 1000, Γ_2, regime_data[:,1], υ_hyper, σ2_hyper, μ_hyper, κ_hyper)

#=
ξ is the time point for a regime change in data.
The Gibbs sampler allows for a posterior distribution of this quanttity.
=#
ξ = new_test[5][1001:11000]  #Discard the initial 1000 draws as "burn-in" values.

#Histogram over the posteior distribution:
histogram(ξ, guidefont=a, titlefont=a, tickfont=a, legendfont=a, xlims=[60000, 105000], label="", xlabel="ξ", ylabel="Frequency", color="blue", bins=50)
```

![grouped](https://github.com/mirzaipatrik/Bayesian_HMM/blob/master/Functions/Stable_state.png)


## References
Following references have been used in the implementation of the Gibbs sampler:

- Cappé, O., Moulines, E. and Rydén, T. (2005), Inference in Hidden Markov Models, New York: Springer-Verlag.

- Gelman, A., Carlin, J., Stern, H., Dunson, D., Vehtari, A. and Rubbin, D. (2014), Bayesian Data Analysis, third edn, Taylor & Francis Group.

- Rydén, T. (2008), ‘EM versus Markov chain Monte Carlo for estimation of hidden Markov models: a computational perspective’, Bayesian Analysis, Vol. 3.

- Zucchini, W. and MacDonald, I. (2009), Hidden Markov Models for Time Series:
An Introduction Using R, Chapman & Hall/CRC.
