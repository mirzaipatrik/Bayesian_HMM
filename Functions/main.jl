import Pkg;
Pkg.add("CSV");
Pkg.add("Distributions");
Pkg.add("Statistics");
Pkg.add("Random");
Pkg.add("LinearAlgebra");
Pkg.add("DataFrames");
Pkg.add("Plots");

using CSV
using Distributions
using Statistics
using Random
using LinearAlgebra
using DataFrames
using Plots

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

dat = CSV.read("Functions/My_data.csv", DataFrame, delim=';', decimal=',', ignoreemptylines=true)
# plot(dat[:, 4], line = (:line), label = "Volume")


n_states = 3  #Number of states
n_iter = 200  #Number of iterations
Γ = [0.1 0.7 0.2; 0.3 0.4 0.3; 0.2 0.2 0.6];  #Transition prob. matrix
#These are our hyperparameters:
υ_hyper = [1,1,1]
σ2_hyper = [1,1,1]
μ_hyper = [4, 5, 6]
κ_hyper = [1,1,1]
lps_iterations = 20


burn_in = 200

MCMC_sim = main_function(
    n_states,
    n_iter,
    burn_in,
    Γ,
    dat[:,4],
    υ_hyper,
    σ2_hyper,
    μ_hyper,
    κ_hyper,
    lps_iterations
)
