# import Pkg;
# Pkg.add("CSV");
# Pkg.add("Distributions");
# Pkg.add("Statistics");
# Pkg.add("Random");
# Pkg.add("LinearAlgebra");
# Pkg.add("DataFrames");
# Pkg.add("Plots");

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
# plot(dat[1:49000, 4], line = (:line), label = "Volume")

filtered_dat = dat[dat[:, 4].<8.3, :];
data = filtered_dat[4000:length(filtered_dat[:, 4]), 4]

n_states = 2  #Number of states
n_iter = 400  #Number of iterations
# Γ = [0.1 0.7 0.2; 0.3 0.4 0.3; 0.2 0.2 0.6];  #Transition prob. matrix
Γ = [0.1 0.9; 0.2 0.8]
#These are our hyperparameters:
# υ_hyper = [1,1,1]
# σ2_hyper = [1,1,1]
# μ_hyper = [4, 5, 6]
# κ_hyper = [1,1,1]
υ_hyper = [1, 1]
σ2_hyper = [1, 1]
μ_hyper = [4, 5]
κ_hyper = [1, 1]
lps_iterations = 3000


burn_in = 200

MCMC_sim = main_function(
    n_states,
    n_iter,
    burn_in,
    Γ,
    data,
    υ_hyper,
    σ2_hyper,
    μ_hyper,
    κ_hyper,
    lps_iterations
)

n_states_sim = [3, 4, 5, 6, 7, 8]
Γ_sim = [
    [
        0.4 0.2 0.4;
        0.4 0.2 0.4;
        0.4 0.2 0.4
    ],
    [
        0.1 0.7 0.1 0.1;
        0.2 0.6 0.1 0.1;
        0.1 0.1 0.6 0.2;
        0.1 0.1 0.2 0.6
    ],
    [
        0.1 0.6 0.1 0.1 0.1;
        0.1 0.6 0.1 0.1 0.1;
        0.1 0.6 0.1 0.1 0.1;
        0.1 0.6 0.1 0.1 0.1;
        0.1 0.6 0.1 0.1 0.1
    ],
    [
        0.1 0.5 0.1 0.1 0.1 0.1;
        0.1 0.5 0.1 0.1 0.1 0.1;
        0.1 0.5 0.1 0.1 0.1 0.1;
        0.1 0.5 0.1 0.1 0.1 0.1;
        0.1 0.5 0.1 0.1 0.1 0.1;
        0.1 0.5 0.1 0.1 0.1 0.1
    ],
    [
        0.1 0.1 0.1 0.1 0.1 0.1 0.4; 0.1 0.1 0.1 0.1 0.1 0.1 0.4;
        0.1 0.1 0.1 0.1 0.1 0.1 0.4; 0.1 0.1 0.1 0.1 0.1 0.1 0.4;
        0.1 0.1 0.1 0.1 0.1 0.1 0.4; 0.1 0.1 0.1 0.1 0.1 0.1 0.4;
        0.1 0.1 0.1 0.1 0.1 0.1 0.4
    ],
    [
        0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3; 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3;
        0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3; 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3;
        0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3; 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3;
        0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3; 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3
    ]
]
υ_hyper_sim = [
    [1, 1, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1]
]

σ2_hyper_sim = [
    [1, 1, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1]
]
μ_hyper_sim = [
    [4, 5, 6],
    [4, 5, 6, 7],
    [4, 5, 6, 7, 8],
    [4, 5, 6, 7, 8, 9],
    [4, 5, 6, 7, 8, 9, 10],
    [4, 5, 6, 7, 8, 9, 10, 11]
]
κ_hyper_sim = [
    [1, 1, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1]
]

results = Array{Any}(undef, 6)

Threads.@threads for i in 1:6
    results[i] = main_function(
        n_states_sim[i],
        n_iter,
        burn_in,
        Γ_sim[i],
        data,
        υ_hyper_sim[i],
        σ2_hyper_sim[i],
        μ_hyper_sim[i],
        κ_hyper_sim[i],
        lps_iterations
    )
end
