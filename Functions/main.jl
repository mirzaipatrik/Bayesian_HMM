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
using Random

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
include("MCMC_Simulation.jl")

Random.seed!(1234)

data = CSV.read("./solder_paste_data.csv", DataFrame).Volume


pl = plot(data, label = "", line=(:line, :blue),
guidefontsize=10, tickfontsize=10, legendfontsize=10, xlims=(0, 52000), xlabel="Observation", ylabel="Volume [nl]")

savefig(pl, "manufacturing_data.png")

n_iter = 500  #Number of iterations
lps_iterations = 2000
burn_in = 200


n_states_sim = [2, 3, 4, 5, 6, 7]
Γ_sim = [
    [
        0.1 0.9;
        0.2 0.8
    ],
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
    ]
]
υ_hyper_sim = [
    [1, 1],
    [1, 1, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1],
]

σ2_hyper_sim =
    [
        [1, 1],
        [1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1],
    ]
μ_hyper_sim = [
    [4, 5],
    [4, 5, 6],
    [4, 5, 6, 7],
    [4, 5, 6, 7, 8],
    [4, 5, 6, 7, 8, 9],
    [4, 5, 6, 7, 8, 9, 10],
]
κ_hyper_sim = [
    [1, 1],
    [1, 1, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1],
]

results = Array{Any}(undef, 6)


observationss = MCMC_Simulation([0.8 0.1 0.1; 0.1 0.8 0.1; 0.1 0.1 0.8], [1/3, 1/3, 1/3], [4, 6, 8], [1, 1, 1], 10000)

Threads.@threads for i in 1:6
    results[i] = main_function(
        n_states_sim[i],
        1000,
        burn_in,
        Γ_sim[i],
        observationss,
        υ_hyper_sim[i],
        σ2_hyper_sim[i],
        μ_hyper_sim[i],
        κ_hyper_sim[i],
        2000
    )
end

results[1][7]
results[2][7]
results[3][7]
results[4][7]
results[5][7]
results[6][7]

results[1][8]
results[2][8]
results[3][8]
results[4][8]
results[5][8]
results[6][8]

results[1][9]
results[2][9]
results[3][9]
results[4][9]
results[5][9]
results[6][9]

results[1][9]
results[2][9]
results[3][9]
println(results[4][9])
println(results[5][9])
println(results[6][9])

results[2][6]


results[2][6]^1000

results[1][2][201:1000, :]
results[2][2][201:1000, :]
results[3][2][201:1000, :]
results[4][2][201:1000, :]
results[5][2][201:1000, :]
results[6][2][201:1000, :]

#results from a three state model: 

mean(results[2][2][201:1000, 1])
quantile(results[2][2][201:1000, 1], 0.975)
mean(results[2][2][201:1000, 3])
results[2][3]
results[2][4]
results[2][5]
results[2][6]
results[2][7]
results[2][8]
results[2][9]


p = plot(results[2][1][201:1000, 1], label="State 1", xlabel="MCMC Iteration", ylabel="\$\\mu\$",
    line=(:line, :red), guidefontsize=16, tickfontsize=14, legendfontsize=12, ylim=(3.8, 10))
plot!(results[2][1][201:1000, 2], label="State 2", line=(:line, :blue),
    guidefontsize=16, tickfontsize=14, legendfontsize=12)
plot!(results[2][1][201:1000, 3], label="State 3", line=(:line, :green),
    guidefontsize=16, tickfontsize=14, legendfontsize=12)
 
savefig(p, "plot.png")

# Do the same for the variance 
q = plot(results[2][2][201:1000, 1], label="State 1", xlabel="MCMC Iteration", ylabel="\$\\sigma^2\$",
    line=(:line, :red), guidefontsize=16, tickfontsize=14, legendfontsize=12, ylim=(0.1, 1.5))
plot!(results[2][2][201:1000, 2], label="State 2", line=(:line, :blue),
    guidefontsize=16, tickfontsize=14, legendfontsize=12)
plot!(results[2][2][201:1000, 3], label="State 3", line=(:line, :green),
    guidefontsize=16, tickfontsize=14, legendfontsize=12)

savefig(q, "trace_plots_3_states_variance.png")

data_Results = Array{Any}(undef, 6)

lpds_iterations_2 = 6000

Threads.@threads for i in 1:6
    data_Results[i] = main_function(
        n_states_sim[i],
        n_iter,
        burn_in,
        Γ_sim[i],
        data,
        υ_hyper_sim[i],
        σ2_hyper_sim[i],
        μ_hyper_sim[i],
        κ_hyper_sim[i],
        lpds_iterations_2
    )
end

# Write the lpds results to a csv file
lpds_results_from_simulation = DataFrame(
    n_states = n_states_sim,
    lpds = [data_Results[i][7] for i in 1:6]
);

# Write a csv file with the results
CSV.write("lpds_results_from_simulation.csv", lpds_results_from_simulation)

data_Results[1][7]
data_Results[2][7]
data_Results[3][7]
data_Results[4][7]
data_Results[5][7]
data_Results[6][7]

# Now that we have chosen our model, let us simulate from a larger posterior distribution
# We will simulate 10000 observations from the posterior distribution

posterior_sim = main_function(
    3,
    10000,
    burn_in,
    Γ_sim[2],
    data,
    υ_hyper_sim[2],
    σ2_hyper_sim[2],
    μ_hyper_sim[2],
    κ_hyper_sim[2],
    1
)

# Now let's write the results for the mean and variance to a csv file
posterior_sim_mean = DataFrame(
    state_1 = posterior_sim[1][201:10000, 1],
    state_2 = posterior_sim[1][201:10000, 2],
    state_3 = posterior_sim[1][201:10000, 3]
)

posterior_sim_variance = DataFrame(
    state_1 = posterior_sim[2][201:10000, 1],
    state_2 = posterior_sim[2][201:10000, 2],
    state_3 = posterior_sim[2][201:10000, 3]
)

# Write the results to a csv file
CSV.write("posterior_sim_mean.csv", posterior_sim_mean)
CSV.write("posterior_sim_variance.csv", posterior_sim_variance)

# Now let's write the quanties for the mean and the variance to a csv file
posterior_sim_quantiles_mean = DataFrame(
    state_1 = [quantile(posterior_sim[1][201:10000, 1], 0.025), quantile(posterior_sim[1][201:10000, 1], 0.975)],
    state_2 = [quantile(posterior_sim[1][201:10000, 2], 0.025), quantile(posterior_sim[1][201:10000, 2], 0.975)],
    state_3 = [quantile(posterior_sim[1][201:10000, 3], 0.025), quantile(posterior_sim[1][201:10000, 3], 0.975)]
)

posterior_sim_quantiles_variance = DataFrame(
    state_1 = [quantile(posterior_sim[2][201:10000, 1], 0.025), quantile(posterior_sim[2][201:10000, 1], 0.975)],
    state_2 = [quantile(posterior_sim[2][201:10000, 2], 0.025), quantile(posterior_sim[2][201:10000, 2], 0.975)],
    state_3 = [quantile(posterior_sim[2][201:10000, 3], 0.025), quantile(posterior_sim[2][201:10000, 3], 0.975)]
)



# Write the results to a csv file
CSV.write("posterior_sim_quantiles_mean.csv", posterior_sim_quantiles_mean)
CSV.write("posterior_sim_quantiles_variance.csv", posterior_sim_quantiles_variance)


# Return and print the tpm
tpm_post = posterior_sim[6]
CSV.write("tpm_post.csv", DataFrame(tpm_post, :auto))

mean(posterior_sim[2][201:10000, 3])
quantile(posterior_sim[2][201:10000, 3], 0.975)
posterior_sim[2]
posterior_sim[3]
posterior_sim[4]^10000
posterior_sim[5]

p = plot(posterior_sim[1][201:10000, 1], label="State 1", xlabel="MCMC Iteration", ylabel="\$\\mu\$",
    line=(:line, :red), guidefontsize=16, tickfontsize=14, legendfontsize=12, ylim=(5, 7.1), xlim=(0, 10500))
plot!(posterior_sim[1][201:10000, 2], label="State 2", line=(:line, :blue),
    guidefontsize=16, tickfontsize=14, legendfontsize=12)
plot!(posterior_sim[1][201:10000, 3], label="State 3", line=(:line, :green),
    guidefontsize=16, tickfontsize=14, legendfontsize=12)
 
savefig(p, "data_plot.png")

# Do the same for the variance 
q = plot(posterior_sim[2][201:10000, 1], label="State 1", xlabel="MCMC Iteration", ylabel="\$\\sigma^2\$",
    line=(:line, :red), guidefontsize=16, tickfontsize=14, legendfontsize=12, ylim=(0, 0.5), xlim=(0, 10500))
plot!(posterior_sim[2][201:10000, 2], label="State 2", line=(:line, :blue),
    guidefontsize=16, tickfontsize=14, legendfontsize=12)
plot!(posterior_sim[2][201:10000, 3], label="State 3", line=(:line, :green),
    guidefontsize=16, tickfontsize=14, legendfontsize=12)

savefig(q, "trace_plots_3_states_variance_data.png")
