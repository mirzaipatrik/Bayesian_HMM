#=
The packages are imported with the "using" function. However,
the first time the packages are used, they need to be installed.
A simple way of installing the packages is to first use "using Pkg" module
and then use "Pkg.add("package name")".
=#
using CSV
using Distributions
using Statistics
using Plots
using Random
using LinearAlgebra
using DataFrames

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


#=
If pycall is not installed, install python packages by:

Pkg.build("PyCall")

import Pkg
Pkg.add("PyCall")

using PyCall


Pkg.add("Conda")
=#


#=
Description of variables:
Γ:= transition probability matrix (tpm)
γ:= transition probability
μ:= mean vector
σ2 := variance vector
=#
my_dat_2 = CSV.read("My_data.csv", header=true, delim=";", decimal=',')
plot(my_dat_2[:,4])
empty_vec = []
for i in 3802:(length(my_dat_2[:,4])-263)
    if my_dat_2[i,4]<=10
        print(my_dat_2[i,4])
        append!(empty_vec, my_dat_2[i,4])
    end
end

empty_vec

mod(48960, 360)

print(48960/360)

converted_dat = convert(Array{Float64,1}, empty_vec)
plot(empty_vec[1:100])
median(empty_vec)
plot(empty_vec, color="blue", label="", ylabel="Volume [nl]", xlabel="Observation", guidefont=a, titlefont=a, tickfont=a, legendfont=a)
histogram(empty_vec, color="blue", bins=80, normed=true, label="",
 xlabel="Volume [nl]", ylabel="Density", guidefont=a, titlefont=a, tickfont=a, legendfont=a)
savefig("BGA_histogram")
savefig("Traceplot")

#=
#test_datt = CSV.read("Test_data.csv", header=true, delim=";", decimal=',')
my_dat = CSV.read("My_data_1.csv", header=true, delim=";", decimal=',')
my_dat_2 = CSV.read("My_data_5.csv", header=true, delim=";", decimal=',')

empty_vec = []
for i in 1:length(my_dat_2[:,4])
    if my_dat_2[i,4]<=10
        print(my_dat_2[i,4])
        append!(empty_vec, my_dat_2[i,4])
    end
end

plot(empty_vec)
histogram(empty_vec, bins=100, color="blue")



plot(my_dat[:, 4])
plot(my_dat_2[:,4])

=#



#Read and plot data for visualization
mycronic_data = CSV.read("Mycronic_dat.txt", header=false);
dat = mycronic_data[:,1];
pyplot(dpi=300)  #set backend
plot(dat, color="blue")
a = histogram(dat[1:10000], bins=80, normed=true, color="blue", label="",
 xlabel="Volume [nl]", ylabel="Density");
plot(a)

plot(dat[1:50],color="blue", label="", linewidth=2.5, ylabel="Volume [nl]", xlabel="Observation", guidefont=a, titlefont=a, tickfont=a, legendfont=a)
scatter!(dat[1:50], color="blue", label="", markercolor = :blue, markerstrokecolor = :blue, markersize = 5)
a=font(18,"Computer Modern")
b=font(12,"Computer Modern")

savefig("BGA_plot")

histogram(dat[1:10000,1], label="", color="blue", ylabel="Density", xlabel="Volume [nl]", guidefont=a, titlefont=a, tickfont=a, legendfont=a)

savefig("hiiist.png")
#Initialize parameters
Γ = [0.1 0.7 0.2; 0.3 0.4 0.3; 0.2 0.2 0.6];  #Transition prob. matrix
#Define hyperparameters


κ_hyper = [1,1,1]
σ2_hyper = [1,1,1]
υ_hyper = [1,1,1]
μ_hyper = [0.8, 1, 1.2]


Γ_new = [0.5 0.5; 0.5 0.5]
ouu = main_function(2, 400, Γ_new, dat, [1,1], [1,1], [1,1], [1,1])
print(ouu[4])

print([0.864554 0.135446; 0.09155 0.90845]^(1000))

var(dat)
lllllllll = convert(Float64, dat)

using StatsBase

looooo = ecdf(dat)
mean_vector = ouu[2]
a=font(18,"Computer Modern")
b=font(12,"Computer Modern")
len = length(mean_vector[:,1])
plot(mean_vector[100:len,1], ylim=[1.2, 1.7], ylabel="Nanoliter", xlabel="MCMC iteration", label="μ1", color="blue", linewidth=3, guidefont=b, titlefont=b, tickfont=b, legendfont=b)
plot!(mean_vector[100:len,2], label="μ2", color="red", linewidth=3)
plot!(mean_vector[100:len,3], label="μ3", color="black", linewidth=3)
savefig("Traceplot1")


#=
Following two lines are optional in order to write the output to a CSV file:
df =  DataFrame(mean_vector)  #Create a dataframe for outout
CSV.write("MCMC_mean_chain.csv", df)  #Create the file
=#


histogram(mean_vector[1000:len,1], bins=80, normed=false, color="blue")
histogram(mean_vector[1000:len,2], bins=50, normed=false, color="blue")
histogram(mean_vector[1000:len,3], bins=80, normed=false, color="blue")

σ2_vector = ouu[1]
plot(σ2_vector[100:400,1], xlabel="MCMC iteration", label="σ1", color="blue", linewidth=3, guidefont=b, titlefont=b, tickfont=b, legendfont=b)
plot!(σ2_vector[100:400,2], label="σ2", color="red", linewidth=3, guidefont=b, titlefont=b, tickfont=b, legendfont=b)
plot!(σ2_vector[100:300,3], label="σ3", color="black", linewidth=3, guidefont=b, titlefont=b, tickfont=b, legendfont=b)
savefig("Variance_plot")

random_1 = rand(Normal(-2, 0.5), 1000)
random_2 = rand(Normal(0, 0.5), 1000)
random_3 = rand(Normal(2, 0.5), 1000)


arr = vcat(random_1, random_2, random_3)
arr = shuffle(arr)

test_MC = main_function(3, 1000, Γ, arr, υ_hyper, σ2_hyper, μ_hyper, κ_hyper)
μ = test_MC[2]
a=font(18,"Computer Modern")
b=font(12,"Computer Modern")
plot(μ[100:length(μ[:,1]),1], xlabel="MCMC iteration", label="μ1", color="blue", linewidth=1.2, guidefont=b, titlefont=b, tickfont=b, legendfont=b, title="sd=1")
plot!(μ[100:length(μ[:,1]),2], label="μ2", color="red", linewidth=1.2)
plot!(μ[100:length(μ[:,1]),3], label="μ3", color="grey", linewidth=1.2)


χ = test_MC[3]
plot(χ)

#plot(lo, lo_2, layout = (2))
savefig("Trae")


quantile(μ[100:length(μ[:,1]),1], [0.025, 0.975])
quantile(μ[100:length(μ[:,1]),2], [0.025, 0.975])
quantile(μ[100:length(μ[:,1]),3], [0.025, 0.975])

σ = test_MC[1]
quantile(σ[200:length(μ[:,1]),1], [0.025, 0.975])
quantile(σ[200:length(μ[:,1]),2], [0.025, 0.975])
quantile(σ[200:length(μ[:,1]),3], [0.025, 0.975])



histogram(μ[1000:10000,], bins=80)

σ2 = test_MC[1]
quantile(σ2[1000:10000,1], [0.025, 0.975])
quantile(σ2[1000:10000,2], [0.025, 0.975])
quantile(σ2[1000:10000,3], [0.025, 0.975])

print(test_MC[4])

#Read the MCMC plot

mean_long = CSV.read("MCMC_mu.csv")
var_long = CSV.read("MCMC_var.csv")

plot(mean_long[1000:50000,1], color="blue", ylim=[1.25, 1.55])
plot!(mean_long[1000:50000,2], color="red")
plot!(mean_long[1000:50000,3], color="black")

plot(var_long[1000:50000,1], color="blue")
plot!(var_long[1000:50000,2], color="red")
plot!(var_long[1000:50000,3], color="black")



regime_data = CSV.read("regime_data.txt")
first_plot = plot(regime_data[1:99360,1], color="blue", guidefont=a, titlefont=a, tickfont=a, legendfont=a, label="", ylabel="Volume [nl]", xlabel="Observation")
savefig("regime_plot")
Γ_newest = [0.5 0.5; 0.5 0.5]
newest_test = main_function(2, 2000, Γ_newest, regime_data[:,1], [1,1], [1,1], [0.5, 1.8], [1,1])
ξ = newest_test[5]
histogram(ξ)
plot(μ_vec[100:2000,1], color="blue")
plot!(μ_vec[100:2000,2], color="red")
plot!(μ_vec[100:10000,3], color="grey")

for i in 1:length(newest_test[3])
    if newest_test[3][i] == 1
        newest_test[3][i] = 2
    else
        newest_test[3][i] == 1
    end
end
end


first_plot = plot(regime_data[1:99360,1], color="blue", guidefont=a, titlefont=a, tickfont=a, legendfont=a, label="", ylabel="Volume [nl]", xlabel="Observation")
second_plot = scatter(newest_test[3], color="blue", markercolor = :blue, markerstrokecolor = :blue, markersize = 5, guidefont=a, titlefont=a, tickfont=a, legendfont=a, xlabel="Iteration", ylabel="State", yticks=[1,2])
plot(first_plot, second_plot, layout=(2,1), legend=false)
savefig("test")
mean(regime_data[60000:95000,1])


plot(newest_test[3], color="blue")


plot([regime_data[:,1], newest_test[3]], layout = (2,1), color="blue", linewidth=2, guidefont=a, titlefont=a, tickfont=a, legendfont=a)
savefig("regime_plot")
print(newest_test[5])

histogram(newest_test[5][100:400], color="blue", label="hello", bins=50, guidefont=a, titlefont=a, tickfont=a, legendfont=a)
savefig("Stable_state")


plot(regime_data[55000:55100,1])

Pkg.rm("rm -m DebuggerFramework")
Pkg.add("HypothesisTests")

using HypothesisTests
using RecipesBase






Pkg.add("StatsPlots")
using StatsPlots

gr(dpi=300)
using Distributions
x = rand(Normal(), 100)
y = rand(Cauchy(), 100)
qqplot(x, y, qqline = :fit)
plot(
 qqplot(x, y, qqline = :fit), # qqplot of two samples, show a fitted regression line
 qqplot(Cauchy, y),           # compare with a Cauchy distribution fitted to y; pass an instance (e.g. Normal(0,1)) to compare with a specific distribution
 qqnorm(x, qqline = :R)       # the :R default line passes through the 1st and 3rd quartiles of the distribution
)

pyplot(dpi=300)
qqnorm(converted_dat, qqline = :R, xlabel="Theoretical quantiles", ylabel="Sample quantiles", linewidth=3.5, linecolor="red", markercolor = :blue, markerstrokecolor = :blue, guidefont=a, titlefont=a, tickfont=a, legendfont=a)
savefig("QQ_plot")


typeof(x)
typeof(dat)

new_dat = convert(Array{Float64}, dat)

using HypothesisTests
mod(345600, 360)
ExactOneSampleKSTest(converted_dat, Normal(mean(converted_dat), sqrt(var(converted_dat))))

1.36/sqrt(48960)
