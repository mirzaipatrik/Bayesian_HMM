#=
This function draws posterior μ-values
=#

function posterior_mean_draws(μ, σ2)
    mean_draws = []
    for k = 1:length(μ)
        append!(mean_draws, rand(Normal(μ[k], sqrt(σ2[k]))))
    end
    return mean_draws
end
