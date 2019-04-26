#=
We create a function that computes the posterior parameters for
σ^2_j, where j ∈ {1,2,...m}
=#

function σ2_post_par(n_states, n_obs_in_each_state, υ_hyper, σ2_hyper, κ_hyper, μ_hyper, grouped_obs)
    #Vectors for storing the posterior parameters
    post_υ = []
    post_σ2 = []

    for m = 1:n_states
        if n_obs_in_each_state[m] > 1
            υ_update = υ_hyper[m]+n_obs_in_each_state[m]
            append!(post_υ, υ_update)
            σ2_update = (1/post_υ[m])*((υ_hyper[m]*σ2_hyper[m])+(n_obs_in_each_state[m]-1)*var(grouped_obs[m])+((κ_hyper[m]*n_obs_in_each_state[m])/(κ_hyper[m]+n_obs_in_each_state[m]))*((mean(grouped_obs[m])-μ_hyper[m])^2))
            append!(post_σ2, σ2_update)
        elseif n_obs_in_each_state[m] == 1
            υ_update = υ_hyper[m]+n_obs_in_each_state[m]
            append!(post_υ, υ_update)
            σ2_update = (1/post_υ[m])*((υ_hyper[m]*σ2_hyper[m])+((κ_hyper[m]*n_obs_in_each_state[m])/(κ_hyper[m]+n_obs_in_each_state[m]))*((mean(grouped_obs[m])-μ_hyper[m])^2))
            append!(post_σ2, σ2_update)
        else
            υ_update = υ_hyper[m]+n_obs_in_each_state[m]
            append!(post_υ, υ_update)
            σ2_update = (1/post_υ[m])*(υ_hyper[m]*σ2_hyper[m])
            append!(post_σ2, σ2_update)
        end
    end
    return post_υ, post_σ2
end
