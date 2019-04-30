#=
This function "Updates" the parameters that are used
to draw from the posterior mean
=#


function μ_post_par(n_states, κ_hyper, μ_hyper, n_obs_in_each_state, grouped_obs, σ2_post_draw)
    post_μ = []
    post_σ2_output = []
    for m=1:n_states
        n = n_obs_in_each_state[m] #Just for convenience
        if n >= 1
            post_μ_temp = κ_hyper[m]/(κ_hyper[m]+n)*μ_hyper[m] + n/(κ_hyper[m]+n)*mean(grouped_obs[m])
            append!(post_μ, post_μ_temp)
            post_σ2_output_temp = σ2_post_draw[m]/(κ_hyper[m]+n)
            append!(post_σ2_output, post_σ2_output_temp)
        else
            post_μ_temp = μ_hyper[m]
            append!(post_μ, post_μ_temp)
            post_σ2_output_temp = σ2_post_draw[m]/(κ_hyper[m])
            append!(post_σ2_output, post_σ2_output_temp)
        end
    end
    return post_μ, post_σ2_output
end
