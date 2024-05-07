function main_function(
    n_states,
    n_iter,
    burn_in,
    Γ,
    raw_data,
    υ_hyper,
    σ2_hyper,
    μ_hyper,
    κ_hyper,
    lps_iterations
)
    global_counter = 0
    LPS = 0

    for j in 1:lps_iterations
        global_counter = global_counter + 1
        #Allocate memory for μ and σ^2 draws
        σ2_post_draws = zeros(n_iter, n_states)
        μ_post_draws = zeros(n_iter, n_states)

        #Create an array that is used to store latest change of states
        latest_change = 1
        state_register = []

        #Empty vector that stores the likelihood for each sweep
        likelihood_vec = []

        #Create matrix in order to compute posterior tpm, Γ
        Γ_output = zeros(n_states, n_states)

        n_obs = length(raw_data) - lps_iterations + global_counter
        MC_chain = MC_simulation(Γ, n_obs)  #Init a random MC


        #start our gibbs sampler
        for i = 1:n_iter
            print("Iteration: ", i, "\n")
            #Count number of visits to a state
            state_count = state_counts(n_states, MC_chain)
            #Group observations according to the MC_chain
            grouped_obser = grouped_obs(n_states, MC_chain, raw_data)
            #Update posterior parameters
            σ2_post_pars = σ2_post_par(n_states, state_count, υ_hyper, σ2_hyper, κ_hyper, μ_hyper, grouped_obser)

            #Draw from the posterior variance
            υ_n = σ2_post_pars[1]  #Array with posterior υ-parameters
            σ2_n = σ2_post_pars[2] #Array with posterior σ^2-parameters, ∼Scaled-inv χ2
            σ2_post_draws[i, :] = posterior_variance_draw(υ_n, σ2_n)
            latest_σ2 = σ2_post_draws[i, :]

            #"Update" posterior μ-parameters
            μ_post_pars = μ_post_par(n_states, κ_hyper, μ_hyper, state_count, grouped_obser, latest_σ2)
            #Draw posterior mean
            μ_n = μ_post_pars[1]
            μ_var = μ_post_pars[2]
            μ_post_draws[i, :] = posterior_mean_draws(μ_n, μ_var)
            latest_μ = μ_post_draws[i, :]

            #Update Γ (tpm)
            transition_counter_mat = transition_counter(MC_chain, n_states)
            dir_pars = transition_counter_mat .+ 1  #"Update" posterior γ_(ij)
            #Update Γ by drawing from the posterior (Dirichlet) distributoin
            for k = 1:n_states
                Γ[k, :] = rand(Dirichlet(dir_pars[k, :]))
            end  #End the for-loop

            #Update π (initial distribution)
            initial_pars = initial_dist_post(n_states, MC_chain)
            initial_dist = rand(Dirichlet(initial_pars))

            #----------------------------------------------
            #Update the latent chain

            #Compute the backward probabilities
            back_mat = backward_function(n_states, raw_data, Γ, latest_μ, latest_σ2)

            #Append the likelihood to the vector
            append!(likelihood_vec, back_mat[2])

            #Let us now do the forward simulation in the Gibbs sampler

            temp_prob_vec = fill(0.0, n_states)  #A placeholder
            for m = 1:n_states
                temp_prob_vec[m] = initial_dist[m] * (pdf(Normal(latest_μ[m], sqrt(latest_σ2[m])), raw_data[1])) * back_mat[1][1, m]
            end
            #Normalize temp_vec
            temp_prob_vec = temp_prob_vec / sum(temp_prob_vec)

            #Update state 1
            MC_chain[1] = rand(Categorical(temp_prob_vec))

            #Simulate the rest of the states using forward simulation
            for j = 2:n_obs
                latest_state = MC_chain[j-1]  #Use to condition on the current state

                for m = 1:n_states  #Loop through all states
                    temp_prob_vec[m] = Γ[latest_state, m] * (pdf(Normal(latest_μ[m], sqrt(latest_σ2[m])), raw_data[j])) * back_mat[1][j, m]
                end
                #Normalize the temp_prob_vec
                temp_prob_vec = temp_prob_vec / sum(temp_prob_vec)
                #Update MC chain
                MC_chain[j] = rand(Categorical(temp_prob_vec))
            end

            #Compute the complete log-likelihood
            normalized_log_likelihood = transpose(initial_dist) * state_dep_diag(n_states, raw_data[1], latest_μ, latest_σ2) * back_mat[1][1, :]
            unnormalized_log_likelihood = log(normalized_log_likelihood) + back_mat[2]
            append!(likelihood_vec, unnormalized_log_likelihood)
            #Compute the latest time the process visits the calibration period

            for l in 1:(length(MC_chain)-1)
                if MC_chain[l+1] != MC_chain[l]
                    latest_change = l
                end
            end
            #Append the latest time point for a change to state_register
            append!(state_register, latest_change)

            #Add the number of transitions to the Γ_output matrix
            if i > burn_in
                Γ_output = Γ_output + transition_counter_mat
            end

        end

        #Compute the LPS
        latest_state = MC_chain[length(MC_chain)]
        println("Latest state: ", latest_state)
        TPM = Γ[latest_state, :]
        println("TPM: ", TPM)
        sampled_state = rand(Categorical(TPM))
        println("Sampled state: ", sampled_state)
        # compute the LPS
        LPS = LPS + log(
            pdf(
                Normal(
                    μ_post_draws[n_iter, sampled_state],
                    sqrt(σ2_post_draws[n_iter, sampled_state])
                ),
                raw_data[n_obs]
            )
        )

        #Normalize the posterior transition probabilities
        Γ_output = convert(Array{Float64}, Γ_output)
        for l = 1:n_states
            Γ_output[l, :] = Γ_output[l, :] / sum(Γ_output[l, :])
        end
        if global_counter == lps_iterations
            LPS = LPS / lps_iterations
            return μ_post_draws, σ2_post_draws, MC_chain, Γ, state_register, Γ_output, LPS
        end
    end
end

# if (i <= n_iter)
#     n_obs = length(raw_data[1:(length(raw_data)-lps_iterations)])
# else
#     iteration_greater_than_n_iter = i - n_iter
#     n_obs = length(raw_data[1:(n_iter+iteration_greater_than_n_iter)])
#     #Simulate Hidden state
#     latest_state = MC_chain[length(MC_chain)]
#     TPM = Γ[latest_state, :]
#     sampled_state = rand(Categorical(temp_prob_vec))
#     append!(MC_chain, sampled_state)
#     #compute the LPS
#     LPS = LPS + log(
#         pdf(
#             Normal(
#                 μ_post_draws[i-1, sampled_state],
#                 sqrt(σ2_post_draws[i-1, sampled_state])
#             ),
#             raw_data[n_obs]
#         )
#     )
# end