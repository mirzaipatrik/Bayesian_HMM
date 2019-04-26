#=
This function groupes the observations according to the latent states
in the Markov chain. Hence, all observations that belong to state k,
k ∈ {1,2,...,m} are gouped into one array
=#

function grouped_obs(n_states, MC_chain, raw_data)
    #Create an array for storing arrays of observations
    list_of_vecs = []

    #Loop through all states
    for m = 1:n_states
        temp_vec = []

        #Loop through all observations
        for j = 1:length(MC_chain)
            if MC_chain[j] == m
                append!(temp_vec, raw_data[j])
            end
        end
        push!(list_of_vecs, temp_vec)
    end
    return list_of_vecs
end
