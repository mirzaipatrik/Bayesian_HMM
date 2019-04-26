function state_counts(n_states, MC_Chain)
    #=
    This function counts the number of visits to state k, k ∈ {1,2,...,m}
    where m denotes the total number of states
    =#

    #Allocate memory
    output_vector = fill(0, n_states)

    for k = 1:n_states
        output_vector[k] = sum(MC_Chain .== k)
    end

    return output_vector

end
