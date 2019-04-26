#=
This function computes the backward probabilities and summerizes
them into a matrix, where the each row represent the backward probabilities
β_k, k ∈ {1,2,...n}. Throughout, we will normalize the
backward probabilities. See e.g. MacDonald and Zucchini (2009) or
Cappé et al. (2005) for more information.
=#

function backward_function(n_states, raw_data, Γ, μ_draws, σ2_draws)
    n = length(raw_data)
    back_mat = zeros(n, n_states)  #n_rows:= rows, n_states:=columns

    #Initialize the nth backward probability corresponding to a vector of 1's
    back_mat[n,:] = fill(1.0, n_states)

    #Normalize the nth backward probability
    back_mat[n,:] = back_mat[n,:] / sum(back_mat[n,:])

    #=
    Now the fun begins... Let's loop through the raw data and compute the
    backward probs for k = n-1, n-2,...,1. Later, the backward probs are
    normalized. Otherwise, the probs will tend to zero or infininty
    exponentially fast depending on the limit of summation.
    See e.g. Cappé et al. (2005) for a thorough discussion on this matter.
    =#
    for k = (n-1):-1:1
        #Compute the kth backward probabilities
        back_mat[k,:] = Γ*state_dep_diag(n_states, raw_data[k+1], μ_draws, σ2_draws)*back_mat[k+1,:]
        #Normalize tha backward probs
        back_mat[k,:] = back_mat[k,:]/sum(back_mat[k,:])
    end
    return back_mat

end
