#=
This function counts the number of transitions from state
i to state j, ∀ i,j ∈{1,2,...,m}, where m denotes the total number of states,
m ∈ Z^+
=#

function transition_counter(X, n_states)
    p = zeros(n_states, n_states)  #Create an empty matrix
    for t = 1:(length(X) - 1)
        p[X[t], X[t + 1]] = p[X[t], X[t + 1]] + 1
    end
    return p
end
