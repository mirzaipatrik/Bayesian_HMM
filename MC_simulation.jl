#=
Input: Γ:= Transition probability matrix, and
the length of the MC chain
=#

 using Distributions
 using Random
 using Statistics



function MC_simulation(Γ, len_of_chain)

    #Allocate memory for chain
    states = fill(0, len_of_chain)

    #Initialize chain
    states[1] = 1

    for n = 2:len_of_chain
        probs = Categorical(Γ[states[n-1],:])
        states[n] = rand(probs) #Generate next state
    end
    return(states)
end
