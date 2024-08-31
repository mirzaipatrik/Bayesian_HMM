using Distributions
using Random
using Statistics

function MCMC_Simulation(Γ, π, μ, σ, len_of_chain)
    #Allocate memory for chain
    states = fill(0, len_of_chain)
    observations = fill(0.0, len_of_chain)

    #Initialize chain
    states[1] = rand(Categorical(π))

    #sample the first observation
    observations[1] = rand(Normal(μ[states[1]], σ[states[1]]))


    for n = 2:len_of_chain
        probs = Categorical(Γ[states[n-1], :])
        states[n] = rand(probs) #Generate next state
        observations[n] = rand(Normal(μ[states[n]], sqrt(σ[states[n]])))
    end
    return observations
end
