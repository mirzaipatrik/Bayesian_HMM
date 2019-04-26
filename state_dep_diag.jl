#=
This function computes P(y_k), i.e. diagonal matrix with diagonal elements
corresponding to the emission probabilities f(y_k|μ_j, σ_j),
j ∈ {1,2,...,m} where m denotes the number of states
=#

function state_dep_diag(n_states, raw_obs, μ_draws, σ2_draws)
  #Create diag-matrix
  diag_mat = Diagonal(fill(1, n_states))
  #Just a simple conversion to Float64 types
  diag_mat = convert(Array{Float64}, diag_mat)
  for m=1:n_states
    diag_mat[m, m] = (pdf(Normal(μ_draws[m], sqrt(σ2_draws[m])), raw_obs))
  end
  return diag_mat
end
