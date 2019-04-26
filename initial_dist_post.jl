#=
This function computes the posterior parameters that are used
to sample parameters for the initial distribution, π.
Recall that π|X ∼ Dir (I{X_(1)=1}+1, I{X_(1)=2}+1,...,I{X_(1)=m}+1).
Thus, we want to add 1 where the initial state is present
=#

function initial_dist_post(n_states, MC_chain)
  #Almost our posterior parameters for the initial dist π
  post_parameters = fill(1, n_states)

  for m=1:n_states
    if MC_chain[1] == m
      post_parameters[m] = 2
    end
  end
  return post_parameters


end
