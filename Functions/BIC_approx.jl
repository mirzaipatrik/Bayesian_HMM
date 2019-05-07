function BIC_approx(log_likelihood, n_parameters, n_obs)
    return -2*log_likelihood+n_parameters*log(n_obs)
end
