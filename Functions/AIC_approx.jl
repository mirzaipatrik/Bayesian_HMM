function AIC_approx(log_likelihood, n_parameters)
    return -2*log_likelihood+2*n_parameters
end
