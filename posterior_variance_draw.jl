#=
This function draws random variables from the
scaled-inverse χ^2 distribution, where υ_n and σ2_n
are the "updated posterior parameters"
=#

function posterior_variance_draw(υ, σ)
    σ2_vec = []
    n_draws = length(υ)
    for k=1:n_draws
        χ2 = rand(Chisq(υ[k]))  #Draw from χ2 distirbution
        σ2 = υ[k]*σ[k]/χ2  #This is our posterior variance draw
        append!(σ2_vec, σ2)
    end
    return σ2_vec


end
