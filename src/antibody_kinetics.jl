function zero_cutoff(x)
    return max(0, x)
end

function titre_component(
    mu_long,
    mu_short, omega,
    sigma_long, sigma_short,
    tau,

    infection_number,

    dist,
    t
)
    seniority = zero_cutoff(1 - tau * (infection_number - 1))

    long_term = mu_long * zero_cutoff(1 - sigma_long * dist)
    short_term = mu_short * zero_cutoff(1 - omega * t) * zero_cutoff(1 - sigma_short * dist)

    return seniority * (long_term + short_term)
end

function waning_curve(
    mu_long,
    mu_short, omega,
    sigma_long, sigma_short,
    tau,


    dist_matrix,
    infections, n_strain, 
    t_steps
)
    y = zeros(
        typeof(mu_long),
        n_strain,
        length(t_steps)
    )

    for ix_inf in eachindex(infections), ix_cross_strain in 1:n_strain

        t_inf, strain_inf = Tuple(infections[ix_inf])
        dist = dist_matrix[strain_inf, ix_cross_strain]

        for ix_t in eachindex(t_steps)
            t = t_steps[ix_t]

            if t >= t_inf
                y[ix_cross_strain, ix_t] += titre_component(
                    mu_long,
                    mu_short, omega,
                    sigma_long, sigma_short,
                    tau,

                    ix_inf,
                
                    dist,
                    
                    t - t_inf
                )
            end
        end
    end

    return y
end