


function titre_component(
    mu, omega,
    t
)
    return mu  * max(0, 1 - omega * t)
end

function waning_curve(
    mu, omega,
    t_infected,

    t_steps
)
    y = zeros(typeof(mu), length(t_steps))

    for ix_inf in eachindex(t_infected)
        for ix_t in eachindex(t_steps)
            t = t_steps[ix_t]
            t_inf = t_infected[ix_inf]

            if t >= t_inf
                y[ix_t] += titre_component(mu, omega, t - t_inf)
            end
        end
    end

    return y
end