
function proposal_jitter_multi(
    rng, theta::AbstractVector{Bool}, ix_subject::Int, n_t_steps::Int
)
    n_multi = 2

    ix_start = (ix_subject - 1) * n_t_steps + 1
    ix_end = ix_start + n_t_steps - 1
    theta_view = @view theta[ix_start:ix_end]

    r_sample = rand(rng)

    if r_sample < 1/3
        # Case 1 - remove an infection
        inf_indices = findall(theta_view)

        if length(inf_indices) >= n_multi
            ix_remove = sample(rng, inf_indices, n_multi)

            # Forward transition probability
            # 1/3 * 1 / n_I(s)
            p_forward = 1/3 * 1 / length(inf_indices)

            # Reverse transition probability, from case 2
            # 1/3 * 1 / n_U(s+1) = 1/3 * 1 / (n_t_steps - n_I(s) + 1)
            p_reverse = 1/3 * 1 / (n_t_steps - length(inf_indices) + n_multi)

            log_hastings_ratio = log(p_reverse) - log(p_forward)

            return SVector{n_multi}(ix_remove), log_hastings_ratio
        end
    elseif r_sample < 2/3
        # Case 2 - add an infection
        not_inf_indices = findall(.!theta_view)

        if length(not_inf_indices) >= n_multi
            ix_add = sample(rng, not_inf_indices, n_multi)

            # Forward transition probability
            # 1/3 * 1 / n_U(s)
            p_forward = 1/3 * 1 / length(not_inf_indices)

            # Reverse transition probability, from case 1
            # 1/3 * 1 / n_I(s+1) = 1/3 * 1 / (n_t_steps - n_U(s) + 1)
            p_reverse = 1/3 * 1 / (n_t_steps - length(not_inf_indices) + n_multi)

            log_hastings_ratio = log(p_reverse) - log(p_forward)

            return SVector{n_multi}(ix_add), log_hastings_ratio
        end
    else
        # Case 3 - move an infection
        inf_indices = findall(theta_view)

        if length(inf_indices) > 0
            ix_t_from = sample(rng, inf_indices)
            ix_t_to = ix_t_from + (rand(rng, Bool) ? -1 : 1)

            if ix_t_to > 0 && ix_t_to <= n_t_steps && !theta_view[ix_t_to]
                # Forward and reverse transition probabilities are equal
                log_hastings_ratio = 0.0
                
                return SA[ix_t_from, ix_t_to], log_hastings_ratio
            end
        end
    end


    return SVector{0, Int}(), 0.0 # Nothing has occurred, no hastings ratio
end