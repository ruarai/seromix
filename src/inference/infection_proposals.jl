function proposal_original_corrected(
    rng, theta::AbstractVector{Bool}, ix_subject::Int, n_t_steps::Int, subject_birth_ix::Int
)
    ix_start = (ix_subject - 1) * n_t_steps + max(1, subject_birth_ix)
    ix_end = ix_subject * n_t_steps
    theta_view = @view theta[ix_start:ix_end]

    t_alive = n_t_steps - max(1, subject_birth_ix) + 1

    r_sample = rand(rng)

    if r_sample < 1/3
        # Case 1 - remove an infection

        inf_indices = findall(theta_view)

        if length(inf_indices) > 0
            ix_remove = sample(rng, inf_indices)

            # Forward transition probability
            # 1/3 * 1 / n_I(s)
            p_forward = 1/3 * 1 / length(inf_indices)

            # Reverse transition probability, from case 2
            # 1/3 * 1 / n_U(s+1) = 1/3 * 1 / (t_alive - n_I(s) + 1)
            p_reverse = 1/3 * 1 / (t_alive - length(inf_indices) + 1)

            log_hastings_ratio = log(p_reverse) - log(p_forward)

            return SA[ix_remove], log_hastings_ratio
        end
    elseif r_sample < 2/3
        # Case 2 - add an infection

        not_inf_indices = findall(.!theta_view)

        if length(not_inf_indices) > 0
            ix_remove = sample(rng, not_inf_indices)

            # Forward transition probability
            # 1/3 * 1 / n_U(s)
            p_forward = 1/3 * 1 / length(not_inf_indices)

            # Reverse transition probability, from case 1
            # 1/3 * 1 / n_I(s+1) = 1/3 * 1 / (t_alive - n_U(s) + 1)
            p_reverse = 1/3 * 1 / (t_alive - length(not_inf_indices) + 1)

            log_hastings_ratio = log(p_reverse) - log(p_forward)

            return SA[ix_remove], log_hastings_ratio
        end
    else
        # Case 3 - move an infection

        inf_indices = findall(theta_view)
        not_inf_indices = findall(.!theta_view)

        if length(inf_indices) > 0 && length(not_inf_indices) > 0
            ix_t_from = sample(rng, inf_indices)
            ix_t_to = sample(rng, not_inf_indices)

            # Forward and reverse transition probabilities are equal
            log_hastings_ratio = 0.0
            
            return SA[ix_t_from, ix_t_to], log_hastings_ratio
        end
    end

    return SVector{0, Int}(), 0.0 # Nothing has occurred, no hastings ratio
end

function proposal_original_uncorrected(
    rng, theta::AbstractVector{Bool}, ix_subject::Int, n_t_steps::Int, subject_birth_ix::Int
)
    ix_start = (ix_subject - 1) * n_t_steps + max(1, subject_birth_ix)
    ix_end = ix_subject * n_t_steps
    theta_view = @view theta[ix_start:ix_end]

    r_sample = rand(rng)

    if r_sample < 1/3
        # Case 1 - remove an infection

        inf_indices = findall(theta_view)

        if length(inf_indices) > 0
            ix_remove = sample(rng, inf_indices)
            
            return SA[ix_remove], 0.0
        end
    elseif r_sample < 2/3
        # Case 2 - add an infection

        not_inf_indices = findall(.!theta_view)

        if length(not_inf_indices) > 0
            ix_remove = sample(rng, not_inf_indices)

            return SA[ix_remove], 0.0
        end
    else
        # Case 3 - move an infection

        inf_indices = findall(theta_view)
        not_inf_indices = findall(.!theta_view)

        if length(inf_indices) > 0 && length(not_inf_indices) > 0
            ix_t_from = sample(rng, inf_indices)
            ix_t_to = sample(rng, not_inf_indices)

            return SA[ix_t_from, ix_t_to], 0.0
        end
    end

    return SVector{0, Int}(), 0.0 # Nothing has occurred, no hastings ratio
end