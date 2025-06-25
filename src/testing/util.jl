
function make_test_data(rng, n_t_steps, n_subjects)
    modelled_years = collect(1:n_t_steps)

    time_diff_matrix = make_time_diff_matrix(modelled_years)
    antigenic_distance = abs.(time_diff_matrix)

    subject_birth_data = DataFrame(
        ix_subject = 1:n_subjects, 
        ix_t_birth = floor.(Int, reverse(1:n_subjects) .* 0.7)
    )

    p = FixedModelParameters(
        n_t_steps, n_subjects,
        antigenic_distance,
        time_diff_matrix,
        subject_birth_data.ix_t_birth
    )

    model_params = (
        mu_long = 2.0,
        mu_short = 2.0,
        omega = 0.75,
        sigma_long = 0.15,
        sigma_short = 0.05,
        tau = 0.05,
        obs_sd = 1.5
    )

    infections = rand(rng, Bernoulli(0.2), (n_t_steps, n_subjects))
    mask_infections_birth_year!(infections, p.subject_birth_ix)

    return (
        model_params = model_params,
        p = p,
        infections = infections
    )
end