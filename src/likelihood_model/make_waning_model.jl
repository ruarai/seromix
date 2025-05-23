
struct FixedModelParameters
    n_t_steps::Int
    n_subjects::Int

    antigenic_distances::Matrix{Float64}

    time_diff_matrix::Matrix{Float64}

    subject_birth_ix::Vector{Int}
end

function make_waning_model(
    model_parameters::FixedModelParameters,
    obs_df::DataFrame;
    prior_infection_dist::Distribution,
    turing_model = waning_model_kucharski
)
    all(diff(obs_df.ix_subject) .>= 0) || throw(ArgumentError("ix_subject in obs_df must be sorted in ascending order."))

    n_max_ind_obs = maximum(length.(make_obs_views(obs_df)))

    individual_titre_obs = [obs_df.observed_titre[v] for v in make_obs_views(obs_df)]
    return turing_model(
        model_parameters,
        prior_infection_dist,

        make_obs_lookup(obs_df),
        make_obs_views(obs_df),
        n_max_ind_obs,
        individual_titre_obs
    );
end