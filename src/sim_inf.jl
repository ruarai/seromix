

include("dependencies.jl")

model_data = load("data/model_data.hdf5")


antigenic_distances = model_data["antigenic_distances"]
modelled_years = model_data["modelled_years"]
subject_birth_ix::Vector{Int64} = model_data["subject_birth_ix"]

time_diff_matrix = make_time_diff_matrix(modelled_years)

n_t_steps = length(modelled_years)
n_subjects = length(subject_birth_ix)

mu_long = 2.0
mu_short = 2.7
omega = 0.8
sigma_long = 0.13
sigma_short = 0.03

tau = 0.04


infections = rand(Bernoulli(0.2), (n_t_steps, n_subjects))

for ix_subject in 1:n_subjects
    if subject_birth_ix[ix_subject] > 0
        infections[1:subject_birth_ix[ix_subject], ix_subject] .= false
    end
end

heatmap(infections)


infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
write_parquet("data/infections_df.parquet", infections_df)

complete_obs = expand_grid(
    ix_t_obs = 1:n_t_steps, ix_strain = 1:n_t_steps, ix_subject = 1:n_subjects,
    observed_titre = 0.0
)


waning_curve!(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,

    antigenic_distances, time_diff_matrix,
    subject_birth_ix,
    infections,

    make_obs_lookup(complete_obs),
    make_obs_views(complete_obs),
    complete_obs.observed_titre
)

write_parquet("data/complete_obs.parquet", complete_obs)


obs_df = filter(:ix_t_obs => ix_t_obs -> ix_t_obs % 1 == 0, complete_obs)
obs_df.observed_titre = round.(max.(0, obs_df.observed_titre .+ rand(Normal(0, 1.29), nrow(obs_df))))
write_parquet("data/obs.parquet", obs_df)

model = waning_model(
    antigenic_distances, time_diff_matrix,
    subject_birth_ix,

    n_subjects, n_t_steps,
    
    make_obs_lookup(obs_df),
    make_obs_views(obs_df), nrow(obs_df),
    obs_df.observed_titre
);

symbols_not_inf = model_symbols_apart_from(model, :infections)

# HMC step size must be tuned to balance efficiency of the two samplers
gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(),
    symbols_not_inf => HMC(0.002, 10) # Must be reduced with number of individuals?
)

# chain = @time sample(model, gibbs_sampler, 100);

chain = @time sample(
    model, gibbs_sampler, 
    MCMCThreads(), 5000, 6,
    callback = log_callback
);

# plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
# plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
# plot(chain, [:tau], seriestype = :traceplot)

save_draws(chain, "data/chain.parquet")




ppd_start_draw_ix = 00
n_ppd_samples = 100
ppd_subjects = 1:40

draw_dfs = Array{DataFrame}(undef, n_ppd_samples)

for ix_draw in 1:n_ppd_samples
    chain_sample = sample(chain[ppd_start_draw_ix:length(chain)], 1)

    draw_dfs[ix_draw] = expand_grid(
        ix_t_obs = 1:n_t_steps, 
        ix_strain = 1:n_t_steps, 
        ix_subject = ppd_subjects,
        draw = ix_draw,
        observed_titre = 0.0
    )

    obs_lookup = make_obs_lookup(draw_dfs[ix_draw])
    obs_views = make_obs_views(draw_dfs[ix_draw])


    n_obs_subset = length(obs_views[1])

    p = get_params(chain_sample)
    
    pred_inf = ntuple_to_matrix(p.infections, n_t_steps, n_subjects)


    for ix_subject in ppd_subjects

        y_pred = zeros(n_obs_subset)

        waning_curve_individual!(
            p.mu_long[1], p.mu_sum[1] - p.mu_long[1], 0.8,
            p.sigma_long[1], p.sigma_short[1], p.tau[1],

            antigenic_distances, time_diff_matrix,
            subject_birth_ix[ix_subject],
            AbstractArray{Bool}(view(pred_inf, :, ix_subject)),

            obs_lookup[ix_subject], n_t_steps,

            y_pred
        )

        for (i, ix_obs) in enumerate(obs_views[ix_subject])
            draw_dfs[ix_draw].observed_titre[ix_obs] = y_pred[i]
        end
    end
end

ppd_obs = vcat(draw_dfs...)
write_parquet("data/ppd_obs.parquet", ppd_obs)