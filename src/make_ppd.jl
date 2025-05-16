


include("dependencies.jl")


data_code = "sim_study_simple_1"

run_dir = "runs/$(data_code)/"


model_data = load("$run_dir/model_data.hdf5")

p = read_model_parameters(model_data)


chain = QuackIO.read_parquet(DataFrame, "$run_dir/chain.parquet")

# draw_n(i, c) = (i .- minimum(c)) .+ (c .- 1) .* (maximum(i) - minimum(i) + 1)
# chain.draw .= draw_n(chain.iteration, chain.chain)


chain_filt = filter(:iteration => iteration -> iteration > 2500, chain)
n_ppd_subjects = 1
n_draws = 100

obs_df = expand_grid(
    ix_t_obs = 1:p.n_t_steps, 
    ix_strain = 1:p.n_t_steps, 
    ix_subject = 1:n_ppd_subjects,
    observed_titre = 0.0,
    ix_draw = 1:n_draws,
)

obs_df_grouped = groupby(obs_df, :ix_draw)

# Share across draws
infections = Matrix{Bool}(undef, p.n_t_steps, p.n_subjects)

for ix_draw in obs_df.ix_draw
    ix_sample = sample(1:nrow(chain_filt))

    draw = chain_filt[ix_sample,:]
    
    obs_lookup = make_obs_lookup(obs_df_grouped[ix_draw])

    for ix_subject in 1:p.n_subjects, ix_t in 1:p.n_t_steps
        infections[ix_t, ix_subject] = chain_filt[ix_sample, "infections[$ix_t, $ix_subject]"]
    end

    for ix_subject in 1:n_ppd_subjects
        waning_curve_individual!(
            draw.mu_long, draw.mu_sum - draw.mu_long, 0.75,
            draw.sigma_long, draw.sigma_short, draw.tau,

            p.antigenic_distances, p.time_diff_matrix,
            p.subject_birth_ix[ix_subject],

            AbstractArray{Bool}(view(infections, :, ix_subject)),

            obs_lookup[ix_subject], 
            obs_df_grouped[ix_draw].observed_titre
        )
    end
end

ppd_df = DataFrame(obs_df_grouped)


save_draws(ppd_df, "$run_dir/ppd.parquet")