include("dependencies.jl")


model_data = load("data/model_data.hdf5")

obs_df = read_obs_df_R(model_data["observations_df"])

antigenic_distances = model_data["antigenic_distances"]
modelled_years = model_data["modelled_years"]

subject_birth_ix::Vector{Int64} = model_data["subject_birth_ix"]

n_ind = maximum(obs_df.ix_subject)
n_strain = size(antigenic_distances, 1)
n_t_steps = length(modelled_years)

time_diff_matrix = make_time_diff_matrix(modelled_years)

model = waning_model(
    antigenic_distances, time_diff_matrix,
    subject_birth_ix,

    n_ind, n_t_steps,
    
    make_obs_lookup(obs_df),
    make_obs_views(obs_df), nrow(obs_df),
    obs_df.observed_titre
);

symbols_not_inf = model_symbols_apart_from(model, :infections)

# Must somehow balance the level of exploration of the MH sampler
# with that of the HMC sampler -- so repeating MH or changing HMC step size
gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(),
    symbols_not_inf => HMC(0.002, 10) # Must be reduced with number of individuals?
)

chain = @time sample(model, gibbs_sampler, 100);

chain = @time sample(
    model, gibbs_sampler, 
    MCMCThreads(), 6000, 6,
    callback = log_callback
);

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

save_draws(chain, "data/2018/chain.parquet")

