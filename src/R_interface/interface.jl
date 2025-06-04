

include("../dependencies.jl")


function fit_model(
    run_dir;

    infection_prior = (name = "MatrixBetaBernoulli", alpha = 1.0, beta = 1.0),
    proposal_name = "original_corrected",
    initial_params_name = "kucharski_sim_study",

    n_samples = 2_000,
    n_thinning = 1,
    n_chain = 4,

    fixed_params = nothing,

    rng_seed = 1
)
    rng = Random.Xoshiro(rng_seed)

    model_data = load("$run_dir/model_data.hdf5")

    obs_df = DataFrame(model_data["observations"])

    p = read_model_parameters(model_data)

    prior_infection_dist = select_infection_prior(infection_prior, p)
    proposal_function = select_proposal_function(proposal_name)

    initial_params = select_initial_params(initial_params_name, n_chain, p, model_data, obs_df, rng)

    model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);

    if !isnothing(fixed_params)
        fixed_params_tuple = NamedTuple(fixed_params)
        model = model | fixed_params_tuple

        # Remove any now fixed params from the initial_params Vector
        initial_params = [Base.structdiff(p, fixed_params_tuple) for p in initial_params]
    end


    gibbs_sampler = make_gibbs_sampler(model, p, proposal_function)

    chain = sample_chain(
        model, initial_params, gibbs_sampler, rng;
        n_sample = n_samples, n_thinning = n_thinning, n_chain = n_chain
    );

    return DataFrame(chain)
end


function select_infection_prior(infection_prior, p)
    infection_prior = NamedTuple(infection_prior)
    prior_name = infection_prior.name

    if prior_name == "MatrixBernoulli"
        return MatrixBernoulli(infection_prior.p, p.n_t_steps, p.n_subjects)
    elseif prior_name == "MatrixBetaBernoulli"
        return MatrixBetaBernoulli(infection_prior.alpha, infection_prior.beta, p.n_t_steps, p.n_subjects)
    end

    error("Incorrect infection prior specified")
end

function select_proposal_function(proposal_name)
    if proposal_name == "original_uncorrected"
        return proposal_original_uncorrected
    elseif proposal_name == "original_corrected"
        return proposal_original_corrected
    elseif proposal_name == "jitter"
        return proposal_jitter
    end

    error("Incorrect proposal function specified")
end

function select_initial_params(initial_params_name, n_chain, p, model_data, obs_df, rng)
    if initial_params_name == "kucharski_sim_study"
        return make_initial_params_kucharski_sim_study(p, obs_df, n_chain, rng)
    elseif initial_params_name == "kucharski_data_study"
        return make_initial_params_kucharski_data_study(n_chain, model_data["initial_infections_manual"], rng)
    elseif initial_params_name == "kucharski_data_study_fluscape"
        return make_initial_params_kucharski_data_study_fluscape(p, obs_df, n_chain, rng)
    end


    error("Incorrect initial params specified")
end