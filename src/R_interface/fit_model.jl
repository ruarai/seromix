


function fit_model(
    model_data_R;

    infection_prior = (name = "BetaBernoulli", alpha = 1.0, beta = 1.0),
    proposal_name = "corrected",
    initial_params_name = "kucharski_sim_study",
    sampler_name = "default",
    turing_model_name = "kucharski",

    n_samples = 2_000,
    n_thinning = 1,
    n_chain = 4,

    fixed_params = nothing,

    use_corrected_titre = true,
    mixture_importance_sampling = false,

    rng_seed = 1
)
    rng = Random.Xoshiro(rng_seed)

    model_data = convert_model_data(model_data_R)

    obs_df = DataFrame(model_data["observations"])

    sp = read_static_parameters(model_data)

    prior_infection_dist = select_infection_prior(infection_prior, sp)
    proposal_function = select_proposal_function(proposal_name)

    initial_params = select_initial_params(initial_params_name, n_chain, sp, model_data, obs_df, rng)

    turing_model = select_turing_model(turing_model_name)

    model = make_waning_model(
       sp, obs_df;
        prior_infection_dist = prior_infection_dist,
        use_corrected_titre = use_corrected_titre,
        turing_model = turing_model,
        mixture_importance_sampling = mixture_importance_sampling
    );

    if !isnothing(fixed_params)
        fixed_params_tuple = NamedTuple(fixed_params)
        model = model | fixed_params_tuple

        # Remove any now fixed params from the initial_params Vector
        initial_params = [Base.structdiff(p, fixed_params_tuple) for p in initial_params]
    end

    sampler = select_sampler(sampler_name, model, sp, proposal_function)

    chain = sample_chain(
        model, initial_params, sampler, sp, rng;
        n_sample = n_samples, n_thinning = n_thinning, n_chain = n_chain
    );

    set_lp!(model, chain)

    return DataFrame(chain)
end

function select_sampler(sampler_name, model, sp, proposal_function)
    if sampler_name == "default"
        return make_gibbs_sampler(model, sp, proposal_function)
    elseif sampler_name == "original"
        return make_gibbs_sampler_original(model, sp, proposal_function)
    elseif sampler_name == "slice_sampler"
        return make_gibbs_sampler_slice(model, sp, proposal_function)
    end

    error("Invalid sampler specified")
end


function select_infection_prior(infection_prior, p)
    infection_prior = NamedTuple(infection_prior)
    prior_name = infection_prior.name

    if prior_name == "Bernoulli"
        return MatrixBernoulli(infection_prior.p, p)
    elseif prior_name == "BetaBernoulli"
        return MatrixBetaBernoulli(infection_prior.alpha, infection_prior.beta, p)
    elseif prior_name == "BetaBernoulliTimeVarying"
        return MatrixBetaBernoulliTimeVarying(infection_prior.alpha, infection_prior.beta, p)
    elseif prior_name == "BetaBernoulliSubjectVarying"
        return MatrixBetaBernoulliSubjectVarying(infection_prior.alpha, infection_prior.beta, p)
    end

    error("Invalid infection prior specified")
end

function select_proposal_function(proposal_name)
    if proposal_name == "uncorrected"
        return proposal_original_uncorrected
    elseif proposal_name == "corrected"
        return proposal_original_corrected
    elseif proposal_name == "jitter"
        return proposal_jitter
    end

    error("Invalid proposal function specified")
end

function select_initial_params(initial_params_name, n_chain, sp, model_data, obs_df, rng)
    if initial_params_name == "kucharski_sim_study"
        return make_initial_params_kucharski_sim_study(sp, obs_df, n_chain, rng)
    elseif initial_params_name == "kucharski_data_study"
        return make_initial_params_kucharski_data_study(sp, n_chain, model_data["initial_infections_manual"], rng)
    elseif initial_params_name == "kucharski_data_study_fluscape"
        return make_initial_params_kucharski_data_study_fluscape(sp, obs_df, n_chain, rng)
    elseif initial_params_name == "broad"
        return make_initial_params_broad(sp, n_chain, rng)
    elseif initial_params_name == "age_effect"
        return make_initial_params_age(sp, obs_df, n_chain, rng)
    end

    error("Invalid initial params specified")
end

function select_turing_model(turing_model_name)
    if turing_model_name == "kucharski"
        return waning_model_kucharski
    elseif turing_model_name == "age_effect"
        return waning_model_age_effect
    elseif turing_model_name == "nonlinear"
        return waning_model_non_linear
    end
    error("Invalid turing model specified")
end
