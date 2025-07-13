
function pointwise_likelihood(
    chain_df, model_data_R, 
    turing_model_name, infection_prior;
    fixed_params = missing,
    use_corrected_titre = true
)
    model_data = convert_model_data(model_data_R)

    obs_df = DataFrame(model_data["observations"])

   sp = read_static_parameters(model_data)

    prior_infection_dist = select_infection_prior(infection_prior, sp)

    turing_model = select_turing_model(turing_model_name)

    model = make_waning_model(
       sp, obs_df;
        prior_infection_dist = prior_infection_dist,
        use_corrected_titre = use_corrected_titre,
        turing_model = turing_model
    );


    if !isnothing(fixed_params)
        fixed_params_tuple = NamedTuple(fixed_params)
        model = model | fixed_params_tuple
    else
        fixed_params_tuple = missing
    end

    return model_pointwise_likelihood(
        chain_df,
       sp, obs_df,
        model;
        fixed_params = fixed_params_tuple
    )
end