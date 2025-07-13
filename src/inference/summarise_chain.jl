

function summarise_chain(chain_file, iterations_drop, model_data)
    chain = read_parquet(DataFrame, chain_file) # TODO make faster :(

    p = read_static_parameters(model_data)

    chain_name = match(r"chain_(.+)(?=\.parquet)", chain_file).match
    run_ix = match(r"\d*(?=/chain)", chain_file).match

    n_chains = length(unique(chain.chain))

    chain_filt = filter(:iteration => i -> i > iterations_drop, chain)
    chain_filt_grouped = groupby(chain_filt, :chain)


    variables = ["mu_long", "mu_short", "omega", "sigma_long", "sigma_short", "tau", "obs_sd"]


    chain_summaries = NamedTuple[]

    for v in variables, c in 1:n_chains
        values = chain_filt_grouped[c][:,v]

        values_summary = (
            run_ix = run_ix,
            chain_name = chain_name,
            variable = v,
            chain = c,
            mean = mean(values),
            sd = std(values),
            q90_lower = quantile(values, 0.05), q90_upper = quantile(values, 0.95),
            q95_lower = quantile(values, 0.025), q95_upper = quantile(values, 0.975)
        )

        push!(chain_summaries, values_summary)
    end


    ix_infections = findall(s -> startswith(s, "infections"), names(chain))

    n_draws = nrow(chain_filt_grouped[1])

    for c in 1:n_chains
        infections = reshape(Matrix(chain_filt_grouped[c][:,ix_infections]), n_draws, p.n_t_steps, p.n_subjects)
        # [mask_infections_birth_year!(view(infections, i, :, :), p.subject_birth_ix) for i in 1:n_draws]

        values = [sum(view(infections, i, :, :)) for i in 1:n_draws]

        values_summary = (
            run_ix = run_ix,
            chain_name = chain_name,
            variable = "total_infections",
            chain = c,
            mean = mean(values),
            sd = std(values),
            q90_lower = quantile(values, 0.05), q90_upper = quantile(values, 0.95),
            q95_lower = quantile(values, 0.025), q95_upper = quantile(values, 0.975)
        )

        push!(chain_summaries, values_summary)

    end

    diagnostics = NamedTuple[]

    for v in variables
        values_matrix = Matrix{Float64}(
            unstack(chain_filt[:, [:iteration, :chain, Symbol(v)]], :chain, Symbol(v))[:,Not(:iteration)]
        )

        diagnostics_summary = (
            run_ix = run_ix,
            chain_name = chain_name,
            variable = v,
            ess = ess(values_matrix),
            ess_tail = ess(values_matrix; kind = :tail),
            rhat = rhat(values_matrix)
        )

        push!(diagnostics, diagnostics_summary)
    end

    return (
        chain_summaries = DataFrame(chain_summaries),
        diagnostics = DataFrame(diagnostics)
    )
end