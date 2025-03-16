

function save_draws(draws, filename)
    write_parquet(filename, DataFrame(draws))
end

function make_obs_views(obs_df)
    n_ind = maximum(complete_obs.i)
    [
        findfirst(obs_df.i .== i):findlast(obs_df.i .== i)
        for i in 1:n_ind
    ]
end

function make_obs_lookup(obs_df)
    n_ind = maximum(complete_obs.i)

    obs_lookup = Vector{Dict{Int64, Vector{Tuple{Int64,Int64}}}}(undef, n_ind)

    for i in 1:n_ind
        obs_lookup[i] = Dict{Int64, Vector{Tuple{Int64,Int64}}}()
    end

    for ix_obs in 1:nrow(obs_df)
        i = obs_df[ix_obs, :i]
        t_obs = obs_df[ix_obs, :t]
        s_obs = obs_df[ix_obs, :s]

        key = t_obs
        
        if !haskey(obs_lookup[i], key)
            obs_lookup[i][key] = Int[]
        end

        push!(obs_lookup[i][key], (s_obs, ix_obs))
    end

    return obs_lookup
end


function make_obs_lookup_individuals(obs_df)
    n_ind = maximum(complete_obs.i)

    obs_lookup = Vector{Dict{Int64, Vector{Tuple{Int64,Int64}}}}(undef, n_ind)

    obs_ix_offsets = [findfirst(complete_obs.i .== i) for i in 1:n_ind]

    for i in 1:n_ind
        obs_lookup[i] = Dict{Int64, Vector{Tuple{Int64,Int64}}}()
    end

    for ix_obs in 1:nrow(obs_df)
        i = obs_df[ix_obs, :i]
        t_obs = obs_df[ix_obs, :t]
        s_obs = obs_df[ix_obs, :s]

        key = t_obs
        
        if !haskey(obs_lookup[i], key)
            obs_lookup[i][key] = Int[]
        end

        push!(obs_lookup[i][key], (s_obs, ix_obs - obs_ix_offsets[i] + 1))
    end

    return obs_lookup
end

function log_callback(rng, model, sampler, sample, state, iteration; kwargs...)
    if iteration % 50 == 0
        print("$iteration,")
    end
end


# Expand grid via
# https://stackoverflow.com/a/67733908
function expand_grid(; iters...)
    var_names = collect(keys(iters))

    var_itr = [1:length(x) for x in values(iters)]
    var_ix = vcat([collect(x)' for x in Iterators.product(var_itr...)]...)
    out = DataFrame()
    for i = eachindex(var_names)
        out[:,var_names[i]] = collect(iters[i])[var_ix[:,i]]
    end
    return out
end