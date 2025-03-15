

function save_draws(draws, filename)
    write_parquet(filename, DataFrame(draws))
end


function make_obs_lookup(obs_df)
    obs_lookup = Dict{Tuple{Int,Int},Vector{Int}}()

    for ix_obs in 1:nrow(obs_df)
        i = obs_df[ix_obs, :i]
        t = obs_df[ix_obs, :t]
    
        key = (i, t)
        
        if !haskey(obs_lookup, key)
            obs_lookup[key] = Int[]
        end
        push!(obs_lookup[key], ix_obs)
    end

    return obs_lookup
end

function make_obs_matrix(obs_df)
    return Matrix{Int}(obs_df[:,1:3])
end