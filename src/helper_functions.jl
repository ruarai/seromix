

function save_draws(draws, filename)
    write_table(filename, DataFrame(draws), format = :parquet)
end

function make_obs_views(obs_df)
    n_ind = maximum(obs_df.ix_subject)
    [
        findfirst(obs_df.ix_subject .== ix_subject):findlast(obs_df.ix_subject .== ix_subject)
        for ix_subject in 1:n_ind
    ]
end

function make_obs_lookup(obs_df)
    n_ind = maximum(obs_df.ix_subject)

    obs_lookup = Vector{Dict{Int64, Vector{Tuple{Int64,Int64}}}}(undef, n_ind)

    obs_ix_offsets = [findfirst(obs_df.ix_subject .== ix_subject) for ix_subject in 1:n_ind]

    for ix_subject in 1:n_ind
        obs_lookup[ix_subject] = Dict{Int64, Vector{Tuple{Int64,Int64}}}()
    end

    for ix_obs in 1:nrow(obs_df)
        ix_subject = obs_df[ix_obs, :ix_subject]
        ix_t_obs = obs_df[ix_obs, :ix_t_obs]
        ix_strain = obs_df[ix_obs, :ix_strain]

        key = ix_t_obs
        
        if !haskey(obs_lookup[ix_subject], key)
            obs_lookup[ix_subject][key] = Int[]
        end

        push!(obs_lookup[ix_subject][key], (ix_strain, ix_obs - obs_ix_offsets[ix_subject] + 1))
    end

    return obs_lookup
end

function make_time_diff_matrix(modelled_years)
    n_t_steps = length(modelled_years)

    time_diff_matrix = zeros(n_t_steps, n_t_steps)

    for ix_t in 1:n_t_steps
        for ix_t_obs in 1:n_t_steps
            time_diff_matrix[ix_t_obs, ix_t] = modelled_years[ix_t_obs] - modelled_years[ix_t]
        end
    end

    return time_diff_matrix
end



function read_obs_df_R(obs_df)
    obs_df = DataFrame(obs_df)

    obs_df[!, :ix_subject] = convert.(Int64,obs_df[!, :ix_subject])
    obs_df[!, :ix_t_obs] = convert.(Int64,obs_df[!, :ix_t_obs])
    obs_df[!, :ix_strain] = convert.(Int64,obs_df[!, :ix_strain])

    obs_df[!, :observed_titre] = convert.(Float64,obs_df[!, :observed_titre])

    return obs_df
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


function ntuple_to_matrix(nt, n_t_steps, n_subjects)
    elements = collect(nt)

    result = Array{Float64}(undef, n_t_steps, n_subjects)

    i = 1
    for ix_subject in 1:n_subjects, ix_t in 1:n_t_steps
        result[ix_t, ix_subject] = elements[i].data[1]
        i += 1
    end

    return result
end


function read_model_parameters(dict)
    antigenic_distances = dict["antigenic_distances"]
    modelled_years = dict["modelled_years"]
    subject_birth_ix::Vector{Int64} = DataFrame(dict["subject_birth_data"]).ix_t_birth
    
    n_t_steps = length(modelled_years)
    n_subjects = length(subject_birth_ix)
    
    time_diff_matrix = make_time_diff_matrix(modelled_years)
    
    params = FixedModelParameters(
        n_t_steps, n_subjects,
        antigenic_distances, time_diff_matrix, subject_birth_ix
    )
    
    return params
end

function df_to_tuple(df)
    NamedTuple.(eachrow(df))
end


function mask_infections_birth_year!(infections, subject_birth_ix)
    n_subjects = length(subject_birth_ix)
    for ix_subject in 1:n_subjects
        if subject_birth_ix[ix_subject] > 1
            infections[1:(subject_birth_ix[ix_subject] - 1), ix_subject] .= false
        end
    end
end


function chain_infections_matrix(chain, ix_iter, ix_chain, params)
    infections = zeros(Bool, params.n_t_steps, params.n_subjects)

    for ix_subject in 1:params.n_subjects
        for ix_t in max(1, params.subject_birth_ix[ix_subject]):params.n_t_steps
            infections[ix_t, ix_subject] = chain[ix_iter, Symbol("infections[$ix_t, $ix_subject]"), ix_chain]
        end
    end

    return infections
end



# TODO add birth year masking
function chain_sum_infections(chain, params)
    # Assuming infections is at end
    chain_sub = chain[:,(end - params.n_t_steps * params.n_subjects + 1):end,:]
    sums = sum(Array(chain_sub)[:, size(chain, 3):end], dims = 2)
    return [sums[(1:size(chain, 1)) .+ i * size(chain, 1)] for i in 0:(size(chain, 3) - 1)]
end


function chain_sum_infections(chain, params, iterations, ix_t_subset)
    n_infections = zeros(length(iterations), size(chain)[3])

    Threads.@threads for ix_c in 1:size(chain)[3]
        for ix_i in eachindex(iterations)
            for ix_subject in 1:params.n_subjects
                for ix_t in ix_t_subset
                    if params.subject_birth_ix[ix_subject] < ix_t
                        n_infections[ix_i, ix_c] += chain[iterations[ix_i], :, ix_c][Symbol("infections[$ix_t, $ix_subject]")][1]
                    end
                end
            end
        end
    end

    return n_infections
end



function chain_infections_prob(chain, params)
    infections = zeros(Float64, params.n_t_steps, params.n_subjects)

    for ix_subject in 1:params.n_subjects
        for ix_t in max(1, params.subject_birth_ix[ix_subject]):params.n_t_steps
            infections[ix_t, ix_subject] = mean(chain[Symbol("infections[$ix_t, $ix_subject]")])
        end
    end

    return infections
end


function drop_nothing(vec)
    filter(x -> !isnothing(x), vec)
end

function split_vector_indices(N::Int, M::Int)
    if N <= 0
        error("N must be a positive integer.")
    end
    if M <= 0
        error("M must be a positive integer.")
    end

    if N < M
        error("Function is constrained to work only when N >= M. Received N=$N, M=$M.")
    end

    result_ranges = Vector{UnitRange{Int}}(undef, M)
    current_start_index = 1

    for i = 1:M
        base_chunk_size = N ÷ M
        remainder_chunks = N % M

        chunk_size = base_chunk_size + (i <= remainder_chunks ? 1 : 0)
        
        current_end_index = current_start_index + chunk_size - 1

        result_ranges[i] = current_start_index:current_end_index
        current_start_index = current_end_index + 1
    end

    return result_ranges
end

function list_files(base_dir, regex)
    [
        joinpath(current_path, file_name)  
        for (current_path, _, files_in_dir) in walkdir(base_dir) 
        for file_name in files_in_dir
        if occursin(regex, file_name)
    ]
end