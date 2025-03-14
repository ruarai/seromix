

function save_draws(draws, filename)
    write_parquet(filename, DataFrame(draws))
end


function obs_df_to_matrix(obs_df)
    return Matrix(Matrix(obs_df[:,1:3])')
end