
include("../dependencies.jl")


function step!(x, rng)
    r = rand()

    if r < 1/3
        true_indices = findall(x)

        if length(true_indices) > 0
            ix_true = sample(rng, true_indices)

            x[ix_true] = false
        end
    elseif r < 2/3
        false_indices = findall(.!x)

        if length(false_indices) > 0
            ix_false = sample(rng, false_indices)

            x[ix_false] = true

        end
    else
        true_indices = findall(x)
        false_indices = findall(.!x)

        if length(true_indices) > 0 && length(false_indices) > 0
            ix_true = sample(rng, true_indices)
            ix_false = sample(rng, false_indices)

            x[ix_true] = false
            x[ix_false] = true


        end
    end
end


rng = Random.Xoshiro(1)

t = 1:100
n_samples = 10000

N_x = zeros(Int, length(t), n_samples)

for j in 1:n_samples

    x = rand(Bernoulli(0.5), 100)
    for i in eachindex(t)
        for k in 1:100
            step!(x, rng)
        end
        N_x[i, j] = sum(x)
    end
end


plot(N_x[:, 1:32], lc = "black", alpha = 0.1)

histogram(vec(N_x[end, :]), bins = 0:10:100)