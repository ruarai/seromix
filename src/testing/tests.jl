

include("../dependencies.jl")
include("util.jl")

using Test, StableRNGs, BenchmarkTools
rng = Random.Xoshiro(1)


@testset "Titre-normal distribution" begin
    y_latent = rand(rng, Normal(1.0, 1.5), 100)
    σ = 1.5
    x = TitreArrayNormal(y_latent, σ, const_titre_min, const_titre_max)

    ys = [rand(rng, x) for i in 1:10000]
    y1 = ys[1]

    # Are random values in the support?
    @test all(y1 .>= const_titre_min) 
    @test all(y1 .<= const_titre_max) 

    # Test that SIMD implementation is approx equal to usual approach
    l1 = [apply_logpdf(y, y_latent, σ, const_titre_min, const_titre_max) for y in ys]
    l2 = [apply_logpdf_simd(y, y_latent, σ, const_titre_min, const_titre_max) for y in ys]

    @test l1 ≈ l2

    # Is pdf over sum approx 1?
    y_support = const_titre_min:1.0:const_titre_max
    @test sum(exp.([apply_logpdf(y, [1.5], σ, const_titre_min, const_titre_max) for y in y_support])) ≈ 1.0
end


@testset "Waning model" begin
    n_t_steps = 40
    n_subjects = 40
    model_params, sp, infections = make_test_data(StableRNG(1), n_t_steps, n_subjects)

    complete_obs = expand_grid(
        ix_t_obs = 1:n_t_steps, ix_strain = 1:n_t_steps, ix_subject = 1:n_subjects,
        observed_titre = 0.0
    )

    model_cache = WaningModelCache(complete_obs)

    waning_curve!(
        model_params,
        individual_waning_kucharski!,
        sp,
        infections,
        model_cache,
        complete_obs.observed_titre
    )

    # Basic sanity checks
    @test all(complete_obs.observed_titre .>= 0)

    # Lazy way to store the approx. expected result
    @test sum(complete_obs.observed_titre) ≈ 60025.987499999996 


    b_trial = @benchmark waning_curve!(
        $model_params,
        $individual_waning_kucharski!,
        $sp,
        $infections,
        $model_cache,
        $complete_obs.observed_titre
    )

    # Should take less than 1ms, and memory use should be minimal
    @test (median(b_trial).time / 1000) < 1_000
    @test (median(b_trial).memory) < 200
    @test (median(b_trial).allocs) < 10
end

@testset "Inference model" begin
    model_data = load("runs/hanam_2018/model_data.hdf5")
    obs_df = DataFrame(model_data["observations"])
    sp = read_static_parameters(model_data)

    prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, sp)

    initial_params = make_initial_params_broad(sp, 4, rng)

    model = make_waning_model(sp, obs_df; prior_infection_dist = prior_infection_dist);

    b_trial = @benchmark logjoint($model, x) setup=(x=rand(model))

    @test (median(b_trial).time / 1000) < 5_000
end


model_data = load("runs/hanam_2018/model_data.hdf5")
sp = read_static_parameters(model_data)
model_params, _, infections = make_test_data(StableRNG(1), sp.n_t_steps, sp.n_subjects)

obs_df = DataFrame(model_data["observations"])

model_cache = WaningModelCache(obs_df)

obs_df.observed_titre .= 0.0

b_trial = @benchmark waning_curve!(
    $model_params,
    $individual_waning_kucharski!,
    $sp,
    $infections,
    $model_cache,
    $obs_df.observed_titre
)

# Previous: 217 microseconds