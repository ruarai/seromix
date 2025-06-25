

include("../dependencies.jl")
include("util.jl")

using Test, StableRNGs
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


complete_obs.observed_titre

@testset "Waning model" begin
    n_t_steps = 40
    n_subjects = 40
    model_params, p, infections = make_test_data(StableRNG(1), n_t_steps, n_subjects)

    complete_obs = expand_grid(
        ix_t_obs = 1:n_t_steps, ix_strain = 1:n_t_steps, ix_subject = 1:n_subjects,
        observed_titre = 0.0
    )

    waning_curve!(
        model_params.mu_long, model_params.mu_short, model_params.omega,
        model_params.sigma_long, model_params.sigma_short, model_params.tau,

        p.antigenic_distances, p.time_diff_matrix, p.subject_birth_ix,

        infections,

        make_obs_lookup(complete_obs), make_obs_views(complete_obs),
        complete_obs.observed_titre
    )

    # Basic sanity checks
    @test all(complete_obs.observed_titre .>= 0)

    # Lazy way to store the approx. expected result
    @test sum(complete_obs.observed_titre) ≈ 52164.14374999999

    b_trial = @benchmark waning_curve!(
        $model_params.mu_long, $model_params.mu_short, $model_params.omega,
        $model_params.sigma_long, $model_params.sigma_short, $model_params.tau,

        $p.antigenic_distances, $p.time_diff_matrix, $p.subject_birth_ix,

        $infections,

        $(make_obs_lookup(complete_obs)), $(make_obs_views(complete_obs)),
        $complete_obs.observed_titre
    )

    # Should take less than 5ms, and memory use should be minimal
    @test (median(x).time / 1000) < 5_000
    @test (median(x).memory) < 100
    @test (median(x).allocs) < 10
end




