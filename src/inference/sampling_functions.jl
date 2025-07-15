

function sample_chain(
    model,
    initial_params,
    gibbs_sampler,
    sp,
    rng;
    n_sample::Int,
    n_thinning::Int,
    n_chain::Int,
    progress=true
)
    chain = sample(
        rng,
        model, gibbs_sampler, 
        MCMCThreads(), n_sample ÷ n_thinning, n_chain,

        thinning = n_thinning,
        callback = log_callback,
        initial_params = initial_params,
        progress = progress
    )

    check_inf_prob_birth_year(chain_infections_prob(chain, sp), sp)

    return chain
end

function model_symbols_apart_from(model, syms; as_vector = false)
    symbols = DynamicPPL.syms(DynamicPPL.VarInfo(model))
    for s in syms
        symbols = symbols[findall(symbols .!= s)]
    end

    if as_vector
        return [i for i in symbols]
    else
        return symbols
    end
end


function make_gibbs_sampler(model, sp, proposal_function)
    symbols_not_inf = model_symbols_apart_from(model, [:infections])
    
    gibbs_sampler = Gibbs(
        :infections => make_mh_infection_sampler(sp, proposal_function),
        symbols_not_inf => make_mh_parameter_sampler()
    )
    
    return gibbs_sampler
end

function make_gibbs_sampler_original(model, sp, proposal_function)
    symbols_not_inf = model_symbols_apart_from(model, [:infections])
    
    gibbs_sampler = Gibbs(
        :infections => make_mh_infection_sampler(sp, proposal_function),
        symbols_not_inf => make_mh_parameter_sampler_original()
    )
    
    return gibbs_sampler
end


function make_gibbs_sampler_slice(model, sp, proposal_function)
    symbols_not_inf = model_symbols_apart_from(model, [:infections])

    gibbs_sampler = Gibbs(
        :infections => make_mh_infection_sampler(sp, proposal_function; prop_sample = 1.0, n_repeats = 10),
        symbols_not_inf => externalsampler(RandPermGibbs(SliceSteppingOut(2.0)))
    )
    
    return gibbs_sampler
end


function log_callback(rng, model, sampler, sample, state, iteration; kwargs...)
    if iteration % 50 == 0
        if length(sampler.alg.samplers) > 1

            inf_state = state.states[1].state
            # param_state = state.states[2].state

            pr_accept_inf = inf_state.n_accepted / (inf_state.n_accepted + inf_state.n_rejected)
            pr_accept_param = 0.0

            sigma_covar = 0.0
            # if param_state isa ParameterSamplerState
            #     pr_accept_param = param_state.n_accepted / (param_state.n_accepted + param_state.n_rejected)
            #     sigma_covar = param_state.sigma_covar
            # end

            # Multiply by two to account for gibbs sampling
            sample_time = (inf_state.time_B - inf_state.time_A) * 1000 * 2.0

            # Time per 10,000 samples
            sampling_rate = sample_time / 6

            @printf(
                "%6d; %6.2f; %6.3f; %10.8f; %7.2fms/sample;%7.2fmin/10,000 samples\n", 
                iteration, 
                pr_accept_inf, 
                pr_accept_param, 
                sigma_covar,
                sample_time, sampling_rate
            )
        else
            @printf("%5d;\n", iteration)
        end
    end
end

function get_logp(theta, model; context = DynamicPPL.DefaultContext())
    f = model.logdensity
    varinfo_prev = DynamicPPL.unflatten(f.varinfo, theta)
    return DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))
end


# Hack to make sure logprob is carried through the Gibbs sampler
# function Turing.Inference.varinfo(state::Turing.Inference.TuringState)
#     if state.state.transition isa InfectionSamplerState
#         θ = state.state.transition.θ
#         vi = DynamicPPL.unflatten(state.ldf.varinfo, θ)
#         vi = setlogp!!(vi, state.state.transition.lp)
#     else
#         θ = Turing.Inference.getparams(state.ldf.model, state.state)
#         vi = DynamicPPL.unflatten(state.ldf.varinfo, θ)
#         return vi
#     end
# end

# Hack to reduce model evaluations
function Turing.Inference.transition_to_turing(f::DynamicPPL.LogDensityFunction, transition)
    if transition isa AdvancedMH.Transition
        θ = Turing.Inference.getparams(f.model, transition)
        varinfo = DynamicPPL.unflatten(f.varinfo, θ)
        return Turing.Inference.Transition(f.model, varinfo, transition)
    end

    θ = if transition isa InfectionSamplerTransition
        transition.θ
    elseif transition isa ParameterSamplerTransition
        transition.θ
    elseif transition isa SliceSampling.Transition
        Turing.Inference.getparams(f.model, transition) # TODO this may not work
    end

    varinfo = DynamicPPL.unflatten(f.varinfo, θ)
    return Turing.Inference.Transition(θ, transition.lp, Turing.Inference.getstats(transition))
end

function set_lp!(model, chain)
    chain.value.data[:, findfirst(chain.value.axes[2] .== :lp),:] .= logjoint_threaded(model, chain)
end

function logjoint_threaded(model::Model, chain::AbstractMCMC.AbstractChains)
    var_info = VarInfo(model) # extract variables info from the model
    logp = zeros(size(chain, 1), size(chain, 3))

    @showprogress Threads.@threads for iteration_idx in 1:size(chain, 1)
        for chain_idx in 1:size(chain, 3)
            argvals_dict = OrderedDict(
                vn_parent =>
                    DynamicPPL.values_from_chain(var_info, vn_parent, chain, chain_idx, iteration_idx) for
                vn_parent in keys(var_info)
            )
            logp[iteration_idx, chain_idx] = logjoint(model, argvals_dict)
        end
    end

    return logp
end