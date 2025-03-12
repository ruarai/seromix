


using AdvancedMH

struct FunctionProposal{issymmetric,P} <: AdvancedMH.Proposal{P}
    proposal::P
end
function FunctionProposal{issymmetric}(proposal) where {issymmetric}
    return FunctionProposal{issymmetric,typeof(proposal)}(proposal)
end

import AdvancedMH: propose

function propose(
    rng::Random.AbstractRNG,
    proposal::FunctionProposal,
    model::AdvancedMH.DensityModelOrLogDensityModel,
    t
)
    return proposal.proposal(rng, t)
end


# import AdvancedMH: q

# # TODO Verify this makes sense ...!
# function q(
#     proposal::FunctionProposal,
#     t,
#     t_cond
# )
#     println("wat")
#     return 1.0
# end

# Symmetric -> ratio of proposal densities is zero
import AdvancedMH: logratio_proposal_density
logratio_proposal_density(::FunctionProposal{true}, state, candidate) = 0
