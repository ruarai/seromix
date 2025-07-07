

# pt_jld = JLD2.load("runs/hanam_2018/pt_pigeons_5.jld2")
# pt = pt_jld["pt"]


plot(pt.reduced_recorders.index_process[3], linewidth = 2)


chain = Chains(pt);


plot(chain, [:mu_long], seriestype = :density)
vline!([2], lc = "black", linewidth = 2)
plot(chain, [:mu_short], seriestype = :density)

plot(chain, [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain, Symbol("infections[20]"), seriestype = :traceplot)

plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:obs_sd], seriestype = :traceplot)
plot(chain, [:omega], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

plot(chain, [:log_density], seriestype = :traceplot)



heatmap(chain_infections_prob_2(chain, p)')
# # heatmap(model_data["infections_matrix"]')

n_inf = chain_sum_infections(chain, p)
plot(n_inf)

scatter(chain[:mu_long], n_inf)
scatter(chain[:mu_long], chain[:mu_short])
scatter(chain[:mu_long], chain[:tau])
scatter(chain[:tau], n_inf)
scatter(chain[:mu_short], chain[:omega])

using StatsPlots
@df pt.shared.reports.swap_prs StatsPlots.plot(:round, :mean, group = :first, legend = false)


plot(pt.shared.tempering.communication_barriers.localbarrier)



X = chain_infections_prob_2(chain, p)

reindex = [
    vcat(
        zeros(max.(0, -(p.subject_birth_ix[i]))),
        X[max.(p.subject_birth_ix[i], 1):p.n_t_steps, i]
    ) 
    for i in 1:p.n_subjects
]

plot(cumsum.(reindex), legend = false, lc = "black")

plot(cumsum(X, dims = 1), legend = false)