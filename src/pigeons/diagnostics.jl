

chain_name = "pigeons_5_mixis"
pt = JLD2.load("runs/hanam_2018/pt_$chain_name.jld2")["pt"]


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



heatmap(chain_infections_prob_2(chain, sp)')
# # heatmap(model_data["infections_matrix"]')

n_inf = chain_sum_infections(chain, sp)
plot(n_inf)

scatter(chain[:mu_long], n_inf)
scatter(chain[:mu_long], chain[:mu_short])
scatter(chain[:mu_long], chain[:tau])
scatter(chain[:tau], n_inf)
scatter(chain[:mu_short], chain[:omega])

using StatsPlots
@df pt.shared.reports.swap_prs StatsPlots.plot(:round, :mean, group = :first, legend = false)


plot(pt.shared.tempering.communication_barriers.localbarrier)


lpp = model_sum_mixIS(DataFrame(chain), sp, obs_df, model)
