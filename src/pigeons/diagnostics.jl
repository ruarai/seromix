

# chain_name = "pigeons_6_2"
# pt = JLD2.load("runs/hanam_2018/pt_$chain_name.jld2")["pt"]


plot(pt.reduced_recorders.index_process, linewidth = 1, size = (1200, 1200))
plot(pt.reduced_recorders.index_process[4], linewidth = 2)
plot(pt.shared.tempering.communication_barriers.localbarrier)

plot(pt.shared.tempering.communication_barriers.localbarrier)


plot([
    maximum(pt.reduced_recorders.index_process[i]) -
    minimum(pt.reduced_recorders.index_process[i]) 
    for i in 1:length(pt.reduced_recorders.index_process)
])

chain = Chains(pt);
plot(chain, [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain, Symbol("infections[20]"), seriestype = :traceplot)

plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:obs_sd], seriestype = :traceplot)
plot(chain, [:omega], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

plot(chain, [:log_density], seriestype = :traceplot)




# @gif for i in 1:64
#     heatmap(chain_infections_prob_2(chain[:,:,i], sp)', clim = (0, 1))
# end fps = 8

heatmap(chain_infections_prob_2(chain, sp)')
heatmap(rand(prior_infection_dist)')

a = [rand(prior_infection_dist) for i in 1:50]
heatmap([mean(a)[ix_t, ix_subject] for ix_t in 1:sp.n_t_steps, ix_subject in 1:sp.n_subjects]')
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




lpp = model_sum_mixIS(DataFrame(chain), sp, obs_df, model)
