using AgentBasedCells
using OrdinaryDiffEq
solver_opts = (stol=1e-4, atol=1e-6, rtol=1e-6, method=Reinsert(), solver=TRBDF2())

using CairoMakie
using Colors
using ColorSchemes
using StatsBase
using LinearAlgebra

include("utils.jl")
include("models.jl")

bopt_wt = skopt.load("data/sim/checkpoint_sos_wt.pkl")
params_wt = pyconvert.(Float64, bopt_wt.x)

bopt_2palsel = skopt.load("data/sim/checkpoint_sos_2palsel.pkl")
params_2palsel = pyconvert.(Float64, bopt_2palsel.x)

sel_funhillfs = AgentBasedCells.gen_division_rate_function(hillfs, bursty_rn)

model_wt = MotherCellModel(bursty_rn, DivisionRateBounded(γτ_wt, 1.0, bursty_rn), BinomialKernel(0.5))
analyticals_wt = run_analytical_single(model_wt, experiment_setup(model_parameters=(params_wt..., 0.0, 0.0, 0.0, 0.0), trn=working_trn); solver_opts...)
division_dist_ancest!(analyticals_wt)
div_dist_wt = division_time_ancest(analyticals_wt)

model_2pal = MotherCellModel(bursty_rn, DivisionRateBounded(γτ_wt*hillfs, 1.0, bursty_rn), BinomialKernel(0.5))
analyticals_2pal_sel = run_analytical_single(model_2pal, experiment_setup(model_parameters=(params_wt..., params_2palsel...), trn=working_trn); solver_opts...)
division_dist_ancest!(analyticals_2pal_sel)
div_dist_2pal_sel = division_time_ancest(analyticals_2pal_sel)

model_naive = MotherCellModel(bursty_rn, DivisionRateBounded(γτ_2pal, 1.0, bursty_rn), BinomialKernel(0.5))
analyticals_naive = run_analytical_single(model_naive, 
    experiment_setup(model_parameters=(params_wt..., 0.0, 0.0, 0.0, 0.0), trn=working_trn); solver_opts...)
division_dist_ancest!(analyticals_naive)
div_dist_naive = division_time_ancest(analyticals_naive)

size_inches = (17, 7)
size_pt = 72 .* size_inches
fig = Figure(resolution=size_pt, fontsize=19)
xs = collect(1:working_trn)
ts = 0.0:0.1:300.0
colors = ColorSchemes.Egypt.colors
transp = 0.4
stairstransp = 0.4

axinter = Axis(fig[1,1]; xlabel="Time", ylabel="Probability density", title="Interdivision time distribution")
axdiv = Axis(fig[1,2]; xlabel="Protein count at division", ylabel="Probability density", title="Division protein distribution")
axjoint = Axis(fig[1,3], xlabel="Cell age at division", ylabel="Protein count at division", title="Division distribution")
axhill = Axis(fig[1,1]; 
    width=Relative(0.30), height=Relative(0.30), halign=0.93, valign=0.9, 
    xlabel="Protein counts x", ylabel="Selection f(x)",
    xlabelsize=16,
    xticksize=10,
    xticklabelsize=16,
    yticklabelsize=16,
    ylabelsize=16,
    yticksize=10)

translate!(axhill.blockscene, 0, 0, 1000)
xlims!(axhill, (0.0, 150))
ylims!(axhill, low=0.0)

xlims!(axdiv, (0, 150))
ylims!(axdiv, low=0.0)
div_hist_wt = normalize(fit(Histogram, df_div_wt[:, :Column3], 1:200); mode=:pdf)
stairs!(axdiv, collect(midpoints(div_hist_wt.edges[1])), div_hist_wt.weights; color=(colors[2], stairstransp), step=:center)
barplot!(axdiv, collect(midpoints(div_hist_wt.edges[1])), div_hist_wt.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(axdiv, xs, analyticals_wt.results[:division_dist]; color=colors[2], linewidth=3.0)

div_hist_2pal = normalize(fit(Histogram, df_div_2pal[:, :Column3], 1:200); mode=:pdf)
stairs!(axdiv, collect(midpoints(div_hist_2pal.edges[1])), div_hist_2pal.weights; color=(colors[1], stairstransp), step=:center)
barplot!(axdiv, collect(midpoints(div_hist_2pal.edges[1])), div_hist_2pal.weights; 
    color=(colors[1], transp), strokecolor=(colors[1], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(axdiv, xs, analyticals_2pal_sel.results[:division_dist]; color=colors[1], linewidth=3.0)
#lines!(axdiv, xs, analyticals_naive.results[:division_dist]; color=colors[2], linewidth=3.0, linestyle=:dot)

interdiv_hist_wt = normalize(fit(Histogram, interdiv_times_wt .* 8, 1:8:maximum(interdiv_times_wt .* 8)); mode=:pdf)
stairs!(axinter, collect(midpoints(interdiv_hist_wt.edges[1])), interdiv_hist_wt.weights; color=(colors[2], stairstransp), step=:center)
barplot!(axinter, collect(midpoints(interdiv_hist_wt.edges[1])), interdiv_hist_wt.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(axinter, ts, div_dist_wt.(ts); color=colors[2], linewidth=3.0)

interdiv_hist_2pal = normalize(fit(Histogram, interdiv_times_2pal .* 8, 1:8:maximum(interdiv_times_2pal .* 8)); mode=:pdf)
stairs!(axinter, collect(midpoints(interdiv_hist_2pal.edges[1])), interdiv_hist_2pal.weights; color=(colors[1], stairstransp), step=:center)
barplot!(axinter, collect(midpoints(interdiv_hist_2pal.edges[1])), interdiv_hist_2pal.weights; 
    color=(colors[1], transp), strokecolor=(colors[1], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(axinter, ts, div_dist_2pal_sel.(ts); color=colors[1], linewidth=3.0)
#lines!(axinter, ts, div_dist_naive.(ts); color=colors[2], linewidth=3.0, linestyle=:dot)
#lines!(axinter, ts, div_dist_2pal_adapt.(ts); color=colors[3], linewidth=3.0, linestyle=:dot)
xlims!(axinter, (0, 200))
ylims!(axinter, low=0.0)

lines!(axhill, xs, sel_funhillfs.(xs, Ref([params_wt..., params_2palsel...]), Ref(0.0)), linewidth=3.0, color=colors[1])
#lines!(axhill, xs, fx_.(xs; ps=ps_2pal_2palparams[2:end], n=2); color=colors[3], linestyle=:dot, linewidth=3.0)

hitting_wt = hcat([[length(traj) * tslice - 0.5 * tslice, traj[end] ] 
                   for traj in cycles_wt]...)
hitting_2pal = hcat([[length(traj) * tslice - 0.5 * tslice, traj[end] ] 
                     for traj in cycles_2pal]...)

#hist_2pal = fit(Histogram, Tuple(eachrow(hitting_2pal)), nbins=100)

hist_2pal = kde(Tuple(eachrow(hitting_2pal)); boundary=((0,400),(0, 200)), bandwidth=(11, 3))
hist_wt = kde(Tuple(eachrow(hitting_wt)))

cmap2 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 0.8), length=15)
contourf!(axjoint, hist_2pal; colormap=cmap2, levels=4)
#contourf!(axjoint, hist_wt; colormap=cmap2, levels=5)

cmap2 = range(alphacolor(colors[1], 0.3), stop=ARGB(colors[1], 1.0), length=15)
zs_2pal_sel = hcat(fpt_dist_ancest(analyticals_2pal_sel).(0.0:1.0:200.0)...)
zs_wt = hcat(fpt_dist_ancest(analyticals_wt).(0.0:1.0:200.0)...)
contour!(axjoint, 0.0:1.0:200.0, 1:working_trn, zs_2pal_sel'; colormap=cmap2, linewidth=3, levels=4)
#contourf!(axjoint, 0.0:1.0:200.0, 1:working_trn, zs_2pal_sel'; colormap=cmap2, linewidth=3, levels=0.3:0.1:1, mode = :relative)
#contour!(axjoint, 0.0:1.0:200.0, 1:200, zs_wt'; colormap=cmap2, linewidth=2, levels=4)
#
#cmap3 = range(alphacolor(colors[3], 0.0), stop=ARGB(colors[3], 1.0), length=15)
#zs_naive = hcat(fpt_dist(analyticals_naive).(0.0:1.0:200.0)...)
#contour!(axjoint, 0.0:1.0:200.0, 1:200, zs_naive'; colormap=cmap3, linewidth=2, levels=4, linestyle=:dot)

idxs_2pal = rand(1:length(cycles_2pal), 30)
#for idx in idxs_2pal[1:5]
#    lines!(axjoint, (collect(1:length(cycles_2pal[idx])) .- 0.5) .* tslice , cycles_2pal[idx]; color=(:black, 0.8), linewidth=2)
#    scatter!(axjoint, (length(cycles_2pal[idx]) -0.5) * tslice, cycles_2pal[idx][end]; color=(:black, 0.8), marker=:cross)
#end
for idx in idxs_2pal
    lines!(axjoint, (collect(1:length(cycles_2pal[idx])) .- 0.5) .* tslice , cycles_2pal[idx]; color=(:black, 0.2), linewidth=2)
    scatter!(axjoint, (length(cycles_2pal[idx]) -0.5) * tslice, cycles_2pal[idx][end]; color=(:black, 0.2), marker=:cross)
end


xlims!(axjoint, (25.0, 200.0))
ylims!(axjoint, (0.0, 100.0))

data_legend = [PolyElement(color=colors[2]), PolyElement(color=(colors[1]))]
fitting_legend = [LineElement(color=colors[2]), LineElement(color=colors[1])]

Legend(fig[2,:], [data_legend, fitting_legend], [["Wild type", "DNA damage"], ["Wild type", "DNA damage", ]], ["Data", "Model fitting"]; 
    orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true, groupgap=60)

#Legend(fig[2, 1:3], axbirth; orientation = :vertical, tellwidth = false, tellheight = true, framevisible=false)

CairoMakie.activate!()
mkpath("$(plotsdir())/growth-sos")
save("$(plotsdir())/growth-sos/fitting.pdf", fig, pt_per_unit = 1)
save("$(plotsdir())/growth-sos/fitting.svg", fig, pt_per_unit = 1)

size_inches = (17, 6)
size_pt = 72 .* size_inches
fig_fpts_naive = Figure(resolution=size_pt, fontsize=14)
xs = collect(1:working_trn)
ts = 0.0:0.1:300.0
colors = ColorSchemes.Egypt.colors
transp = 0.4
stairstransp = 0.4

ax_inter_naive = Axis(fig_fpts_naive[1,1]; xlabel="Time", ylabel="Probability density", title="Interdivision time distribution")
ax_div_naive = Axis(fig_fpts_naive[1,2]; xlabel="Protein counts at division", ylabel="Probability density", title="Division protein distribution")
ax_wt_naive = Axis(fig_fpts_naive[1,3]; xlabel="Cell age at division", ylabel="Protein count at division", title="Division distribution")
ax_fpt_naive = Axis(fig_fpts_naive[1,4]; xlabel="Cell age at division", ylabel="Protein count at division", title="Division distribution")

stairs!(ax_div_naive, collect(midpoints(div_hist_wt.edges[1])), div_hist_wt.weights; color=(colors[1], stairstransp), step=:center)
barplot!(ax_div_naive, collect(midpoints(div_hist_wt.edges[1])), div_hist_wt.weights; 
    color=(colors[1], transp), strokecolor=(colors[1], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_div_naive, xs, analyticals_wt.results[:division_dist]; color=colors[1], linewidth=3.0)

stairs!(ax_div_naive, collect(midpoints(div_hist_2pal.edges[1])), div_hist_2pal.weights; color=(colors[2], stairstransp), step=:center)
barplot!(ax_div_naive, collect(midpoints(div_hist_2pal.edges[1])), div_hist_2pal.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_div_naive, xs, analyticals_naive.results[:division_dist]; color=colors[2], linewidth=3.0)

xlims!(ax_div_naive, (0, 150))
ylims!(ax_div_naive, low=0.0)

stairs!(ax_inter_naive, collect(midpoints(interdiv_hist_wt.edges[1])), interdiv_hist_wt.weights; color=(colors[1], stairstransp), step=:center)
barplot!(ax_inter_naive, collect(midpoints(interdiv_hist_wt.edges[1])), interdiv_hist_wt.weights; 
    color=(colors[1], transp), strokecolor=(colors[1], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_inter_naive, ts, div_dist_wt.(ts); color=colors[1], linewidth=3.0)

stairs!(ax_inter_naive, collect(midpoints(interdiv_hist_2pal.edges[1])), interdiv_hist_2pal.weights; color=(colors[2], stairstransp), step=:center)
barplot!(ax_inter_naive, collect(midpoints(interdiv_hist_2pal.edges[1])), interdiv_hist_2pal.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_inter_naive, ts, div_dist_naive.(ts); color=colors[2], linewidth=3.0)

xlims!(ax_inter_naive, (0, 200))
ylims!(ax_inter_naive, low=0.0)

cmap1 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 0.8), length=15)
cmap2 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 0.8), length=15)

contourf!(ax_wt_naive, hist_wt; colormap=cmap1, levels=5)
contour!(ax_wt_naive, 0.0:1.0:200.0, 1:working_trn, zs_wt'; colormap=cmap1, linewidth=3, levels=4)

contourf!(ax_fpt_naive, hist_2pal; colormap=cmap2, levels=5)
zs_naive = hcat(fpt_dist_ancest(analyticals_naive).(0.0:1.0:200.0)...)
contour!(ax_fpt_naive, 0.0:1.0:200.0, 1:working_trn, zs_naive'; colormap=cmap2, linewidth=3, levels=4)

xlims!(ax_wt_naive, (25.0, 200.0))
ylims!(ax_wt_naive, (0.0, 100.0))

xlims!(ax_fpt_naive, (25.0, 200.0))
ylims!(ax_fpt_naive, (0.0, 100.0))

Legend(fig_fpts_naive[2,:], [data_legend, fitting_legend], [["Wild type", "DNA damage"], ["Wild type", "No selection", ]], ["Data", "Model fitting"]; 
    orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true, groupgap=60)

save("$(plotsdir())/growth-sos/apx_fpts_naive.pdf", fig_fpts_naive, pt_per_unit = 1)
