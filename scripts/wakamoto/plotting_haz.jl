using AgentBasedFSP
using OrdinaryDiffEq
solver_opts = (stol=1e-4, atol=1e-6, rtol=1e-6, method=Reinsert(), solver=TRBDF2())
using CairoMakie
using ColorSchemes
using Colors
using StatsBase
using LinearAlgebra

include("models.jl")
include("utils.jl")

params₊Sm = [0.7978187613870786, 1.326304795656268]
model_spline₊Sm = CellPopulationModel(bursty_rn_base, DivisionRateBounded(γτ₋Sm*itp₊Sm_(Protein),1.0,bursty_rn_base), BinomialKernel(0.5))
#analyticals₊Sm_spl = run_analytical_single(model_spline₊Sm, experiment_setup(model_parameters=params₊Sm, trn=working_trn); solver_opts...)
division_dist_ancest!(analyticals₊Sm_spl)
div_dist₊Sm_spl = division_time_ancest(analyticals₊Sm_spl)

params_w₋Sm = [ps_bursty₋Sm..., [1.5758894318598402, 92.08532187610685, 0.4704446128065254, 2]...]
model_bursty_hillfs = CellPopulationModel(bursty_rn, DivisionRateBounded(γτ₋Sm*hillfs,1.0,bursty_rn), BinomialKernel(0.5))
#analyticals₊Sm_sel = run_analytical_single(model_bursty_hillfs, experiment_setup(model_parameters=params_w₋Sm, trn=working_trn); solver_opts...)
division_dist_ancest!(analyticals₊Sm_sel)
div_dist₊Sm_sel = division_time_ancest(analyticals₊Sm_sel)

params_w₊Smhillfs = [0.6938520229435957, 2.1203328315585255, 1.9464579944214924, 101.62383803442542, 0.014570213174131217, 5]
analyticals₊Sm_adaptfs = run_analytical_single(model_bursty_hillfs, experiment_setup(model_parameters=params_w₊Smhillfs, trn=working_trn); solver_opts...)
division_dist_ancest!(analyticals₊Sm_adaptfs)
div_dist₊Sm_adaptfs = division_time_ancest(analyticals₊Sm_adaptfs)
joint_fpt_ancest_cdist!(analyticals₊Sm_adaptfs; step=tslice)

size_inches = (17, 7)
size_pt = 72 .* size_inches
fig = Figure(resolution=size_pt, fontsize=19)
xs = collect(1:working_trn)
ts = 0.0:0.1:150.0
colors = ColorSchemes.Egypt.colors
transp = 0.4
stairstransp = 0.4

ax_interdiv = Axis(fig[1,1], xlabel="Division time", ylabel="Probability density", title="Interdivision time distribution")
ax_div = Axis(fig[1,2], xlabel="Protein counts", ylabel="Probability density", title="Division distribution")
ax_joint = Axis(fig[1,3], xlabel="Cell age", ylabel="Protein count", title="First passage time distribution")
ax_sel = Axis(fig[1,1]; 
    width=Relative(0.3), height=Relative(0.3), halign=0.93, valign=0.92, 
    xlabel="Protein counts x", ylabel="Selection f(x)", xticks=0:100:200,
    xlabelsize=16,
    xticksize=10,
    xticklabelsize=16,
    yticklabelsize=16,
    ylabelsize=16,
    yticksize=10)
translate!(ax_sel.blockscene, 0, 0, 1000)
xlims!(ax_sel, (0.0, 200))
ylims!(ax_sel, low=0.0)

xs = collect(1:working_trn)

div_hist₊Sm = normalize(fit(Histogram, df_div₊Sm[:, :Column3], 1:4:200); mode=:pdf)
stairs!(ax_div, collect(midpoints(div_hist₊Sm.edges[1])), div_hist₊Sm.weights; color=(colors[2], stairstransp), step=:center)
barplot!(ax_div, collect(midpoints(div_hist₊Sm.edges[1])), div_hist₊Sm.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)

#lines!(ax_div, xs, analyticals₊Sm_spl.results[:division_dist_ancest]; color=colors[2], linewidth=3.0, linestyle=:dash)
#lines!(ax_div, xs, analyticals₊Sm_sel.results[:division_dist_ancest]; color=colors[2], linewidth=3.0)
lines!(ax_div, xs, analyticals₊Sm_adaptfs.results[:division_dist_ancest]; color=colors[1], linewidth=3.0)
#lines!(ax_div, xs, analyticals₊Sm_adaptfs.results[:division_dist_ancest]; color=colors[2], linewidth=3.0)
#lines!(ax_div, xs, analyticals₊Sm_adaptfsrs.results[:division_dist_ancest]; color=colors[2], linewidth=3.0)
#lines!(ax_div, xs, analyticals₊Sm_naive.results[:division_dist]; color=colors[2], linewidth=3.0)
xlims!(ax_div, (50, 200))
ylims!(ax_div, low=0.0)

ts = 0.0:0.1:150.0

interdiv_hist₊Sm = normalize(fit(Histogram, interdiv_times₊Sm .* tslice .- 0.5 * tslice, 1:5:maximum(interdiv_times₊Sm .* 5)); mode=:pdf)
stairs!(ax_interdiv, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; color=(colors[2], stairstransp), step=:center)
barplot!(ax_interdiv, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)

#inter_cdist = sum.(analyticals₊Sm_adaptfs.results[:joint_fpt_ancest_cdist])

#lines!(ax_interdiv, range(0.0, length=length(inter_cdist), step=tslice), inter_cdist)
lines!(ax_interdiv, ts, div_dist₊Sm_spl.(ts); color=colors[2], linewidth=3.0, linestyle=:dash)
#lines!(ax_interdiv, ts, div_dist₊Sm_sel.(ts); color=colors[2], linewidth=3.0)
lines!(ax_interdiv, ts, div_dist₊Sm_adaptfs.(ts); color=colors[1], linewidth=3.0)
#lines!(ax_interdiv, ts, div_dist₊Sm_adaptfsrs.(ts); color=colors[2], linewidth=3.0)
xlims!(ax_interdiv, (25, 150))
ylims!(ax_interdiv, low=0.0)

#lines!(ax_sel, xs, sel_funhillfs.(xs, Ref(params_w₋Sm), Ref(0.0)), linewidth=3.0, color=colors[2], linestyle=:dash)
#lines!(ax_sel, xs, sel_funhillfsrs.(xs, Ref(params_w₊Smhillfs), Ref(0.0)), linewidth=3.0, color=colors[2])

cycles₋Sm = vcat(make_cell_cycle.(lineages₋Sm)...); 
cycles₊Sm = vcat(make_cell_cycle.(lineages₊Sm)...); 

fluors₋Sm = fluorescence_trajectory.(cycles₋Sm)
fluors₊Sm = fluorescence_trajectory.(cycles₊Sm)

hitting₋Sm = hcat([[length(traj) * tslice - 0.5 * tslice, traj[end] ./ cf₋Sm] for traj in fluors₋Sm]...)
hitting₊Sm = hcat([[length(traj) * tslice - 0.5 * tslice, traj[end] ./ cf₋Sm] for traj in fluors₊Sm]...)

hist₋Sm = kde(Tuple(eachrow(hitting₋Sm))) #fit(Histogram, Tuple(eachrow(hitting₋Sm)), nbins=50)
hist₊Sm = kde(Tuple(eachrow(hitting₊Sm)))

cmap2 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 0.8), length=15)
contourf!(ax_joint, hist₊Sm; colormap=cmap2, levels=5)

#zs₊Sm_adaptfsrs = hcat(fpt_dist_ancest(analyticals₊Sm_adaptfsrs).(0.0:1.0:125.0)...)
#cmap3 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 1.0), length=15)
#contour!(ax_joint, 0.0:1.0:125.0, 1:length(analyticals₊Sm_adaptfsrs.results[:birth_dist]), zs₊Sm_adaptfsrs'; colormap=cmap3, linewidth=3, levels=4)

zs₊Sm_adaptfs = hcat(fpt_dist_ancest(analyticals₊Sm_adaptfs).(0.0:1.0:125.0)...)
cmap3 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 1.0), length=15)
contour!(ax_joint, 0.0:1.0:125.0, 1:length(analyticals₊Sm_adaptfs.results[:birth_dist]), zs₊Sm_adaptfs'; colormap=cmap3, linewidth=3, levels=4)

#zs₊Sm_spl = hcat(fpt_dist_ancest(analyticals₊Sm_spl).(0.0:1.0:125.0)...)
#cmap1 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 1.0), length=15)
#contour!(ax_joint, 0.0:1.0:125.0, 1:length(analyticals₊Sm_spl.results[:birth_dist]), zs₊Sm_spl'; colormap=cmap1, linewidth=3, levels=4)

xlims!(ax_joint, (25.0, 125.0))
ylims!(ax_joint, (50.0, 150.0))

lines!(ax_sel, xs, itp₊Sm.(xs); colors=colors[2])
lines!(ax_sel, xs, sel_funhillfs.(xs, Ref(params_w₊Smhillfs), Ref(0.0)); colors=colors[2], linestyle=:dash)

mkpath("$(plotsdir())/Wakamoto")
CairoMakie.activate!()
save("$(plotsdir())/Wakamoto/fitting_haz.pdf", fig, pt_per_unit = 1)

