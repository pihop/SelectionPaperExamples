using AgentBasedCells
using OrdinaryDiffEq
solver_opts = (stol=1e-4, atol=1e-6, rtol=1e-6, method=Reinsert(), solver=TRBDF2())
using CairoMakie
using ColorSchemes
using Colors
using StatsBase
using LinearAlgebra

include("models.jl")
include("utils.jl")

tslice = 5
bopt_bursty₋Sm = skopt.load("data/sim/checkpoint_$(strain)_notreatment.pkl")
params_bursty₋Sm = pyconvert.(Float64, bopt_bursty₋Sm.x)

bopt_bursty₊Sm_mult_hill = skopt.load("data/sim/checkpoint_$(strain)_hillsel.pkl")
params_bursty₊Sm_mult_hill = pyconvert.(Float64, bopt_bursty₊Sm_mult_hill.x)

bopt_bursty₊Sm_mult_adapt = skopt.load("data/sim/checkpoint_$(strain)_hilladapt.pkl")
params_bursty₊Sm_adapt = pyconvert.(Float64, bopt_bursty₊Sm_mult_adapt.x)
params_bursty₊Sm_adapt = [ 0.6999621927855469, 2.072504215985731, 1.464313090526479, 92.79936824735046, 0.0, 4.0]

model_bursty₋Sm = CellPopulationModel(
    bursty_rn, DivisionRateBounded(γτ₋Sm,1.0,bursty_rn), BinomialKernel(0.5))
model_bursty_hillfs = CellPopulationModel(
    bursty_rn, DivisionRateBounded(γτ₋Sm*hillfs,1.0,bursty_rn), BinomialKernel(0.5))
model_bursty_naive = CellPopulationModel(
    bursty_rn, DivisionRateBounded(γτ₊Sm,1.0,bursty_rn), BinomialKernel(0.5))
#model_bursty_hillfsrs = CellPopulationModel(bursty_rn, DivisionRateBounded(γτ₋Sm*hillfs*hillfr,1.0,bursty_rn), BinomialKernel(0.5))

#sel_funhillfsrs = AgentBasedCells.gen_division_rate_function(hillfs*hillfr, bursty_rn)
#sel_funhillfs = AgentBasedCells.gen_division_rate_function(hillfs, bursty_rn)

# With -Sm parameters (selection)
#params_w₋Sm = [ps_bursty₋Sm..., 1.1837173382640955, 78.04605744215044, 0.6082567620480603, 20, 0, 0, 0, 0]
#params_w₋Sm = [ps_bursty₋Sm..., [1.5758894318598402, 92.08532187610685, 0.4704446128065254, 2]..., 0, 0, 0, 0]
#params_w₊Smhillfs = [0.7324267886429983, 1.5935188957765232, 0.7352902257893541, 36.08820301303741, 0.49711252456840127, 3, 0, 0, 0, 0]
#params_w₊Smhillfs = [0.7016368240570192, 2.0666994047047123, 1.3617244511243656, 87.55734637539929, 0.0, 5, 0, 0, 0, 0]

#params_w₊Smhillfs = [0.6871638588495048, 2.058707195201312, 1.8888388575562163, 102.77794532513391, 0.07959460183465, 5, 0,0,0,0]
# With +Sm parameters (adaptation)
#params_w₊Smhillfsrs = [0.7178546086984883, 2.256696475214855, 2.341219937202911, 120.93142980273701, 0.08257371504139617, 4, 0.14323538395606195, 173.09755995183392, 7]
#params_w₊Smhillfsrs = [0.7108217445435515, 2.576964399305984, 1.0880363309363426, 94.06189046498318, 0.34065787736155206, 12, 0.665034709628048, 207.15352119634306, 4]
#params_w₊Smhillfsrs = [0.6995410934354995, 1.9195087814493856, 2.1193521256770906, 109.78267964419724, 0.0, 4, 1.0, 250.0, 50]

params_w₋Sm = [params_bursty₋Sm..., params_bursty₊Sm_mult_hill...]
params_w₊Smhillfs = params_bursty₊Sm_adapt

#NO TREATMENT
analyticals₋Sm = run_analytical_single(
    model_bursty₋Sm, 
    experiment_setup(
        model_parameters=(params_bursty₋Sm..., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 
        trn=(working_trn,));
    solver_opts...)
division_dist_ancest!(analyticals₋Sm)
div_dist₋Sm = division_time_ancest(analyticals₋Sm)

# TREATMENT SELECTION
analyticals₊Sm_sel = run_analytical_single(
    model_bursty_hillfs, 
    experiment_setup(model_parameters=(params_w₋Sm...,), trn=(working_trn,)); 
    solver_opts...)
division_dist_ancest!(analyticals₊Sm_sel)
div_dist₊Sm_sel = division_time_ancest(analyticals₊Sm_sel)

# TREATMENT ADAPTATION
analyticals₊Sm_adaptfs = run_analytical_single(
    model_bursty_hillfs, 
    experiment_setup(model_parameters=(params_w₊Smhillfs...,), trn=(working_trn, )); 
    solver_opts...)
division_dist_ancest!(analyticals₊Sm_adaptfs)
div_dist₊Sm_adaptfs = division_time_ancest(analyticals₊Sm_adaptfs)

## TREATMENT ADAPTATION
#analyticals₊Sm_adaptfsrs = run_analytical_single(model_bursty_hillfsrs, experiment_setup(model_parameters=(params_w₊Smhillfsrs...,), trn=(working_trn,)); solver_opts...)
#division_dist_ancest!(analyticals₊Sm_adaptfsrs)
#div_dist₊Sm_adaptfsrs = division_time_ancest(analyticals₊Sm_adaptfsrs)

## TREATMENT NAIVE
#analyticals₊Sm_naive = run_analytical_single(model_bursty_naive, experiment_setup(model_parameters=(params_w₋Sm..., 0.0, 0.0, 0.0, 0.0), trn=(working_trn,)); solver_opts...)
##analyticals₊Sm_naive = run_analytical_single(model_bursty_naive, experiment_setup(model_parameters=[params₊Sm..., 0, 0, 0, 0 0, 0, 0, 0], trn=working_trn); solver_opts...)
#division_dist_ancest!(analyticals₊Sm_naive)
#div_dist₊Sm_naive = division_time_ancest(analyticals₊Sm_naive)

size_inches = (17, 7)
size_pt = 72 .* size_inches
fig = Figure(resolution=size_pt, fontsize=19)
xs = collect(1:working_trn)
ts = 0.0:0.1:150.0
colors = ColorSchemes.Egypt.colors
transp = 0.4
stairstransp = 0.4

ax_interdiv = Axis(fig[1,1], xlabel="Time", ylabel="Probability density", title="Interdivision time distribution")
#ax_birth = Axis(fig[1,2], xlabel="Protein counts", ylabel="Probability density", title="Birth distribution")
ax_div = Axis(fig[1,2], xlabel="Protein count at division", ylabel="Probability density", title="Division protein distribution")
ax_joint = Axis(fig[1,3], xlabel="Cell age at division", ylabel="Protein count at division", title="Division distribution")
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
div_hist₋Sm = normalize(fit(Histogram, df_div₋Sm[:, :Column3], 1:4:200); mode=:pdf)
stairs!(ax_div, collect(midpoints(div_hist₋Sm.edges[1])), div_hist₋Sm.weights; color=(colors[2], stairstransp), step=:center)
barplot!(ax_div, collect(midpoints(div_hist₋Sm.edges[1])), div_hist₋Sm.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_div, xs, analyticals₋Sm.results[:division_dist_ancest]; color=colors[2], linewidth=3.0)

div_hist₊Sm = normalize(fit(Histogram, df_div₊Sm[:, :Column3], 1:4:200); mode=:pdf)
stairs!(ax_div, collect(midpoints(div_hist₊Sm.edges[1])), div_hist₊Sm.weights; color=(colors[1], stairstransp), step=:center)
barplot!(ax_div, collect(midpoints(div_hist₊Sm.edges[1])), div_hist₊Sm.weights; 
    color=(colors[1], transp), strokecolor=(colors[1], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_div, xs, analyticals₊Sm_sel.results[:division_dist_ancest]; color=colors[1], linewidth=3.0, linestyle=:dash)
lines!(ax_div, xs, analyticals₊Sm_adaptfs.results[:division_dist_ancest]; color=colors[1], linewidth=3.0)
#lines!(ax_div, xs, analyticals₊Sm_adaptfsrs.results[:division_dist_ancest]; color=colors[2], linewidth=3.0)
#lines!(ax_div, xs, analyticals₊Sm_naive.results[:division_dist]; color=colors[2], linewidth=3.0)
xlims!(ax_div, (50, 200))
ylims!(ax_div, low=0.0)

ts = 0.0:0.1:150.0
interdiv_hist₋Sm = normalize(fit(Histogram, interdiv_times₋Sm .* 5, 1:5:maximum(interdiv_times₋Sm .* 5)); mode=:pdf)
stairs!(ax_interdiv, collect(midpoints(interdiv_hist₋Sm.edges[1])), interdiv_hist₋Sm.weights; 
    color=(colors[2], stairstransp), step=:center)
barplot!(ax_interdiv, collect(midpoints(interdiv_hist₋Sm.edges[1])), interdiv_hist₋Sm.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_interdiv, ts, div_dist₋Sm.(ts); 
    color=colors[2], linewidth=3.0)

interdiv_hist₊Sm = normalize(fit(Histogram, interdiv_times₊Sm .* 5, 1:5:maximum(interdiv_times₊Sm .* 5)); mode=:pdf)
stairs!(ax_interdiv, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; 
    color=(colors[1], stairstransp), step=:center)
barplot!(ax_interdiv, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; 
    color=(colors[1], transp), strokecolor=(colors[1], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_interdiv, ts, div_dist₊Sm_sel.(ts); 
    color=colors[1], linewidth=3.0, linestyle=:dash)
lines!(ax_interdiv, ts, div_dist₊Sm_adaptfs.(ts); 
    color=colors[1], linewidth=3.0)
xlims!(ax_interdiv, (25, 150))
ylims!(ax_interdiv, low=0.0)

#lines!(ax_sel, xs, fx_.(xs; ps=params_w₋Sm[3:end]), linewidth=3.0, linestyle=:dash, color=colors[2])
#lines!(ax_sel, xs, fx_.(xs; ps=params_w₊Sm[3:end]), linewidth=3.0, color=colors[2])
lines!(ax_sel, xs, sel_funhillfs.(xs, Ref(params_w₋Sm), Ref(0.0)), linewidth=3.0, color=colors[1], linestyle=:dash)
lines!(ax_sel, xs, sel_funhillfs.(xs, Ref(params_w₊Smhillfs), Ref(0.0)), linewidth=3.0, color=colors[1])
xlims!(ax_sel, (0, 200))

cycles₋Sm = rmnothing(vcat(make_cell_cycle.(lineages₋Sm)...));
cycles₊Sm = rmnothing(vcat(make_cell_cycle.(lineages₊Sm)...));

fluors₋Sm = fluorescence_trajectory.(cycles₋Sm)
fluors₊Sm = fluorescence_trajectory.(cycles₊Sm)

hitting₋Sm = hcat([[length(traj) * tslice .- 0.5 * tslice, traj[end] ./ cf₋Sm] for traj in fluors₋Sm]...)
hitting₊Sm = hcat([[length(traj) * tslice .- 0.5 * tslice, traj[end] ./ cf₊Sm] for traj in fluors₊Sm]...)

hist₋Sm = kde(Tuple(eachrow(hitting₋Sm))) #fit(Histogram, Tuple(eachrow(hitting₋Sm)), nbins=50)
hist₊Sm = kde(Tuple(eachrow(hitting₊Sm)))

cmap2 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 0.8), length=15)
contourf!(ax_joint, hist₊Sm; colormap=cmap2, levels=5)

zs₊Sm_adaptfs = hcat(fpt_dist_ancest(analyticals₊Sm_adaptfs).(0.0:1.0:125.0)...)
cmap3 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 1.0), length=15)
contour!(ax_joint, 0.0:1.0:125.0, 1:length(analyticals₊Sm_adaptfs.results[:birth_dist]), zs₊Sm_adaptfs'; colormap=cmap3, linewidth=3, levels=4)

#zs₊Sm_adaptfsrs = hcat(fpt_dist_ancest(analyticals₊Sm_adaptfs).(0.0:1.0:125.0)...)
#cmap3 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 1.0), length=15)
#contour!(ax_joint, 0.0:1.0:125.0, 1:length(analyticals₊Sm_adaptfsrs.results[:birth_dist]), zs₊Sm_adaptfsrs'; colormap=cmap3, linewidth=3, levels=4)

idxs₊Sm = rand(1:length(fluors₊Sm), 30)
#for idx in idxs₊Sm[1:5] 
#    lines!(ax_joint, (collect(1:length(fluors₊Sm[idx])) .- 0.5) .* tslice , fluors₊Sm[idx] ./ cf₋Sm; color=(:black, 0.8), linewidth=2)
#    scatter!(ax_joint, (length(fluors₊Sm[idx]) -0.5) * tslice, fluors₊Sm[idx][end] ./ cf₋Sm; color=(:black, 0.8), marker=:cross)
#end

for idx in idxs₊Sm[1:end] 
    lines!(ax_joint, (collect(1:length(fluors₊Sm[idx])) .- 0.5) .* tslice , fluors₊Sm[idx] ./ cf₋Sm; color=(:black, 0.2), linewidth=2)
    scatter!(ax_joint, (length(fluors₊Sm[idx]) -0.5) * tslice, fluors₊Sm[idx][end] ./ cf₋Sm; color=(:black, 0.2), marker=:cross)
end


#zs₊Sm_sel = hcat(fpt_dist(analyticals₊Sm_sel).(0.0:1.0:125.0)...)
#cmap3 = range(alphacolor(colors[3], 0.0), stop=ARGB(colors[3], 1.0), length=15)
#contour!(ax_joint, 0.0:1.0:125.0, 1:length(analyticals₊Sm_sel.results[:birth_dist]), zs₊Sm_sel'; colormap=cmap3, linewidth=2, levels=4)

#zs₊Sm_naive = hcat(fpt_dist(analyticals₊Sm_naive).(0.0:1.0:125.0)...)
#cmap3 = range(alphacolor(colors[3], 0.0), stop=ARGB(colors[3], 1.0), length=15)
#contour!(ax_joint, 0.0:1.0:125.0, 1:length(analyticals₊Sm_naive.results[:birth_dist]), zs₊Sm_naive'; colormap=cmap2, linewidth=2, levels=4)

xlims!(ax_joint, (25.0, 125.0))
ylims!(ax_joint, (50.0, 150.0))

data_legend = [PolyElement(color=colors[2]), PolyElement(color=(colors[1]))]
fitting_legend = [LineElement(color=colors[2]), LineElement(color=colors[1]), LineElement(color=colors[1], linestyle=:dash)]

Legend(fig[2,:], [data_legend, fitting_legend], [["No treatment", "Treatment"], ["No treatment", "Selection + adaptation", "Selection"]], ["Data", "Model fitting"]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true, groupgap=60)

#Legend(fig[2, 1:3], ax_birth, orientation = :vertical, tellwidth = false, tellheight = true, framevisible=false)
mkpath("$(plotsdir())/Wakamoto")
CairoMakie.activate!()
save("$(plotsdir())/Wakamoto/fitting.pdf", fig, pt_per_unit = 1)
save("$(plotsdir())/Wakamoto/fitting.svg", fig, pt_per_unit = 1)

size_inches = (17, 6)
size_pt = 72 .* size_inches
fig_naive = Figure(resolution=size_pt, fontsize=14)
xs = collect(1:working_trn)
ts = 0.0:0.1:150.0
colors = ColorSchemes.Egypt.colors
transp = 0.4
stairstransp = 0.4

ax_interdiv_naive = Axis(fig_naive[1,1], xlabel="Time", ylabel="Probability density", title="Interdivision time distribution")
#ax_birth = Axis(fig[1,2], xlabel="Protein counts", ylabel="Probability density", title="Birth distribution")
ax_div_naive = Axis(fig_naive[1,2], xlabel="Protein count at division", ylabel="Probability density", title="Division protein distribution")
ax_joint₋Sm_naive = Axis(fig_naive[1,3], xlabel="Cell age at division", ylabel="Protein count at division", title="Division distribution")
ax_joint₊Sm_naive = Axis(fig_naive[1,4], xlabel="Cell age at division", ylabel="Protein count at division", title="Division distribution")

stairs!(ax_div_naive, collect(midpoints(div_hist₋Sm.edges[1])), div_hist₋Sm.weights; 
    color=(colors[2], stairstransp), step=:center)
barplot!(ax_div_naive, collect(midpoints(div_hist₋Sm.edges[1])), div_hist₋Sm.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_div_naive, xs, analyticals₋Sm.results[:division_dist_ancest]; 
    color=colors[2], linewidth=3.0)

stairs!(ax_div_naive, collect(midpoints(div_hist₊Sm.edges[1])), div_hist₊Sm.weights; 
    color=(colors[1], stairstransp), step=:center)
barplot!(ax_div_naive, collect(midpoints(div_hist₊Sm.edges[1])), div_hist₊Sm.weights; 
    color=(colors[1], transp), strokecolor=(colors[1], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_div_naive, xs, analyticals₊Sm_naive.results[:division_dist_ancest]; 
    color=colors[1], linewidth=3.0)
xlims!(ax_div_naive, (50, 200))
ylims!(ax_div_naive, low=0.0)

stairs!(ax_interdiv_naive, collect(midpoints(interdiv_hist₋Sm.edges[1])), interdiv_hist₋Sm.weights; 
    color=(colors[2], stairstransp), step=:center)
barplot!(ax_interdiv_naive, collect(midpoints(interdiv_hist₋Sm.edges[1])), interdiv_hist₋Sm.weights; 
    color=(colors[2], transp), strokecolor=(colors[2], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_interdiv_naive, ts, div_dist₋Sm.(ts); 
    color=colors[2], linewidth=3.0)

stairs!(ax_interdiv_naive, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; 
    color=(colors[1], stairstransp), step=:center)
barplot!(ax_interdiv_naive, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; 
    color=(colors[1], transp), strokecolor=(colors[1], transp), strokewidth=0.0, gap=0.0, dodge_gap=0.0)
lines!(ax_interdiv_naive, ts, div_dist₊Sm_naive.(ts); 
    color=colors[1], linewidth=3.0)
xlims!(ax_interdiv_naive, (25, 150))
ylims!(ax_interdiv_naive, low=0.0)

zs₋Sm = hcat(fpt_dist_ancest(analyticals₋Sm).(0.0:1.0:125.0)...)
cmap1 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 0.8), length=15)
contourf!(ax_joint₋Sm_naive, hist₋Sm; colormap=cmap1, levels=5)
cmap1 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 1.0), length=15)
contour!(ax_joint₋Sm_naive, 0.0:1.0:125.0, 1:length(analyticals₋Sm.results[:birth_dist]), zs₋Sm'; colormap=cmap1, linewidth=3, levels=4)

zs₊Sm_naive = hcat(fpt_dist_ancest(analyticals₊Sm_naive).(0.0:1.0:125.0)...)
cmap2 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 0.8), length=15)
contourf!(ax_joint₊Sm_naive, hist₊Sm; colormap=cmap2, levels=5)
cmap2 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 1.0), length=15)
contour!(ax_joint₊Sm_naive, 0.0:1.0:125.0, 1:length(analyticals₊Sm_naive.results[:birth_dist]), zs₊Sm_naive'; colormap=cmap2, linewidth=2, levels=4)


xlims!(ax_joint₋Sm_naive, (25.0, 125.0))
ylims!(ax_joint₋Sm_naive, (50.0, 150.0))

xlims!(ax_joint₊Sm_naive, (25.0, 125.0))
ylims!(ax_joint₊Sm_naive, (50.0, 150.0))

Legend(fig_naive[2,:], [data_legend, fitting_legend], [["No treatment", "Treatment"], ["No treatment", "No selection"]], ["Data", "Model fitting"]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true, groupgap=60)

save("$(plotsdir())/Wakamoto/apx_fpts_wk.pdf", fig_naive, pt_per_unit = 1)
save("$(plotsdir())/Wakamoto/apx_fpts_wk.svg", fig_naive, pt_per_unit = 1)

#birth_div = collect.(collect(zip(analyticals₋Sm.results[:birth_dist], analyticals₋Sm.results[:division_dist])))
#birth_div = reduce(hcat, vec.(birth_div))'
#CSV.write("$(datadir())/exp_pro/Wakamoto/F3NW-Sm/birth_div_theory_noSm.csv",  
#        Tables.table(birth_div), writeheader=false)
#CSV.write("$(datadir())/exp_pro/Wakamoto/F3NW-Sm/pdf_interdiv_noSm.csv",
#          Tables.table(div_dist₋Sm.(ts)), writeheader=false)


#birth_div₊Sm = collect.(collect(zip(analyticals_mult₊Sm.results[:birth_dist], analyticals_mult₊Sm.results[:division_dist])))
#birth_div₊Sm = reduce(hcat, vec.(birth_div₊Sm))'
#CSV.write("$(datadir())/exp_pro/Wakamoto/F3NW+Sm/birth_div_theory_yesSm_reaction_rates_changed.csv",  
#        Tables.table(birth_div₊Sm), writeheader=false)
#CSV.write("$(datadir())/exp_pro/Wakamoto/F3NW+Sm/pdf_interdiv_yesSm_reaction_rates_changed.csv",
#          Tables.table(div_dist₊Sm.(ts)), writeheader=false)
