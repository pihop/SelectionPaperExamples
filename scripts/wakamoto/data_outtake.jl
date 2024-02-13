#GLMakie.activate!()
#size_inches = (16.2, 11.2)
#size_pt = 72 .* size_inches
#figure_hitting = Figure(resolution=size_pt, fontsize=12, figure_padding=(2, 65, 2, 2))
#ax_hitting = Axis(figure_hitting[2,1], xlabel="Cell age", ylabel="Protein counts", title="First passage time distribution")
#ax_hist_fluor = Axis(figure_hitting[2, 2], xlabel="Probability density", xticks=0.0:0.01:0.02, title="Division protein distribution")
#ax_hist_div = Axis(figure_hitting[1, 1], ylabel="Probability density", xticks=0.0:50:150, title="Division age distribution")
#
#cycles₋Sm = vcat(make_cell_cycle.(lineages₋Sm)...); 
#cycles₊Sm = vcat(make_cell_cycle.(lineages₊Sm)...); 
#
#fluors₋Sm = fluorescence_trajectory.(cycles₋Sm)
#fluors₊Sm = fluorescence_trajectory.(cycles₊Sm)
#
#hitting₋Sm = hcat([[length(traj) * 5, traj[end] ./ cf₋Sm] for traj in fluors₋Sm]...)
#hitting₊Sm = hcat([[length(traj) * 5, traj[end] ./ cf₋Sm] for traj in fluors₊Sm]...)
#
#hist₋Sm = fit(Histogram, Tuple(eachrow(hitting₋Sm)), nbins=50)
#hist₊Sm = fit(Histogram, Tuple(eachrow(hitting₊Sm)), nbins=50)
#
#kde₋Sm = kde(Tuple(eachrow(hitting₋Sm)))
#kde₊Sm = kde(Tuple(eachrow(hitting₊Sm)))
#
#cmap1 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 0.8), length=15)
#contourf!(ax_hitting, kde₋Sm; colormap=cmap1)
#
#cmap2 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 0.8), length=15)
#contourf!(ax_hitting, kde₊Sm; colormap=cmap2)
#
#xlims!(ax_hitting, (0.0, 125.0))
#ylims!(ax_hitting, (30.0, 170.0))
#
#idxs₋Sm = rand(1:length(fluors₋Sm), 10)
#idxs₊Sm = rand(1:length(fluors₊Sm), 10)
#
##for idx in idxs₋Sm 
##    lines!(ax_hitting, collect(1:length(fluors₋Sm[idx])) .* 5, fluors₋Sm[idx] ./ cf₋Sm; color=(colors[1], 0.5), linewidth=3)
##    scatter!(ax_hitting, length(fluors₋Sm[idx]) * 5, fluors₋Sm[idx][end] ./ cf₋Sm; color=(colors[1], 0.8), strokewidth=2)
##end
#
#for idx in idxs₊Sm 
#    lines!(ax_hitting, collect(1:length(fluors₊Sm[idx])) .* 5, fluors₊Sm[idx] ./ cf₋Sm; color=(colors[2], 0.7), linewidth=2)
#    scatter!(ax_hitting, length(fluors₊Sm[idx]) * 5, fluors₊Sm[idx][end] ./ cf₋Sm; color=(colors[2], 0.8), strokewidth=2, marker=:xcross)
#end
#
##barplot!(ax_hist_div, collect(midpoints(interdiv_hist₋Sm.edges[1])), interdiv_hist₋Sm.weights; color=(colors[1], transp), gap=0.0)
##stairs!(ax_hist_div, collect(midpoints(interdiv_hist₋Sm.edges[1])), interdiv_hist₋Sm.weights; color=colors[1], step=:center)
#
#barplot!(ax_hist_div, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; color=(colors[2], transp), gap=0.0)
#stairs!(ax_hist_div, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; color=colors[2], step=:center)
#xlims!(ax_hist_div, (0.0, 125.0))
#ylims!(ax_hist_div, low=0.0)
#
##barplot!(ax_hist_fluor, collect(midpoints(divn_hist₋Sm.edges[1])), divn_hist₋Sm.weights; color=(colors[1], transp), gap=0.0, direction=:x)
##stairs!(ax_hist_fluor, divn_hist₋Sm.weights, collect(midpoints(divn_hist₋Sm.edges[1])) .- 0.5; color=colors[1])
#
#barplot!(ax_hist_fluor, collect(midpoints(divn_hist₊Sm.edges[1])), divn_hist₊Sm.weights; color=(colors[2], transp), gap=0.0, direction=:x)
#stairs!(ax_hist_fluor, divn_hist₊Sm.weights, collect(midpoints(divn_hist₊Sm.edges[1])) .- 0.5; color=colors[2])
#
#ylims!(ax_hist_fluor, (30.0, 170.0))
#xlims!(ax_hist_fluor, low=0.0)
#
#colsize!(figure_hitting.layout, 1, Relative(3.2 / 5))
#rowsize!(figure_hitting.layout, 1, Relative(1.5 / 5))
#
#hideydecorations!(ax_hist_fluor, ticks=false, grid=false)
#hidexdecorations!(ax_hist_div, ticks=false, grid=false)
#
##xticks!(ax_hist_div, xtickranges=0.0:0.1:0.2)
#
#colgap!(figure_hitting.layout, 70)
#rowgap!(figure_hitting.layout, 10)
#
#CairoMakie.activate!()
#save("$(plotsdir())/wakamoto/first_passage_$strain.pdf", figure_hitting)
#
#ratios₋Sm = []
#for t in 5:15
#    valsa = filter(x -> length(x) >= t, fluors₋Sm)
#    valsb = filter(x -> length(x) == t, fluors₋Sm)
#    a = getindex.(valsa, t)
#    b = getindex.(valsb, t)
#    r = densratiofunc(a ./ cf₋Sm, b ./ cf₋Sm, LSIF(10, 10, 0.001, rng))
#    push!(ratios₋Sm, r)
#end
#
#ratios₊Sm = [] 
#for t in 5:15
#    valsa = filter(x -> length(x) >= t, fluors₊Sm)
#    valsb = filter(x -> length(x) == t, fluors₊Sm)
#    a = getindex.(valsa, t)
#    b = getindex.(valsb, t)
#    r = densratiofunc(a ./ cf₋Sm, b ./ cf₋Sm, LSIF(10, 10, 0.001, rng))
#    push!(ratios₊Sm, r)
#end
#fig = Figure()
#ax = Axis(fig[1,1])
#
#trange = 1:5:150
#lines!(ax, trange, ratios₋Sm[10].(trange); color=colors[1])
#lines!(ax, trange, ratios₊Sm[10].(trange); color=colors[2])
#
#save("$(plotsdir())/wakamoto/test_density_ratio.pdf", fig)

#function selection_dists_wasser(age, cycles, cf)
#    bars_cycles = fluorescence_trajectory.(filter(x -> size(x)[1] == age, cycles)) ./ cf
#    nobars_cycles = fluorescence_trajectory.(filter(x -> size(x)[1] >= age, cycles)) ./ cf
#
#    bars_fluor = Float64.(getindex.(bars_cycles, age))
#    nobars_fluor = 1 ./ Float64.(getindex.(nobars_cycles, age))
#
#    bars_hist = normalize(fit(Histogram, getindex.(bars_cycles, age), 50:10:160); mode=:probability)
#    nobars_hist = normalize(fit(Histogram, getindex.(nobars_cycles, age), 50:10:160); mode=:probability)
#
#    dbars = DiscreteNonParametric(midpoints(bars_hist.edges[1]), bars_hist.weights)
#    dnobars = DiscreteNonParametric(midpoints(nobars_hist.edges[1]), nobars_hist.weights)
#
#    return midpoints(bars_hist.edges[1]), wasserstein(dnobars, dbars) 
#    #return midpoints(bars_hist.edges[1]), map(x -> isnan(x) ? missing : x,  bars_hist.weights ./ nobars_hist.weights)
#end
#
#function selection_dists(age, cycles, cf)
#    bars_cycles = fluorescence_trajectory.(filter(x -> size(x)[1] == age, cycles)) ./ cf
#    nobars_cycles = fluorescence_trajectory.(filter(x -> size(x)[1] >= age, cycles)) ./ cf
#
#    bars_fluor = Float64.(getindex.(bars_cycles, age))
#    nobars_fluor = 1 ./ Float64.(getindex.(nobars_cycles, age))
#
#    bars_hist = normalize(fit(Histogram, getindex.(bars_cycles, age), 70:10:120); mode=:probability)
#    nobars_hist = normalize(fit(Histogram, getindex.(nobars_cycles, age), 70:10:120); mode=:probability)
#
#    dbars = DiscreteNonParametric(midpoints(bars_hist.edges[1]), bars_hist.weights)
#    dnobars = DiscreteNonParametric(midpoints(nobars_hist.edges[1]), nobars_hist.weights)
#
#    return midpoints(bars_hist.edges[1]), map(x -> isnan(x) ? missing : x,  bars_hist.weights ./ nobars_hist.weights)
#end
#
#
## Compute the age distribution.
##ts = 4:20
#ts₊Sm = [x.first for x in filter(x -> x.second > 100, countmap(getindex.(size.(cycles₊Sm), 1)))]
#ts₋Sm = [x.first for x in filter(x -> x.second > 100, countmap(getindex.(size.(cycles₋Sm), 1)))]
#ts = intersect(ts₊Sm, ts₋Sm)
#
#age_dist₋Sm = [length(filter(x -> size(x)[1] >= age, cycles₋Sm)) for age in ts]
#age_dist₋Sm = age_dist₋Sm ./ sum(age_dist₋Sm)
##age_dist₋Sm = age_dist₋Sm ./ sum(getindex.(size.(cycles₋Sm), 1))
#age_dist₊Sm = [length(filter(x -> size(x)[1] >= age, cycles₊Sm)) for age in ts]
#age_dist₊Sm = age_dist₊Sm ./ sum(age_dist₊Sm)
##age_dist₊Sm = age_dist₊Sm ./ sum(getindex.(size.(cycles₊Sm), 1))
#
#sel₋Sm = [selection_dists_wasser(age, cycles₋Sm, cf₋Sm) for age in ts]
#sel₊Sm = [selection_dists_wasser(age, cycles₊Sm, cf₋Sm) for age in ts]
#
#wassmean₋Sm = Impute.interp(vec(sum(hcat(getindex.(sel₋Sm, 2)...)' .* age_dist₋Sm, dims=1)))
#wassmean₊Sm = Impute.interp(vec(sum(hcat(getindex.(sel₊Sm, 2)...)' .* age_dist₊Sm, dims=1)))
#
#wassvar₋Sm = sqrt.(Impute.interp(vec(sum((hcat(getindex.(sel₋Sm, 2)...)' .- mean₋Sm') .^2 .* age_dist₋Sm, dims=1))))
#wassvar₊Sm = sqrt.(Impute.interp(vec(sum((hcat(getindex.(sel₊Sm, 2)...)' .- mean₊Sm') .^2 .* age_dist₊Sm, dims=1))))
#
#sel₋Sm = [selection_dists(age, cycles₋Sm, cf₋Sm) for age in ts]
#sel₊Sm = [selection_dists(age, cycles₊Sm, cf₋Sm) for age in ts]
#
#mean₋Sm = Impute.interp(vec(sum(hcat(getindex.(sel₋Sm, 2)...)' .* age_dist₋Sm, dims=1)))
#mean₊Sm = Impute.interp(vec(sum(hcat(getindex.(sel₊Sm, 2)...)' .* age_dist₊Sm, dims=1)))
#
#var₋Sm = sqrt.(Impute.interp(vec(sum((hcat(getindex.(sel₋Sm, 2)...)' .- mean₋Sm') .^2 .* age_dist₋Sm, dims=1))))
#var₊Sm = sqrt.(Impute.interp(vec(sum((hcat(getindex.(sel₊Sm, 2)...)' .- mean₊Sm') .^2 .* age_dist₊Sm, dims=1))))
#
#figure_selection_wk = Figure()
#ax_traj_wk = Axis(figure_selection_wk[1,1];)
#ax_hist_wk = Axis(figure_selection_wk[2,1];)
#
#rangebars!(ax_traj_wk, sel₋Sm[1][1], mean₋Sm .- 2*var₋Sm, mean₋Sm .+ 2*var₋Sm; color=colors[1], whiskerwidth = 10)
#rangebars!(ax_traj_wk, sel₋Sm[1][1], mean₊Sm .- 2*var₊Sm, mean₊Sm .+ 2*var₊Sm; color=colors[2], whiskerwidth = 10)
#hlines!(ax_traj_wk, 1.0)
#ylims!(ax_traj_wk, (-1.0, 3.0))
#
##GLMakie.activate!()
#
##ts₊Sm = [x.first for x in filter(x -> x.second > 400, countmap(getindex.(size.(cycles₊Sm), 1)))]
##ts₋Sm = [x.first for x in filter(x -> x.second > 400, countmap(getindex.(size.(cycles₋Sm), 1)))]
#
##for (i, t) in enumerate(ts) #enumerate(intersect(ts₋Sm, ts₊Sm))
##    k₋Smbar, k₋Smnobar = selection_dists(t, cycles₋Sm, cf₋Sm)
##    k₊Smbar, k₊Smnobar = selection_dists(t, cycles₊Sm, cf₋Sm)
##
###    mean₋Smbar = sum(midpoints(k₋Smbar.edges[1]) .* k₋Smbar.weights) ./ sum(midpoints(k₋Smnobar.edges[1]) .* k₋Smnobar.weights)
###    mean₊Smbar = sum(midpoints(k₊Smbar.edges[1]) .* k₊Smbar.weights) ./ sum(midpoints(k₊Smnobar.edges[1]) .* k₊Smnobar.weights)
##
###    rangebars!(ax_traj_wk, t, [k₋Smbar[1] - sqrt(k₋Smbar[2]),], [k₋Smbar[1] + sqrt(k₋Smbar[2]),]; color=colors[1])
###    rangebars!(ax_traj_wk, t, [k₊Smbar[1] - sqrt(k₊Smbar[2]),], [k₊Smbar[1] + sqrt(k₊Smbar[2]),]; color=colors[2])
##
##
###    lines!(ax_traj_wk, 0.0:1.0:150, t->k₋Smbar(t)/k₋Smnobar(t); color=colors[1])
###    lines!(ax_traj_wk, 0.0:1.0:150, t->k₊Smbar(t)/k₊Smnobar(t); color=colors[2])
###
####    scatter!(ax_traj_wk, t, kldivergence(k₋Smbar.weights, k₋Smnobar.weights); color=colors[1])
####    scatter!(ax_traj_wk, t, kldivergence(k₊Smbar.weights, k₊Smnobar.weights); color=colors[2])
##    scatterlines!(ax_traj_wk, midpoints(k₋Smbar.edges[1]), k₋Smbar.weights ./ k₋Smnobar.weights .-1; color=colors[1])
##    scatterlines!(ax_traj_wk, midpoints(k₊Smbar.edges[1]), k₊Smbar.weights ./ k₊Smnobar.weights .-1; color=colors[2], linestyle=:dash)
##end
#
#
##idx = 9
##ylims!(ax_traj_wk, (0.0, 2.0))
##k1, k2 = selection_dists(idx, cycles₋Sm, cf₋Sm)
##lines!(ax_traj_wk, 0.0:1.0:150, t->k1(t)/k2(t); color=colors[1])
##k1, k2 = selection_dists(idx, cycles₊Sm, cf₋Sm)
##lines!(ax_traj_wk, 0.0:1.0:150, t->k1(t)/k2(t); color=colors[2])
#
##lines!(ax_traj_wk, 0.0:1.0:150, k2; color=colors[2])
#
##stairs!(ax_hist_wk, collect(midpoints(k1.edges[1])), k1.weights; color=colors[1], step=:center) 
##barplot!(ax_hist_wk, collect(midpoints(k1.edges[1])), k1.weights; color=(colors[1], transp), gap=0.0) 
##stairs!(ax_hist_wk, collect(midpoints(k2.edges[1])), k2.weights; color=colors[2], step=:center) 
##barplot!(ax_hist_wk, collect(midpoints(k2.edges[1])), k2.weights; color=(colors[2], transp), gap=0.0) 
##stairs!(ax_hist_wk, collect(midpoints(divn_hist₊Sm.edges[1])), divn_hist₊Sm.weights; color=colors[2], step=:center)
##barplot!(ax_hist_wk, collect(midpoints(divn_hist₊Sm.edges[1])), divn_hist₊Sm.weights; color=(colors[2], transp), gap=0.0)
##stairs!(ax_hist_wk, collect(midpoints(divn_hist₋Sm.edges[1])), divn_hist₋Sm.weights; color=colors[1], step=:center)
##barplot!(ax_hist_wk, collect(midpoints(divn_hist₋Sm.edges[1])), divn_hist₋Sm.weights; color=(colors[1], transp), gap=0.0)
##xlims!(ax_hist_wk, (15,150))
##xlims!(ax_traj_wk, (15,150))
##ylims!(ax_traj_wk, (0.0, 3.0))
##ylims!(ax_hist_wk, low=0)
#
#CairoMakie.activate!()
#save("$(plotsdir())/wakamoto/selection_$strain.pdf", figure_selection_wk)


