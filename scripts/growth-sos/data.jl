using DataFrames
using CSV
using GLMakie
using Glob
using StatsBase
using Colors
using ColorSchemes
using CairoMakie
using LinearAlgebra
using KernelDensity
using Distributions
using JLD2

include("utils.jl")


function birth_div_dist(cf, md_fluors, str)
    birth_fluor = collect.(collect(zip(getindex.(md_fluors, 2), getindex.(md_fluors, 2))))
    birth_fluor = reduce(hcat, vec.(birth_fluor))'
    birth_fluor[:, 2] .= birth_fluor[:, 2] ./ cf
    div_joint = collect.(collect(zip(getindex.(md_fluors, 1), getindex.(md_fluors, 3), getindex.(md_fluors, 3))))

    mkpath("$(datadir())/exp_pro/growth_sos/$str/")
    CSV.write("$(datadir())/exp_pro/growth_sos/$str/birth.csv",  
        Tables.table(birth_fluor), writeheader=false)

    div_joint = reduce(hcat, vec.(div_joint))'
    div_joint[:, 3] .= div_joint[:, 3] ./ cf

    CSV.write("$(datadir())/exp_pro/growth_sos/$str/div.csv",  
        Tables.table(div_joint), writeheader=false)

#    CSV.write("$(datadir())/exp_pro/growth_sos/$str/div.csv",  
#        Tables.table(div_joint), writeheader=false)

    return birth_fluor[:,1], birth_fluor[:, 2], div_joint[:,2], div_joint[:,3], div_joint[:, 1]
end

strainWT = "WT"
strain2pal = "2pal"
env = "Glu"
tslice = 8

# Input files.
dirs_wt = glob("*_$(strainWT)_$(env)_Analysis_tables/", "$(datadir())/exp_raw/growth-sos-msb-2021/mmc_segmentation-and-lineage-data/Data tables/$strainWT/")
dirs_2pal = glob("*_$(strain2pal)_$(env)_Analysis_tables/", "$(datadir())/exp_raw/growth-sos-msb-2021/mmc_segmentation-and-lineage-data/Data tables/$strain2pal/")

files_lengths_wt = vcat([glob("*_areamats.csv", dir) for dir in dirs_wt]...)
files_lengths_2pal = vcat([glob("*_areamats.csv", dir) for dir in dirs_2pal]...)

lengths_wt = Matrix.(DataFrame.(CSV.File.(files_lengths_wt, header=false);))
lengths_2pal = Matrix.(DataFrame.(CSV.File.(files_lengths_2pal, header=false);))

files_gfp_wt = vcat([glob("*_gfpmat_scaled.csv", dir) for dir in dirs_wt]...)
files_gfp_2pal = vcat([glob("*_gfpmat_scaled.csv", dir) for dir in dirs_2pal]...)

gfpmat_areanorm_wt = Matrix.(DataFrame.(CSV.File.(files_gfp_wt, header=false);))
gfpmat_areanorm_2pal = Matrix.(DataFrame.(CSV.File.(files_gfp_2pal, header=false);))

gfpmat_wt = make_gfp_matrix.(lengths_wt, gfpmat_areanorm_wt)
gfpmat_2pal = make_gfp_matrix.(lengths_2pal, gfpmat_areanorm_2pal)

#div_times_wt = division_times.(lengths_wt)
#div_times_2pal = division_times.(lengths_2pal)

div_times_wt = div_time_analysis.(lengths_wt; window=0, threshold=(-0.7, -0.3))
div_times_2pal = div_time_analysis.(lengths_2pal; window=0, threshold=(-0.7, -0.3)) 

md_div_times_wt = getindex.(div_times_wt, 1)
bm_div_times_wt = getindex.(div_times_wt, 2)

md_div_times_2pal = getindex.(div_times_2pal, 1)
bm_div_times_2pal = getindex.(div_times_2pal, 2)

#fluor_div_wt = vcat([fluorecence_at_division(div_idxs, gfp) for (div_idxs, gfp) in zip(getindex.(div_times_wt, 2), wt_gfp)]...)
#fluor_div_del = vcat([fluorecence_at_division(div_idxs, gfp) for (div_idxs, gfp) in zip(getindex.(div_times_del, 2), del_gfp)]...)

md_fluor_wt = vcat([md_fluorescence(div_idxs, gfp; ftime=(1, 100)) for (div_idxs, gfp) in zip(md_div_times_wt, gfpmat_wt)]...)
md_fluor_2pal = vcat([md_fluorescence(div_idxs, gfp; ftime=(1, 100)) for (div_idxs, gfp) in zip(md_div_times_2pal, gfpmat_2pal)]...)

divfluor_time_wt = vcat([birthdiv_fluorescence_time(div_idxs, gfp; ftime=(1, 100)) for (div_idxs, gfp) in zip(bm_div_times_wt, gfpmat_wt)]...)
divfluor_time_2pal = vcat([birthdiv_fluorescence_time(div_idxs, gfp; ftime=(1, 100)) for (div_idxs, gfp) in zip(bm_div_times_2pal, gfpmat_2pal)]...)

df_wt, ols_wt, cf_wt = conversion_fact_regression(md_fluor_wt; fltr=(1.0, 3.5))
sort!(df_wt, [:X])
df_2pal, ols_2pal, cf_2pal = conversion_fact_regression(md_fluor_2pal; fltr=(1.0, 4.5))
sort!(df_2pal, [:X])

birth_dist_wt, birthn_dist_wt, div_dist_wt, divn_dist_wt, interdiv_times_wt = birth_div_dist(cf_wt, divfluor_time_wt, "wt")
birth_dist_2pal, birthn_dist_2pal, div_dist_2pal, divn_dist_2pal, interdiv_times_2pal = birth_div_dist(cf_wt, divfluor_time_2pal, "2pal")

interdiv_dist_wt = interdiv_dist_kde(interdiv_times_wt; time_slice=tslice)
interdiv_dist_2pal = interdiv_dist_kde(interdiv_times_2pal; time_slice=tslice)

#fluor_div_wt = vcat([fluorecence_at_division(div_idxs, gfp) for (div_idxs, gfp) in zip(getindex.(div_times_wt, 2), gfpmat_wt)]...)
#fluor_div_2pal = vcat([fluorecence_at_division(div_idxs, gfp) for (div_idxs, gfp) in zip(getindex.(div_times_2pal, 2), gfpmat_2pal)]...)
#
#md_fluor_wt = vcat([fluorescence_tuples(div_idxs, gfp) for (div_idxs, gfp) in zip(div_times_wt, gfpmat_wt)]...)
#md_fluor_2pal = vcat([fluorescence_tuples(div_idxs, gfp) for (div_idxs, gfp) in zip(div_times_2pal, gfpmat_2pal)]...)

colors = ColorSchemes.Egypt.colors
transp = 0.5
# Plot
size_inches = (17, 18)
size_pt = 72 .* size_inches
fig = Figure(resolution=size_pt, fontsize=18)

ax_mdfluor = Axis(fig[1,1]; xlabel="Mother cell fluorescence", ylabel="Daugther cell \n fluorescence")
ax_mdfluor2 = Axis(fig[1,2]; xlabel="Mother cell fluorescence", ylabel="Daugther cell \n fluorescence")
#ax_divtime = Axis(fig[1,2]; xlabel="Division time", ylabel="Probability density")
ax_var = Axis(fig[2,1:2]; xlabel="Mother cell fluorescence", ylabel="Variance of daughter \n cell fluorescence")
ax_birth = Axis(fig[3,1]; xlabel="Protein fluorescence", ylabel="Probability density", title="Birth fluorescence")
ax_div = Axis(fig[3,2]; xlabel="Protein fluorescence", ylabel="Probability density", title="Division fluorescence")

ax_birthn = Axis(fig[4,1]; xlabel="Protein numbers", ylabel="Probability density", title="Birth protein numbers")
ax_divn = Axis(fig[4,2]; xlabel="Protein numbers", ylabel="Probability density", title="Division protein numbers")

#lines!(ax, 1:length(lengths_wt[1][1,:]), lengths_wt[1][1,:])
#scatter!(ax, 1:length(lengths_wt[1][1,:]), lengths_wt[1][1,:])
#scatter!(ax, div_times_wt[1][2][1], lengths_wt[1][1,:][div_times_wt[1][2][1]])
#scatter!(ax, div_times_wt[1][4][1], lengths_wt[1][1,:][div_times_wt[1][4][1]])

#hist!(ax2, fluor_div_wt; bins = 50, normalization=:pdf, color=colors[1])
#hist!(ax2, fluor_div_2pal; bins = 50, normalization=:pdf, color=(colors[3], 0.8))
#xlims!(ax2, (0.0, 20))

scatter!(ax_mdfluor, getindex.(md_fluor_wt, 1), getindex.(md_fluor_wt, 2); 
    color=(colors[1], transp), strokewidth=1.0, strokecolor = (colors[1], 0.9), label="Wild type")
scatter!(ax_mdfluor2, getindex.(md_fluor_2pal, 1), getindex.(md_fluor_2pal, 2); 
    color=(colors[2], transp), strokewidth=1.0, strokecolor = (colors[2], 0.9), label="2pal")
str_line(x) = 0.5*x
lines!(ax_mdfluor, 0.0:0.1:20.0, str_line.(0.0:0.1:20.0); color=(:gray, 0.8), linewidth=3.0)
lines!(ax_mdfluor2, 0.0:0.1:20.0, str_line.(0.0:0.1:20.0); color=(:gray, 0.8), linewidth=3.0)
xlims!(ax_mdfluor, (0, 10))
xlims!(ax_mdfluor2, (0, 10))
ylims!(ax_mdfluor2, (0, 10))
ylims!(ax_mdfluor, (0, 5))

scatter!(ax_var, df_wt[!,:X], df_wt[!,:Z]; color=colors[1])
errorbars!(ax_var, df_wt[!,:X], df_wt[!,:Z], 1.96.*sqrt.(df_wt[!,:W]); color=colors[1])
lines!(ax_var, df_wt[!,:X], predict(ols_wt, df_wt); color=colors[1])

scatter!(ax_var, df_2pal[!,:X], df_2pal[!,:Z]; color=colors[2])
errorbars!(ax_var, df_2pal[!,:X], df_2pal[!,:Z], 1.96.*sqrt.(df_2pal[!,:W]); color=colors[2])
lines!(ax_var, df_2pal[!,:X], predict(ols_2pal, df_2pal); color=colors[2])

tspan = collect(0.0:1.0:150.0)
interdiv_hist_wt = normalize(fit(Histogram, interdiv_times_wt .* 8, 1:8:maximum(interdiv_times_wt .* 9)); mode=:pdf)
stairs!(ax_divtime, collect(midpoints(interdiv_hist_wt.edges[1])), interdiv_hist_wt.weights; color=colors[1], step=:center)
barplot!(ax_divtime, collect(midpoints(interdiv_hist_wt.edges[1])), interdiv_hist_wt.weights; color=(colors[1], transp), gap=0.0)
dist_wt = lines!(ax_divtime, tspan, pdf.(Ref(interdiv_dist_wt), tspan); color=colors[1], linewidth=2, linestyle=:dash)

interdiv_hist_2pal = normalize(fit(Histogram, interdiv_times_2pal .* 8, 1:8:maximum(interdiv_times_2pal .* 8)); mode=:pdf)
stairs!(ax_divtime, collect(midpoints(interdiv_hist_2pal.edges[1])), interdiv_hist_2pal.weights; color=colors[2], step=:center)
barplot!(ax_divtime, collect(midpoints(interdiv_hist_2pal.edges[1])), interdiv_hist_2pal.weights; color=(colors[2], transp), gap=0.0)
dist_2pal = lines!(ax_divtime, tspan, pdf.(Ref(interdiv_dist_2pal), tspan); color=colors[2], linewidth=2, linestyle=:dash)
xlims!(ax_divtime, (0, 250))
#axislegend(ax_divtime, [LineElement(color = (:black, 0.8), linestyle = :dash),], ["Kernel density estimate", ], framevisible = false)

birth_hist_wt = normalize(fit(Histogram, birth_dist_wt, 0:0.1:4); mode=:pdf)
stairs!(ax_birth, collect(midpoints(birth_hist_wt.edges[1])), birth_hist_wt.weights; color=colors[1], step=:center)
barplot!(ax_birth, collect(midpoints(birth_hist_wt.edges[1])), birth_hist_wt.weights; color=(colors[1], transp), gap=0.0)
birth_hist_2pal = normalize(fit(Histogram, birth_dist_2pal, 0:0.1:4); mode=:pdf)
stairs!(ax_birth, collect(midpoints(birth_hist_2pal.edges[1])), birth_hist_2pal.weights; color=colors[2], step=:center)
barplot!(ax_birth, collect(midpoints(birth_hist_2pal.edges[1])), birth_hist_2pal.weights; color=(colors[2], transp), gap=0.0)
xlims!(ax_birth, (0.0, 4.))

div_hist_wt = normalize(fit(Histogram, div_dist_wt, 1:0.2:7); mode=:pdf)
stairs!(ax_div, collect(midpoints(div_hist_wt.edges[1])), div_hist_wt.weights; color=colors[1], step=:center)
barplot!(ax_div, collect(midpoints(div_hist_wt.edges[1])), div_hist_wt.weights; color=(colors[1], transp), gap=0.0)
div_hist_2pal = normalize(fit(Histogram, div_dist_2pal, 1:0.2:7); mode=:pdf)
stairs!(ax_div, collect(midpoints(div_hist_2pal.edges[1])), div_hist_2pal.weights; color=colors[2], step=:center)
barplot!(ax_div, collect(midpoints(div_hist_2pal.edges[1])), div_hist_2pal.weights; color=(colors[2], transp), gap=0.0)
xlims!(ax_div, (1.0, 7.))

birthn_hist_wt = normalize(fit(Histogram, birthn_dist_wt, 1:200); mode=:pdf)
stairs!(ax_birthn, collect(midpoints(birthn_hist_wt.edges[1])), birthn_hist_wt.weights; color=colors[1], step=:center)
barplot!(ax_birthn, collect(midpoints(birthn_hist_wt.edges[1])), birthn_hist_wt.weights; color=(colors[1], transp), gap=0.0)
birthn_hist_2pal = normalize(fit(Histogram, birthn_dist_2pal, 1:200); mode=:pdf)
stairs!(ax_birthn, collect(midpoints(birthn_hist_2pal.edges[1])), birthn_hist_2pal.weights; color=colors[2], step=:center)
barplot!(ax_birthn, collect(midpoints(birthn_hist_2pal.edges[1])), birthn_hist_2pal.weights; color=(colors[2], transp), gap=0.0)
xlims!(ax_birthn, (10.0, 60.))

divn_hist_wt = normalize(fit(Histogram, divn_dist_wt, 1:200); mode=:pdf)
stairs!(ax_divn, collect(midpoints(divn_hist_wt.edges[1])), divn_hist_wt.weights; color=colors[1], step=:center)
barplot!(ax_divn, collect(midpoints(divn_hist_wt.edges[1])), divn_hist_wt.weights; color=(colors[1], transp), gap=0.0)
divn_hist_2pal = normalize(fit(Histogram, divn_dist_2pal, 1:200); mode=:pdf)
stairs!(ax_divn, collect(midpoints(divn_hist_2pal.edges[1])), divn_hist_2pal.weights; color=colors[2], step=:center)
barplot!(ax_divn, collect(midpoints(divn_hist_2pal.edges[1])), divn_hist_2pal.weights; color=(colors[2], transp), gap=0.0)

xlims!(ax_divn, (25.0, 125.))
#Legend(fig[5, 1:2], ax_mdfluor; orientation = :horizontal, tellwidth = false, tellheight = true, framevisible = false)
#
data_legend = [PolyElement(color=colors[1]), PolyElement(color=(colors[2]))]
Legend(fig[5,1:2], data_legend, ["Wild type", "DNA damage"]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true)

CairoMakie.activate!()
mkpath("$(plotsdir())/growth-sos")
save("$(plotsdir())/growth-sos/cf.pdf", fig)

GLMakie.activate!()
size_inches = (16.2, 11.2)
size_pt = 72 .* size_inches
figure_hitting = Figure(resolution=size_pt, fontsize=34, figure_padding=(2, 65, 2, 2))
ax_hitting = Axis(figure_hitting[2,1], xlabel="Cell age", ylabel="Protein counts", title="First passage time distribution")
ax_hist_fluor = Axis(figure_hitting[2, 2], xlabel="Probability density", xticks=0.0:0.01:0.03, title="Division protein distribution")
ax_hist_div = Axis(figure_hitting[1, 1], ylabel="Probability density", xticks=0.0:50:150, title="Division age distribution")

cycles_wt = vcat(make_cell_cycles.(gfpmat_wt, bm_div_times_wt)...)
cycles_2pal = vcat(make_cell_cycles.(gfpmat_2pal, bm_div_times_2pal)...)


hitting_wt = hcat([[length(traj) * tslice, traj[end] ./ cf_wt] for traj in cycles_wt]...)
hitting_2pal = hcat([[length(traj) * tslice, traj[end] ./ cf_wt] for traj in cycles_2pal]...)

hist_wt = fit(Histogram, Tuple(eachrow(hitting_wt)), nbins=150)
hist_2pal = fit(Histogram, Tuple(eachrow(hitting_2pal)), nbins=150)

kde_wt = kde(Tuple(eachrow(hitting_wt)))
kde_2pal = kde(Tuple(eachrow(hitting_2pal)))

cmap1 = range(alphacolor(colors[1], 0.0), stop=ARGB(colors[1], 0.8), length=15)
contourf!(ax_hitting, kde_wt; colormap=cmap1)

cmap2 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 0.8), length=15)
contourf!(ax_hitting, kde_2pal; colormap=cmap2)

xlims!(ax_hitting, (0.0, 200.0))
ylims!(ax_hitting, (0.0, 170.0))

idxs_wt = rand(1:length(cycles_wt), 10)
idxs_2pal = rand(1:length(cycles_2pal), 10)

#for idx in idxs_wt 
#    lines!(ax_hitting, collect(1:length(fluors_wt[idx])) .* 5, fluors_wt[idx] ./ cf_wt; color=(colors[1], 0.5), linewidth=3)
#    scatter!(ax_hitting, length(fluors_wt[idx]) * 5, fluors_wt[idx][end] ./ cf_wt; color=(colors[1], 0.8), strokewidth=2)
#end

for idx in idxs_2pal
    lines!(ax_hitting, collect(1:length(cycles_2pal[idx])) .* tslice, cycles_2pal[idx] ./ cf_wt; color=(colors[2], 0.7), linewidth=2)
    scatter!(ax_hitting, length(cycles_2pal[idx]) * tslice, cycles_2pal[idx][end] ./ cf_wt; color=(colors[2], 0.8), strokewidth=2, marker=:xcross)
end

#barplot!(ax_hist_div, collect(midpoints(interdiv_hist_wt.edges[1])), interdiv_hist_wt.weights; color=(colors[1], transp), gap=0.0)
#stairs!(ax_hist_div, collect(midpoints(interdiv_hist_wt.edges[1])), interdiv_hist_wt.weights; color=colors[1], step=:center)

barplot!(ax_hist_div, collect(midpoints(interdiv_hist_2pal.edges[1])), interdiv_hist_2pal.weights; color=(colors[2], transp), gap=0.0)
stairs!(ax_hist_div, collect(midpoints(interdiv_hist_2pal.edges[1])), interdiv_hist_2pal.weights; color=colors[2], step=:center)
xlims!(ax_hist_div, (0.0, 200.0))
ylims!(ax_hist_div, low=0.0)

#barplot!(ax_hist_fluor, collect(midpoints(divn_hist_wt.edges[1])), divn_hist_wt.weights; color=(colors[1], transp), gap=0.0, direction=:x)
#stairs!(ax_hist_fluor, divn_hist_wt.weights, collect(midpoints(divn_hist_wt.edges[1])) .- 0.5; color=colors[1])

barplot!(ax_hist_fluor, collect(midpoints(divn_hist_2pal.edges[1])), divn_hist_2pal.weights; color=(colors[2], transp), gap=0.0, direction=:x)
stairs!(ax_hist_fluor, divn_hist_2pal.weights, collect(midpoints(divn_hist_2pal.edges[1])) .- 0.5; color=colors[2])

ylims!(ax_hist_fluor, (0.0, 170.0))
xlims!(ax_hist_fluor, low=0.0)

colsize!(figure_hitting.layout, 1, Relative(3.2 / 5))
rowsize!(figure_hitting.layout, 1, Relative(1.5 / 5))

hideydecorations!(ax_hist_fluor, ticks=false, grid=false)
hidexdecorations!(ax_hist_div, ticks=false, grid=false)

#xticks!(ax_hist_div, xtickranges=0.0:0.1:0.2)

colgap!(figure_hitting.layout, 70)
rowgap!(figure_hitting.layout, 10)

CairoMakie.activate!()
save("$(plotsdir())/growth-sos/first_passage.pdf", figure_hitting)


cycles_wt = vcat(make_cell_cycles.(gfpmat_wt, bm_div_times_wt)...)
cycles_2pal = vcat(make_cell_cycles.(gfpmat_2pal, bm_div_times_2pal)...)

jldsave("$(datadir())/exp_pro/growth_sos/cycles.jld2"; cycles_wt, cycles_2pal, cf_wt, cf_2pal)
