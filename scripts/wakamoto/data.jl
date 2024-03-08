using DataFrames
using CSV
using Glob
using LinearAlgebra
using StatsBase
using KernelDensity
using Bootstrap
#using Impute
#using DensityRatioEstimation 

using ColorSchemes
#using GLMakie
using CairoMakie
using Colors
using JLD2

include("utils.jl")

strain = "F3NW"
#strain = "F3pTN001"
files₋Sm=glob("*/Results/*.xls", "$(datadir())/exp_raw/Wakamoto/ExperimentalData/$strain-Sm")
files₊Sm=glob("*/Results/*.xls", "$(datadir())/exp_raw/Wakamoto/ExperimentalData/$strain+Sm")

dfs₋Sm = DataFrame.(CSV.File.(files₋Sm))
dfs₊Sm = DataFrame.(CSV.File.(files₊Sm))

lineages₋Sm = compute_lineages.(dfs₋Sm; slice = [21, 121]);
lineages₊Sm = compute_lineages.(dfs₊Sm; slice = [21, 121]);

# Removes very short and long interdivision times. Ad hoc, is there a better
# gauge?
interdiv_times₋Sm = filter(x -> 1 < x, rmnothing(vcat(interdivision_times.(lineages₋Sm)...)))
interdiv_times₊Sm = filter(x -> 1 < x, rmnothing(vcat(interdivision_times.(lineages₊Sm)...)))
joint_fluor_times₋Sm = filter(x -> 1 < x[2], rmnothing(vcat(joint_fluor_times.(lineages₋Sm)...)))
joint_fluor_times₊Sm = filter(x -> 1 < x[2], rmnothing(vcat(joint_fluor_times.(lineages₊Sm)...)))

interdiv_times₋Sm = getindex.(joint_fluor_times₋Sm, 1)
interdiv_times₊Sm = getindex.(joint_fluor_times₊Sm, 1)

interdiv_dist₋Sm = interdiv_dist_gamma(interdiv_times₋Sm; time_slice=5)
interdiv_dist₊Sm = interdiv_dist_gamma(interdiv_times₊Sm; time_slice=5)

aggr_dist₋Sm = aggr_pdf(interdiv_dist₋Sm; rnge=range(0.0, step=5, stop=150))
aggr_dist₊Sm = aggr_pdf(interdiv_dist₊Sm; rnge=range(0.0, step=5, stop=150))

function birth_div_dist(cf, lineages, str)
    birth_fluor = rmnothing(vcat(birth_fluorescence.(lineages)...))
    birth_fluor = collect.(collect(zip(birth_fluor, birth_fluor)))
    birth_fluor = reduce(hcat, vec.(birth_fluor))'
    birth_fluor[:, 2] .= birth_fluor[:, 2] ./ cf

    div_joint = rmnothing(vcat(joint_fluor_times.(lineages)...))

    mkpath("$(datadir())/exp_pro/Wakamoto/$strain$str/")
    CSV.write("$(datadir())/exp_pro/Wakamoto/$strain$str/birth.csv",  
        Tables.table(birth_fluor), writeheader=false)

    div_joint = reduce(hcat, div_joint)'
    div_joint[:, 3] .= div_joint[:, 3] ./ cf

    CSV.write("$(datadir())/exp_pro/Wakamoto/$strain$str/div.csv",  
        Tables.table(div_joint), writeheader=false)

    jldsave("$(datadir())/exp_pro/Wakamoto/$strain$str/lineages.jld2"; lineages, cf)
    return birth_fluor[:, 2], div_joint[:,3], div_joint[:, 1]
end

# Compute the conversion factor from -Sm data and save the birth and death
# distributions as CSV files.
mother_daughter_fluor₋Sm = rmnothing(vcat(fluorescence_tuples.(lineages₋Sm)...))
df₋Sm, ols₋Sm, cf₋Sm = conversion_fact_regression(mother_daughter_fluor₋Sm)

mother_daughter_fluor₊Sm = rmnothing(vcat(fluorescence_tuples.(lineages₊Sm)...))
df₊Sm, ols₊Sm, cf₊Sm = conversion_fact_regression(mother_daughter_fluor₊Sm)

birth_dist₋Sm, div_dist₋Sm, interdiv_times₋Sm = birth_div_dist(cf₋Sm, lineages₋Sm, "-Sm")
#birth_dist₊Sm, div_dist₊Sm, interdiv_times₊Sm = birth_div_dist(cf₋Sm, lineages₊Sm, "+Sm")
birth_dist₊Sm, div_dist₊Sm, interdiv_times₊Sm = birth_div_dist(cf₋Sm, lineages₊Sm, "+Sm")

colors = ColorSchemes.Egypt.colors
transp = 0.5
# Plot
size_inches = (17, 18)
size_pt = 72 .* size_inches
fig = Figure(resolution=size_pt, fontsize=18, title="$strain")

ax_mdfluor = Axis(fig[1,1]; xlabel="Mother cell fluorescence", ylabel="Daugther cell fluorescence")
ax_mdfluor2 = Axis(fig[1,2]; xlabel="Mother cell fluorescence", ylabel="Daugther cell fluorescence")
#ax_divtime = Axis(fig[1,2]; xlabel="Division age", ylabel="Probability density")
ax_var = Axis(fig[2,1:2]; xlabel="Mother cell fluorescence", ylabel="Variance of daughter \n cell fluorescence")
ax_birth = Axis(fig[3,1]; xlabel="Protein fluorescence", ylabel="Probability density", title="Birth fluorescence")
ax_div = Axis(fig[3,2]; xlabel="Protein fluorescence", ylabel="Probability density", title="Division fluorescence")

ax_birthn = Axis(fig[4,1]; xlabel="Protein numbers", ylabel="Probability density", title="Birth protein numbers")
ax_divn = Axis(fig[4,2]; xlabel="Protein numbers", ylabel="Probability density", title="Division protein numbers")
xs = 0.0:1.0:30.0

scatter!(ax_mdfluor, getindex.(mother_daughter_fluor₋Sm, 1), getindex.(mother_daughter_fluor₋Sm, 2); 
    color=(colors[1], transp), strokewidth=1.0, strokecolor = (colors[1], 0.9), label="No treatment")
scatter!(ax_mdfluor2, getindex.(mother_daughter_fluor₊Sm, 1), getindex.(mother_daughter_fluor₊Sm, 2); 
    color=(colors[2], transp), strokewidth=1.0, strokecolor = (colors[2], 0.9), label="Treatment")
str_line(x) = 0.5*x
lines!(ax_mdfluor, 0.0:0.1:60.0, str_line.(0.0:0.1:60.0); color=(:gray, 0.8), linewidth=3.0)
lines!(ax_mdfluor2, 0.0:0.1:60.0, str_line.(0.0:0.1:60.0); color=(:gray, 0.8), linewidth=3.0)

scatter!(ax_var, df₋Sm[!,:X], df₋Sm[!,:Z]; color=colors[1])
errorbars!(ax_var, df₋Sm[!,:X], df₋Sm[!,:Z], 1.96.*sqrt.(df₋Sm[!,:W]); color=colors[1])
lines!(ax_var, df₋Sm[!,:X], StatsBase.predict(ols₋Sm, df₋Sm); color=colors[1])

scatter!(ax_var, df₊Sm[!,:X], df₊Sm[!,:Z]; color=colors[2])
errorbars!(ax_var, df₊Sm[!,:X], df₊Sm[!,:Z], 1.96.*sqrt.(df₊Sm[!,:W]); color=colors[2])
lines!(ax_var, df₊Sm[!,:X], StatsBase.predict(ols₊Sm, df₊Sm); color=colors[2])

tspan = collect(0.0:1.0:150.0)
#interdiv_hist₋Sm = normalize(fit(Histogram, interdiv_times₋Sm .* 5, 1:5:maximum(interdiv_times₋Sm .* 5)); mode=:pdf)
#stairs!(ax_divtime, collect(midpoints(interdiv_hist₋Sm.edges[1])), interdiv_hist₋Sm.weights; color=colors[1], step=:center)
#barplot!(ax_divtime, collect(midpoints(interdiv_hist₋Sm.edges[1])), interdiv_hist₋Sm.weights; color=(colors[1], transp), gap=0.0)
#dist₋Sm = lines!(ax_divtime, tspan, pdf.(interdiv_dist₋Sm, tspan); color=colors[1], linewidth=2, linestyle=:dash)
#
#interdiv_hist₊Sm = normalize(fit(Histogram, interdiv_times₊Sm .* 5, 1:5:maximum(interdiv_times₊Sm .* 5)); mode=:pdf)
#stairs!(ax_divtime, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; color=colors[2], step=:center)
#barplot!(ax_divtime, collect(midpoints(interdiv_hist₊Sm.edges[1])), interdiv_hist₊Sm.weights; color=(colors[2], transp), gap=0.0)
#dist₊Sm = lines!(ax_divtime, tspan, pdf.(interdiv_dist₊Sm, tspan); color=colors[2], linewidth=2, linestyle=:dash)
#axislegend(ax_divtime, [LineElement(color = (:black, 0.8), linestyle = :dash),], ["Moment matched Gamma distribution", ], framevisible = false)

birth_hist₋Sm = normalize(fit(Histogram, birth_dist₋Sm .* cf₋Sm, 1:40); mode=:pdf)
stairs!(ax_birth, collect(midpoints(birth_hist₋Sm.edges[1])), birth_hist₋Sm.weights; color=colors[1], step=:center)
barplot!(ax_birth, collect(midpoints(birth_hist₋Sm.edges[1])), birth_hist₋Sm.weights; color=(colors[1], transp), gap=0.0)
birth_hist₊Sm = normalize(fit(Histogram, birth_dist₊Sm .* cf₋Sm, 1:40); mode=:pdf)
stairs!(ax_birth, collect(midpoints(birth_hist₊Sm.edges[1])), birth_hist₊Sm.weights; color=colors[2], step=:center)
barplot!(ax_birth, collect(midpoints(birth_hist₊Sm.edges[1])), birth_hist₊Sm.weights; color=(colors[2], transp), gap=0.0)
xlims!(ax_birth, (5.0, 20.))

div_hist₋Sm = normalize(fit(Histogram, div_dist₋Sm .* cf₋Sm, 1:40); mode=:pdf)
stairs!(ax_div, collect(midpoints(div_hist₋Sm.edges[1])), div_hist₋Sm.weights; color=colors[1], step=:center)
barplot!(ax_div, collect(midpoints(div_hist₋Sm.edges[1])), div_hist₋Sm.weights; color=(colors[1], transp), gap=0.0)
div_hist₊Sm = normalize(fit(Histogram, div_dist₊Sm .* cf₋Sm, 1:40); mode=:pdf)
stairs!(ax_div, collect(midpoints(div_hist₊Sm.edges[1])), div_hist₊Sm.weights; color=colors[2], step=:center)
barplot!(ax_div, collect(midpoints(div_hist₊Sm.edges[1])), div_hist₊Sm.weights; color=(colors[2], transp), gap=0.0)
xlims!(ax_div, (10.0, 30.))

birthn_hist₋Sm = normalize(fit(Histogram, birth_dist₋Sm, 1:200); mode=:pdf)
stairs!(ax_birthn, collect(midpoints(birthn_hist₋Sm.edges[1])), birthn_hist₋Sm.weights; color=colors[1], step=:center)
barplot!(ax_birthn, collect(midpoints(birthn_hist₋Sm.edges[1])), birthn_hist₋Sm.weights; color=(colors[1], transp), gap=0.0)
birthn_hist₊Sm = normalize(fit(Histogram, birth_dist₊Sm, 1:200); mode=:pdf)
stairs!(ax_birthn, collect(midpoints(birthn_hist₊Sm.edges[1])), birthn_hist₊Sm.weights; color=colors[2], step=:center)
barplot!(ax_birthn, collect(midpoints(birthn_hist₊Sm.edges[1])), birthn_hist₊Sm.weights; color=(colors[2], transp), gap=0.0)
xlims!(ax_birthn, (25.0, 100.))

divn_hist₋Sm = normalize(fit(Histogram, div_dist₋Sm, 1:200); mode=:pdf)
stairs!(ax_divn, collect(midpoints(divn_hist₋Sm.edges[1])), divn_hist₋Sm.weights; color=colors[1], step=:center)
barplot!(ax_divn, collect(midpoints(divn_hist₋Sm.edges[1])), divn_hist₋Sm.weights; color=(colors[1], transp), gap=0.0)
divn_hist₊Sm = normalize(fit(Histogram, div_dist₊Sm, 1:200); mode=:pdf)
stairs!(ax_divn, collect(midpoints(divn_hist₊Sm.edges[1])), divn_hist₊Sm.weights; color=colors[2], step=:center)
barplot!(ax_divn, collect(midpoints(divn_hist₊Sm.edges[1])), divn_hist₊Sm.weights; color=(colors[2], transp), gap=0.0)

xlims!(ax_divn, (50.0, 200.))
#Legend(fig[5, 1:2], ax_mdfluor; orientation = :horizontal, tellwidth = false, tellheight = true, framevisible = false)
#
data_legend = [PolyElement(color=colors[1]), PolyElement(color=(colors[2]))]

Legend(fig[5,1:2], data_legend, ["No treatment", "Treatment"]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true)

CairoMakie.activate!()
mkpath("$(plotsdir())/wakamoto/")
save("$(plotsdir())/wakamoto/cf_$strain.pdf", fig)


