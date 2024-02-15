using Glob
using DataFrames
using CSV
using JLD2
using MLJ

using Catalyst
using AgentBasedFSP
using KernelDensity
using Distributions
using QuadGK
using Interpolations
import NaNMath

using Integrals
using CairoMakie
using ColorSchemes

colors = ColorSchemes.Egypt.colors

using StatProfilerHTML
include("$(srcdir())/fitting_helpers.jl")
include("utils.jl")

# Model.
# Register the parameters and variables used.
@parameters t k b
@variables Protein(t)

working_trn = 250
tslice = 5

# Load in data.
strain = "F3NW"
files₋Sm=glob("*/Results/*.xls", "$(datadir())/exp_raw/Wakamoto/ExperimentalData/$strain-Sm")
files₊Sm=glob("*/Results/*.xls", "$(datadir())/exp_raw/Wakamoto/ExperimentalData/$strain+Sm")
dfs₋Sm = DataFrame.(CSV.File.(files₋Sm))
dfs₊Sm = DataFrame.(CSV.File.(files₊Sm))

# Read in the fluorecsent protein distribution.
fbirth₋Sm = "$(datadir())/exp_pro/Wakamoto/$strain-Sm/birth.csv"
fdiv₋Sm = "$(datadir())/exp_pro/Wakamoto/$strain-Sm/div.csv"
df_birth₋Sm = DataFrame(CSV.File(fbirth₋Sm, header=false))
df_div₋Sm = DataFrame(CSV.File(fdiv₋Sm, header=false))

fbirth₊Sm = "$(datadir())/exp_pro/Wakamoto/$strain+Sm/birth.csv"
fdiv₊Sm = "$(datadir())/exp_pro/Wakamoto/$strain+Sm/div.csv"
df_birth₊Sm = DataFrame(CSV.File(fbirth₊Sm, header=false))
df_div₊Sm = DataFrame(CSV.File(fdiv₊Sm, header=false))

f₊Sm = jldopen("$(datadir())/exp_pro/Wakamoto/$strain+Sm/lineages.jld2", "r")
lineages₊Sm = f₊Sm["lineages"]
cf₊Sm = f₊Sm["cf"]
close(f₊Sm)

f₋Sm = jldopen("$(datadir())/exp_pro/Wakamoto/$strain-Sm/lineages.jld2", "r")
lineages₋Sm = f₋Sm["lineages"]
cf₋Sm = f₋Sm["cf"]
close(f₋Sm)

filter!(row -> 1 < row.Column1 < 50 && row.Column3 < 175, df_div₊Sm)
filter!(row -> 1 < row.Column1 < 50 && row.Column3 < 175, df_div₋Sm)
#filter!(row -> 1 < row.Column1 < 50 && row.Column3 < working_trn, df_div₊Sm)
#filter!(row -> 1 < row.Column1 < 50 && row.Column3 < working_trn, df_div₋Sm)

interdiv_times₋Sm = df_div₋Sm[!, :Column1]
interdiv_times₊Sm = df_div₊Sm[!, :Column1]

cycles₋Sm = rmnothing(vcat(make_cell_cycle.(lineages₋Sm)...));
fluors₋Sm = fluorescence_trajectory.(cycles₋Sm) ./ cf₋Sm
filter!(x -> 1 < length(x) < 50 && x[end] < 250, fluors₋Sm)
fluors₋Sm_train, fluors₋Sm_test = partition(fluors₋Sm, 0.7; shuffle=true)

cycles₊Sm = rmnothing(vcat(make_cell_cycle.(lineages₊Sm)...));
fluors₊Sm = fluorescence_trajectory.(cycles₊Sm) ./ cf₋Sm
filter!(x -> 1 < length(x) < 50 && x[end] < 250, fluors₊Sm)
fluors₊Sm_train, fluors₊Sm_test = partition(fluors₊Sm, 0.7; shuffle=true)

interdiv_times₊Sm = df_div₊Sm[!, :Column1]
interdiv_times₋Sm = df_div₋Sm[!, :Column1]

interdiv_dist₋Sm = interdiv_dist_kde_pop(interdiv_times₋Sm; time_slice=tslice, bounds=(0.0, 300.0), bandwidth=3)
interdiv_dist₊Sm = interdiv_dist_kde_pop(interdiv_times₊Sm; time_slice=tslice, bounds=(0.0, 300.0), bandwidth=3)

haz₋Sm(t) = interdiv_dist₋Sm.hazard(t)
haz₊Sm(t) = interdiv_dist₊Sm.hazard(t)
@register haz₋Sm(t)
@register haz₊Sm(t)

@parameters t v[1:5]
@variables Protein(t)
@named paramrn = ReactionSystem([], t, [Protein], [v...])

nknots = 3
nstd = 1.5
pbounds = fill((-10.0, 5.0), nknots)
qbounds = [(0.0, 100.0),]
x0 = fill(0.0, nknots)
q0 = (2.0, true) # optimise the ridge regression hyperparameter
priorσ = 4.0
max_evals = 100
restarts = 3
opt_params = (pbounds=pbounds, qbounds=qbounds, x0=x0, q0=q0, tslice=tslice, strain=strain, likelihood=likelihood_hazard_pop, max_evals, restarts=restarts) 

mcsampler = Turing.Gibbs(fill(Turing.MH(), 5)...)
nsamples = 1000
mcmc_params = (sampler=mcsampler, nsamples=nsamples, x0=x0, priorσ=priorσ, strain=strain, tslice=tslice, likelihood=likelihood_hazard_pop)


# RFP -Sm
quant₋Sm = quantile(last.(fluors₋Sm), [0.05, 0.95])
knots₋Sm = range(quant₋Sm[1], stop=quant₋Sm[2], length=nknots) 
itpp_spline₋Sm = Interp(x0, knots₋Sm)
Symbolics.@register_symbolic (itpp_spline₋Sm::Interp)(x, p1, p2, p3, p4, p5)
Symbolics.@register_symbolic (itpp_spline₋Sm::Interp)(x, p1, p2, p3)
Symbolics.@register_symbolic (itpp_spline₋Sm::Interp)(x, p1, p2)

bopt₋Sm = fit_spline(
    haz₋Sm(t)*itpp_spline₋Sm(Protein, [v...][1:nknots]...),
    paramrn; 
    data=fluors₋Sm, 
    name="notreatment",
    opt_params...
) 
#bopt₋Sm = skopt.load("data/sim/checkpoint_hazard_$(strain)_notreatment_5.pkl")

# RFP +Sm
quant₊Sm = quantile(last.(fluors₊Sm), [0.05, 0.95])
knots₊Sm = range(quant₊Sm[1], stop=quant₊Sm[2], length=nknots) 
itpp_spline₊Sm = Interp(x0, knots₊Sm)
Symbolics.@register_symbolic (itpp_spline₊Sm::Interp)(x, p1, p2, p3, p4, p5)
Symbolics.@register_symbolic (itpp_spline₊Sm::Interp)(x, p1, p2, p3)
Symbolics.@register_symbolic (itpp_spline₊Sm::Interp)(x, p1, p2)

bopt₊Sm = fit_spline(
    haz₊Sm(t)*itpp_spline₊Sm(Protein, [v...][1:nknots]...),
    paramrn; 
    data=fluors₊Sm, 
    name="treatment",
    opt_params...
   ) 
#bopt₊Sm = skopt.load("data/sim/checkpoint_hazard_$(strain)_treatment_5.pkl")
#
#bopt₊Sm₋Sm = fit_spline(
#    haz₋Sm(t)*itpp_spline₊Sm(Protein, [v...][1:nknots]...),
#    paramrn; 
#    data=fluors₊Sm, 
#    name="treatmentwt",
#    opt_params...) 
#bopt₊Sm₋Sm = skopt.load("data/sim/checkpoint_hazard_$(strain)_treatmentwt_5.pkl")

#opt₋Sm = fit_spline_mcmc(
#    haz₋Sm(t)*itpp_spline₋Sm(Protein, [v...][1:nknots]...),
#    paramrn; 
#    data=fluors₋Sm, 
#    name="notreatment",
#    mcmc_params...
#    ) 
opt₋Sm = deserialize("data/sim/chainfile_$(strain)_notreatment_$(nsamples).jls")

#opt₊Sm = fit_spline_mcmc(
#    haz₊Sm(t)*itpp_spline₊Sm(Protein, [v...][1:nknots]...),
#    paramrn; 
#    data=fluors₊Sm, 
#    name="treatmentalt",
#    mcmc_params...
#   ) 
opt₊Sm = deserialize("data/sim/chainfile_$(strain)_treatmentalt_$(nsamples).jls")

#opt₊Sm_alt = fit_spline_mcmc(
#    haz₋Sm(t)*itpp_spline₊Sm(Protein, [v...][1:nknots]...),
#    paramrn; 
#    data=fluors₊Sm, 
#    name="treatment",
#    mcmc_params...
#    ) 
opt₊Sm_alt = deserialize("data/sim/chainfile_$(strain)_treatment_$(nsamples).jls")

params₋Sm = pyconvert(Vector, bopt₋Sm.x)[1:nknots]
params₊Sm = pyconvert(Vector, bopt₊Sm.x)[1:nknots]
params₊Sm₋Sm = pyconvert(Vector, bopt₊Sm₋Sm.x)[1:nknots]

itp₋Sm(x) = itpp_spline₋Sm(x, params₋Sm...)
itp₊Sm(x) = itpp_spline₊Sm(x, params₊Sm...)
itp₊Sm₋Sm(x) = itpp_spline₊Sm(x, params₊Sm₋Sm...)

@register itp₋Sm(x)
@register itp₊Sm(x)
@register itp₊Sm₋Sm(x)

params_hillfs = [1.3771485956602671, 78.24118123524796, 0.1039282608597518, 10]
@parameters L K δ n 
hillfs = L*((1-δ)*Protein^n / (K^n + Protein^n) + δ)
@named hillfsparamrn = ReactionSystem([], t, [Protein], [L, K, δ, n])
sel_hillfs = AgentBasedFSP.gen_division_rate_function(hillfs, hillfsparamrn)

fig = Figure()
#ax_marginal = Axis(fig[1,1]; title="Marginal interdivision time", xlabel="Time", ylabel="Probability density")
ax_sel = Axis(fig[1,1]; xlabel="Protein numbers", ylabel="Selection strength")
#ax_marginalx = Axis(fig[3,1]; xlabel="Protein numbers", ylabel="Probability density")

xs = 50.0:1.0:170.0

#mγ₋Sm = AgentBasedFSP.gen_division_rate_function(haz₋Sm(t), paramrn)
mγ₋Sm = AgentBasedFSP.gen_division_rate_function(haz₋Sm(t) * itp₋Sm(Protein), paramrn)
marginalγ₋Sm = empirical_marginal_cycles_γ(mγ₋Sm, params₋Sm, fluors₋Sm; tslice=5)
integ₋Sm(t) = quadgk(marginalγ₋Sm, 0.0, t)[1]
λ₋Sm = growth_rate(haz₋Sm(t)* itp₋Sm(Protein), params₋Sm, paramrn; data=fluors₋Sm, tslice=5, offset=0.0)
#λ₋Sm = growth_rate(haz₋Sm(t), params₋Sm, paramrn; data=fluors₋Sm, tslice=5)
pd₋Sm(t) = 2*exp(-λ₋Sm*t)*marginalγ₋Sm(t) * exp(-integ₋Sm(t))
#pd₋Sm(t) = marginalγ₋Sm(t) * exp(-integ₋Sm(t))

mγ₊Sm = AgentBasedFSP.gen_division_rate_function(haz₊Sm(t) * itp₊Sm(Protein), paramrn)
marginalγ₊Sm = empirical_marginal_cycles_γ(mγ₊Sm, params₊Sm, fluors₊Sm; tslice=5)
integ₊Sm(t) = quadgk(marginalγ₊Sm, 0.0, t)[1]
λ₊Sm = growth_rate(haz₊Sm(t) * itp₊Sm(Protein), params₊Sm, paramrn; data=fluors₊Sm, tslice=5, offset=0.0)
pd₊Sm(t) = 2*exp(-λ₊Sm*t)*marginalγ₊Sm(t) * exp(-integ₊Sm(t))

mγ₊Sm₋Sm = AgentBasedFSP.gen_division_rate_function(haz₋Sm(t) * itp₊Sm₋Sm(Protein), paramrn)
marginalγ₊Sm₋Sm = empirical_marginal_cycles_γ(mγ₊Sm₋Sm, params₊Sm₋Sm, fluors₊Sm; tslice=5)
integ₊Sm₋Sm(t) = quadgk(marginalγ₊Sm₋Sm, 0.0, t)[1]
λ₊Sm₋Sm = growth_rate(haz₋Sm(t) * itp₊Sm₋Sm(Protein), params₊Sm₋Sm, paramrn; data=fluors₊Sm, tslice=5, offset=0.0)
pd₊Sm₋Sm(t) = 2*exp(-λ₊Sm₋Sm*t)*marginalγ₊Sm₋Sm(t) * exp(-integ₊Sm₋Sm(t))

integ(t) = quadgk(haz₋Sm, 0.0, t)[1]
pdfhaz(t) = 2*exp(-interdiv_dist₋Sm.λ*t)*haz₋Sm(t)*exp(-integ(t))

ts = 0.0:1.0:135.0
#hist!(ax_marginal, length.(fluors₋Sm) .* tslice .- 0.5*tslice; bins=0.0:5.0:150.0, normalization=:pdf, color=(colors[1], 0.5))
#hist!(ax_marginal, length.(fluors₋Sm) .* tslice .- 0.5*tslice; bins=2.5:5.0:150.0, normalization=:pdf, color=(colors[1], 0.5))
#hist!(ax_marginal, length.(fluors₊Sm) .* tslice .- 0.5*tslice; bins=2.5:5.0:150.0, normalization=:pdf, color=(colors[2], 0.5))
#lines!(ax_marginal, ts, pdf.(Ref(interdiv_dist₋Sm), ts); color=colors[1])
#lines!(ax_marginal, ts, pdfhaz.(ts); color=colors[1])
#lines!(ax_marginal, ts, haz₊Sm.(ts); color=colors[2])
#lines!(ax_marginal, ts, pd₋Sm.(ts); color=colors[1])
#lines!(ax_marginal, ts, pd₊Sm.(ts); color=colors[2])
#lines!(ax_marginal, ts, pd₊Sm₋Sm.(ts); color=colors[2], linestyle=:dash)

#lines!(ax_sel, ts, marginalγ₋Sm.(ts); color=colors[1])
#lines!(ax_sel, ts, interdiv_dist₋Sm.hazard.(ts); color=colors[2])
lines!(ax_sel, xs, itp₋Sm.(xs); color=colors[1], label="Fitted (no treatment)")
scatter!(ax_sel, knots₋Sm, itp₋Sm.(knots₋Sm); color=colors[1])
lines!(ax_sel, xs, itp₊Sm.(xs); color=colors[2], label="Fitted (treatment)")
scatter!(ax_sel, knots₊Sm, itp₊Sm.(knots₊Sm); color=colors[2])

Legend(fig[2,1], ax_sel, framevisible = false, tellwidth=false, tellheight=true, orientation=:horizontal)

#lines!(ax_sel, xs, itp₊Sm₋Sm.(xs); color=colors[3])
#lines!(ax_sel, xs, sel_hillfs.(xs, Ref(params_hillfs), Ref(0.0)); color=:gray)

#scatter!(ax_sel, knots₋Sm, itp₋Sm.(knots₋Sm); color=colors[1], linestyle=:dash)
#scatter!(ax_sel, knots₊Sm, itp₊Sm.(knots₊Sm); color=colors[2], linestyle=:dash)
#scatter!(ax_sel, knots₊Sm, itp₊Sm₋Sm.(knots₊Sm); color=colors[2], linestyle=:dash)

#sum₋Sm = summarystats(opt₋Sm)
#mcmc_itp₋Sm = make_itpp_spline(knots₋Sm, mean(opt₋Sm).nt[:mean])
#mcmc_mean₋Sm(x) = exp(mcmc_itp₋Sm(x))
#mcmc_meanvar₋Sm = map(x -> mean_and_var(exp.(x)), [eachcol(opt₋Sm.value.data[:, 1:nknots, 1])...])
#lines!(ax_sel, xs, mcmc_mean₋Sm.(xs); color=colors[1], linestyle=:dash)
#rangebars!(ax_sel, knots₋Sm, 
#    mcmc_mean₋Sm.(knots₋Sm) .- 1.96*sqrt.(getindex.(mcmc_meanvar₋Sm, 2)), 
#    mcmc_mean₋Sm.(knots₋Sm) .+ 1.96*sqrt.(getindex.(mcmc_meanvar₋Sm, 2)); 
#    color=colors[1], whiskerwidth = 10)

#sum₊Sm_alt = summarystats(opt₊Sm_alt)
#mcmc_itp₊Sm_alt = make_itpp_spline(knots₊Sm, mean(opt₊Sm_alt).nt[:mean])
#mcmc_mean₊Sm_alt(x) = exp(mcmc_itp₊Sm_alt(x))
#mcmc_meanvar₊Sm_alt = map(x -> mean_and_var(exp.(x)), [eachcol(opt₊Sm_alt.value.data[:, 1:nknots, 1])...])
#lines!(ax_sel, xs, mcmc_mean₊Sm_alt.(xs); color=colors[3], linestyle=:dash)
#rangebars!(ax_sel, knots₊Sm, 
#    mcmc_mean₊Sm_alt.(knots₊Sm) .- 1.96*sqrt.(getindex.(mcmc_meanvar₊Sm_alt, 2)), 
#    mcmc_mean₊Sm_alt.(knots₊Sm) .+ 1.96*sqrt.(getindex.(mcmc_meanvar₊Sm_alt, 2)); 
#    color=colors[3], whiskerwidth = 10)
#
#sum₊Sm = summarystats(opt₊Sm)
#mcmc_itp₊Sm = make_itpp_spline(knots₊Sm, mean(opt₊Sm).nt[:mean])
#mcmc_mean₊Sm(x) = exp(mcmc_itp₊Sm(x))
#mcmc_meanvar₊Sm = map(x -> mean_and_var(exp.(x)), [eachcol(opt₊Sm.value.data[:, 1:nknots, 1])...])
#lines!(ax_sel, xs, mcmc_mean₊Sm.(xs); color=colors[2], linestyle=:dash)
#rangebars!(ax_sel, knots₊Sm, 
#    mcmc_mean₊Sm.(knots₊Sm) .- 1.96*sqrt.(getindex.(mcmc_meanvar₊Sm, 2)), 
#    mcmc_mean₊Sm.(knots₊Sm) .+ 1.96*sqrt.(getindex.(mcmc_meanvar₊Sm, 2)); 
#    color=colors[2], whiskerwidth = 10)

#hist!(ax_marginalx, last.(fluors₋Sm); bins=0.0:5.0:250.0, normalization=:pdf, color=(colors[1], 0.5))
#hist!(ax_marginalx, last.(fluors₊Sm); bins=0.0:5.0:250.0, normalization=:pdf, color=(colors[2], 0.5))
#xlims!(ax_marginalx, (50.0, 170.0))
xlims!(ax_sel, (50.0, 170.0))

mkpath("$(plotsdir())/wakamoto/")
save("$(plotsdir())/wakamoto/selfuns.pdf", fig)


## Pick a slice, 80.
#slice = 80
## Filter all lineages up to tslice=80.
#ages = [] 
#
#for lin in lineages₋Sm
#    filter!.(lin -> lin.Slice <= slice, lin.lineages)
#    filter!(lin -> maximum(x -> x.Slice, eachrow(lin)) >= slice, lin.lineages)
#    mother_idxs = lin.mother_cells[:, :1]
#    df_ = unique(subset.(lin.lineages, :1 => x -> in.(x, Ref(mother_idxs))))
#    ag_ = [slice - maximum(df[!, :Slice]) for df in df_]
#    push!(ages, ag_...)
#end
#
#agesx = [count(x -> length(x) >= i, fluors₋Sm) for i in 1:1:30]
#agesh = Histogram(collect(0:1:30) .* tslice, agesx, :right, false)
#agesh = normalize(agesh, mode=:pdf)
##agesx = agesx ./ sum(agesx)
#
#fig = Figure()
#ts = 0.0:1.0:200.0
#ax1 = Axis(fig[1,1])
#ax2 = Axis(fig[2,1])
#hist!(ax1, ages .* tslice; bins=0.0:5.0:200.0, normalization=:pdf)
#lines!(ax1, ts, interdiv_dist₋Sm.age.(ts))
#stairs!(ax1, collect(midpoints(agesh.edges[1])), agesh.weights) 
#
#hist!(ax2, interdiv_times₋Sm .* tslice .- 0.5 * tslice; bins=0.0:5.0:200.0, normalization=:pdf)
#hist!(ax2, length.(fluors₋Sm) .* tslice .- 0.5 * tslice; bins=0.0:5.0:200.0, normalization=:pdf)
#lines!(ax2, ts, pdf(interdiv_dist₋Sm.kde, ts))
#save("$(plotsdir())/wakamoto/ages.pdf", fig)

#mγ₋Sm_ = AgentBasedFSP.gen_division_rate_function(haz₋Sm(t), paramrn)
#likelihood_hazard_pop(mγ₋Sm, params₋Sm, q0[1]; data=fluors₋Sm_test, inf=1e10, tslice=5, x0=x0) / length(fluors₋Sm_test)
#likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm_test, inf=1e10, tslice=5, x0=x0) / length(fluors₋Sm_test)
#likelihood_hazard_pop(mγ₋Sm, params₋Sm, q0[1]; data=fluors₋Sm_train, inf=1e10, tslice=5, x0=x0) / length(fluors₋Sm_train)
#likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm_train, inf=1e10, tslice=5, x0=x0) / length(fluors₋Sm_train)

#mγ₋Sm_ = AgentBasedFSP.gen_division_rate_function(haz₋Sm(t), paramrn)
#println("Spline")
#5*log(length(fluors₋Sm_test)) - 2*likelihood_hazard_pop(mγ₋Sm, params₋Sm, q0[1]; data=fluors₋Sm_test, inf=1e10, tslice=5, x0=x0, offset=0.5)
#5*log(length(fluors₋Sm_train)) - 2*likelihood_hazard_pop(mγ₋Sm, params₋Sm, q0[1]; data=fluors₋Sm_train, inf=1e10, tslice=5, x0=x0, offset=0.5)
#println("No spline")
#-2*likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm_test, inf=1e10, tslice=5, x0=x0, offset=0.5)
#-2*likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm_train, inf=1e10, tslice=5, x0=x0, offset=0.5)

#mγ₊Sm_ = AgentBasedFSP.gen_division_rate_function(haz₊Sm(t), paramrn)
#likelihood_hazard_pop(mγ₊Sm, params₊Sm, q0[1]; data=fluors₊Sm_test, inf=1e10, tslice=5, x0=x0) / length(fluors₊Sm_test)
#likelihood_hazard_pop(mγ₊Sm_, params₊Sm, q0[1]; data=fluors₊Sm_test, inf=1e10, tslice=5, x0=x0) / length(fluors₊Sm_test)
#likelihood_hazard_pop(mγ₊Sm, params₊Sm, q0[1]; data=fluors₊Sm_train, inf=1e10, tslice=5, x0=x0) / length(fluors₊Sm_train)
#likelihood_hazard_pop(mγ₊Sm_, params₊Sm, q0[1]; data=fluors₊Sm_train, inf=1e10, tslice=5, x0=x0) / length(fluors₊Sm_train)

#mγ₊Sm_ = AgentBasedFSP.gen_division_rate_function(haz₊Sm(t), paramrn)
#5*log(length(fluors₊sm_test)) - 2*likelihood_hazard_pop(mγ₊sm, params₊sm, q0[1]; data=fluors₊sm_test, inf=1e10, tslice=5, x0=x0 )
#5*log(length(fluors₊sm_train)) - 2*likelihood_hazard_pop(mγ₊sm, params₊sm, q0[1]; data=fluors₊sm_train, inf=1e10, tslice=5, x0=x0)
#-2*likelihood_hazard_pop(mγ₊sm_, params₊sm, q0[1]; data=fluors₊sm_test, inf=1e10, tslice=5, x0=x0)
#-2*likelihood_hazard_pop(mγ₊sm_, params₊sm, q0[1]; data=fluors₊sm_train, inf=1e10, tslice=5, x0=x0)
#

mγ₋Sm_ = AgentBasedFSP.gen_division_rate_function(haz₋Sm(t), paramrn)
3*log(length(fluors₋Sm)) - 2*likelihood_hazard_pop(mγ₋Sm, params₋Sm, q0[1]; data=fluors₋Sm, inf=1e10, tslice=5, x0=x0)
-2*likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm, inf=1e10, tslice=5, x0=x0)


mγ₊Sm_ = AgentBasedFSP.gen_division_rate_function(haz₊Sm(t), paramrn)
3*log(length(fluors₊Sm)) - 2*likelihood_hazard_pop(mγ₊Sm, params₊Sm, q0[1]; data=fluors₊Sm, inf=1e10, tslice=5, x0=x0)
-2*likelihood_hazard_pop(mγ₊Sm_, params₊Sm, q0[1]; data=fluors₊Sm, inf=1e10, tslice=5, x0=x0)



#julia> 5*log(length(fluors₊Sm)) - 2*likelihood_hazard_pop(mγ₊Sm, params₊Sm, q0[1]; data=fluors₊Sm, inf=1e10, tslice=5, x0=x0)
#63585.22015755291
#
#julia> -2*likelihood_hazard_pop(mγ₊Sm_, params₊Sm, q0[1]; data=fluors₊Sm, inf=1e10, tslice=5, x0=x0)
#63776.81583924206
#
#julia> 5*log(length(fluors₋Sm)) - 2*likelihood_hazard_pop(mγ₋Sm, params₋Sm, q0[1]; data=fluors₋Sm, inf=1e10, tslice=5, x0=x0)
#72051.70033833991
#
#julia> -2*likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm, inf=1e10, tslice=5, x0=x0)
#72328.01796952363
