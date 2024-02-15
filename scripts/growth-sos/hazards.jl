using Glob
using DataFrames
using CSV
using JLD2

using Catalyst
using AgentBasedFSP
using KernelDensity
using Distributions
using QuadGK
using Interpolations
import NaNMath
using Turing

using Integrals
using CairoMakie
using ColorSchemes

colors = ColorSchemes.Egypt.colors

using StatProfilerHTML
include("$(srcdir())/fitting_helpers.jl")
include("utils.jl")

function interdiv_dist_kde(interdiv_times; time_slice, bounds, bandwidth)
    x_ = interdiv_times .* time_slice 
    return CumulativeKDE(InterpKDE(kde(x_; bandwidth=bandwidth, boundary=bounds)), tend=bounds[end])
end

strain = "sos"
working_trn = 200
tslice = 8

# Data.
# Read in the fluorecsent protein distribution.
fbirth_wt = "$(datadir())/exp_pro/growth_sos/wt/birth.csv"
fdiv_wt = "$(datadir())/exp_pro/growth_sos/wt/div.csv"
df_birth_wt = DataFrame(CSV.File(fbirth_wt, header=false))
df_div_wt = DataFrame(CSV.File(fdiv_wt, header=false))

fbirth_2pal = "$(datadir())/exp_pro/growth_sos/2pal/birth.csv"
fdiv_2pal = "$(datadir())/exp_pro/growth_sos/2pal/div.csv"
df_birth_2pal = DataFrame(CSV.File(fbirth_2pal, header=false))
df_div_2pal = DataFrame(CSV.File(fdiv_2pal, header=false))

f = jldopen("$(datadir())/exp_pro/growth_sos/cycles.jld2", "r")
cycles_2pal = f["cycles_2pal"]
cycles_wt = f["cycles_wt"]
cf_wt = f["cf_wt"]
cf_2pal = f["cf_2pal"]

fluors_2pal = cycles_2pal ./ cf_2pal
fluors_wt = cycles_wt ./ cf_wt

filter!(x -> last(x) < working_trn && length(x) < 40, fluors_2pal)
filter!(x -> last(x) < working_trn && length(x) < 40, fluors_wt)
#filter!(x -> last(x) < working_trn, fluors_2pal)
#filter!(x -> last(x) < working_trn, fluors_wt)

filter!(row -> 0 < row.Column1 < 40 && row.Column3 < working_trn, df_div_wt)
filter!(row -> 0 < row.Column1 < 40 && row.Column3 < working_trn, df_div_2pal)

interdiv_times_2pal = df_div_2pal[!, :Column1]
interdiv_times_wt = df_div_wt[!, :Column1]

interdiv_dist_2pal = interdiv_dist_kde_mother(interdiv_times_2pal; time_slice=tslice, bounds=(0.0, 500.0), bandwidth=8)
interdiv_dist_wt = interdiv_dist_kde_mother(interdiv_times_wt; time_slice=tslice, bounds=(0.0, 500.0), bandwidth=8)

@parameters t v[1:8]
@variables Protein(t) 

haz_wt(t) = interdiv_dist_wt.hazard(t)
haz_2pal(t) = interdiv_dist_2pal.hazard(t)
@register haz_wt(t)
@register haz_2pal(t)

nknots = 5
nstd = 1.5
pbounds = fill((-10.0, 5.0), nknots)
qbounds = [(0.0, 100.0),]
x0 = fill(0.0, nknots)
q0 = (2.0, false) # optimise the ridge regression hyperparameter
σ0 = 0.5 
priorσ = 0.5
max_evals = 100
restarts = 5
opt_params = (pbounds=pbounds, qbounds=qbounds, x0=x0, q0=q0, tslice=tslice, strain=strain, likelihood=likelihood_hazard_mother, max_evals, restarts=restarts) 

mcsampler = Turing.Gibbs(fill(Turing.MH(), 5)...)
nsamples = 20000
mcmc_params = (sampler=mcsampler, nsamples=nsamples, x0=x0, priorσ=priorσ, strain=strain, tslice=tslice, likelihood=likelihood_hazard_mother)

@named paramrn = ReactionSystem([], t, [Protein], [v[1:nknots]...])

# wt
quant_wt = quantile(last.(fluors_wt), [0.05, 0.95])
knots_wt = range(quant_wt[1], stop=quant_wt[2], length=nknots) 
itpp_spline_wt = Interp(x0, knots_wt)
Symbolics.@register_symbolic (itpp_spline_wt::Interp)(x, p1, p2, p3, p4, p5)

#bopt_wt = fit_spline(
#    haz_wt(t)*itpp_spline_wt(Protein, [v...][1:nknots]...),
#    paramrn;
#    data=fluors_wt, 
#    name="wt",
#    opt_params...
#   ) 
#bopt_wt = skopt.load("data/sim/checkpoint_hazard_$(strain)_wt_5.pkl")

const_wt(x, v) = exp(v)
@register const_wt(x, v)
@named paramrn1 = ReactionSystem([], t, [Protein], [v[1],])

# 2pal 
quant_2pal = quantile(last.(fluors_2pal), [0.05, 0.95])
knots_2pal = range(quant_2pal[1], stop=quant_2pal[2], length=nknots) 
itpp_spline_2pal = Interp(x0, knots_2pal)
Symbolics.@register_symbolic (itpp_spline_2pal::Interp)(x, p1, p2, p3, p4, p5)

#bopt_2palwt = fit_spline(
#    haz_wt(t)*itpp_spline_2pal(Protein, [v...][1:nknots]...),
#    paramrn;
#    data=fluors_2pal, 
#    name="2pal",
#    opt_params...
#   )
#bopt_2pal = skopt.load("data/sim/checkpoint_hazard_$(strain)_2pal_5.pkl")

#bopt_2pal = fit_spline(
#    haz_2pal(t)*itpp_spline_2pal(Protein, [v...][1:nknots]...),
#    paramrn;
#    data=fluors_2pal, 
#    name="2pal",
#    opt_params...
#   )
#bopt_2pal = skopt.load("data/sim/checkpoint_hazard_$(strain)_2pal_5.pkl")

#opt_wt_mcmc = fit_spline_mcmc(
#    haz_wt(t)*itpp_spline_wt(Protein, [v...][1:nknots]...),
#    paramrn; 
#    data=fluors_wt, 
#    name="wt",
#    mcmc_params...
#   )
#opt_wt_mcmc = deserialize("data/sim/chainfile_$(strain)_wt_$(nsamples).jls")

#opt_2palwt_mcmc = fit_spline_mcmc(
#    haz_wt(t)*itpp_spline_2pal(Protein, [v...][1:nknots]...),
#    paramrn; 
#    data=fluors_2pal, 
#    name="2palalt",
#    mcmc_params...
#   )
#opt_2palwt_mcmc = deserialize("data/sim/chainfile_$(strain)_2palalt_$(nsamples).jls")

#opt_2pal_mcmc = fit_spline_mcmc(
#    haz_2pal(t)*itpp_spline_2pal(Protein, [v...][1:nknots]...),
#    paramrn; 
#    data=fluors_2pal, 
#    name="2pal",
#    mcmc_params...
#   )
#opt_2pal_mcmc = deserialize("data/sim/chainfile_$(strain)_2pal_$(nsamples).jls")

params_wt = pyconvert(Vector, bopt_wt.x)
params_2pal = pyconvert(Vector, bopt_2pal.x)
params_2palwt = pyconvert(Vector, bopt_2palwt.x)

itp_wt(x) = itpp_spline_wt(x, params_wt...)
itp_2pal(x) = itpp_spline_2pal(x, params_2pal...)
itp_2palwt(x) = itpp_spline_2pal(x, params_2palwt...)

@register itp_wt(x)
@register itp_2pal(x)
@register itp_2palwt(x)

fig = Figure()
ax_marginal = Axis(fig[1,1])
ax_sel = Axis(fig[2,1]; xlabel="Protein counts", ylabel="Selection strength")
ax_marginalx = Axis(fig[3,1])
xs = 0.0:1.0:120.0

transp = 0.4
stairstransp = 0.4

mγ_wt = AgentBasedFSP.gen_division_rate_function(haz_wt(t)*itp_wt(Protein), paramrn)
marginalγ_wt = empirical_marginal_cycles_γ(mγ_wt, params_wt, fluors_wt; tslice=8)

mγ_2pal = AgentBasedFSP.gen_division_rate_function(haz_2pal(t)*itp_2pal(Protein), paramrn)
marginalγ_2pal = empirical_marginal_cycles_γ(mγ_2pal, params_2pal, fluors_2pal; tslice=8)

mγ_2palhaz = AgentBasedFSP.gen_division_rate_function(haz_wt(t)*itp_2palwt(Protein), paramrn)
marginalγ_2palhaz = empirical_marginal_cycles_γ(mγ_2palhaz, params_2palwt, fluors_2pal; tslice=8)

integ_wt(t) = quadgk(marginalγ_wt, 0.0, t)[1]
pd_wt(t) = marginalγ_wt(t) * exp(-integ_wt(t))

integ_2pal(t) = quadgk(marginalγ_2pal, 0.0, t)[1]
pd_2pal(t) = marginalγ_2pal(t) * exp(-integ_2pal(t))

integ_2palhaz(t) = quadgk(marginalγ_2palhaz, 0.0, t)[1]
pd_2palhaz(t) = marginalγ_2palhaz(t) * exp(-integ_2palhaz(t))

ts = 0.0:1.0:200.0
hist!(ax_marginal, length.(fluors_wt) .* tslice .- 0.5*tslice; bins=0.0:8.0:250.0, normalization=:pdf, color=(colors[1], 0.5))
hist!(ax_marginal, length.(fluors_2pal) .* tslice .- 0.5*tslice; bins=0.0:8.0:250.0, normalization=:pdf, color=(colors[2], 0.5))
lines!(ax_marginal, ts, pd_wt.(ts); color=colors[1])
lines!(ax_marginal, ts, pd_2pal.(ts); color=colors[2])
#lines!(ax_marginal, ts, pd_2palhaz.(ts); color=colors[2], linestyle=:dash)

hist!(ax_marginalx, last.(fluors_wt); bins=20.0:2.5:250.0, normalization=:pdf, color=(colors[1], 0.5))
hist!(ax_marginalx, last.(fluors_2pal); bins=20.0:2.5:250.0, normalization=:pdf, color=(colors[2], 0.5))
xlims!(ax_marginalx, (0.0, 120.0))

sum_wt = summarystats(opt_wt_mcmc)
mcmc_itp_wt = make_itpp_spline(knots_wt, mean(opt_wt_mcmc).nt[:mean])
mcmc_mean_wt(x) = exp(mcmc_itp_wt(x))
mcmc_meanvar_wt = map(x -> mean_and_var(exp.(x)), [eachcol(opt_wt_mcmc.value.data[:, 1:5, 1])...])

sum_2pal = summarystats(opt_2pal_mcmc)
mcmc_itp_2pal = make_itpp_spline(knots_2pal, mean(opt_2pal_mcmc).nt[:mean])
mcmc_mean_2pal(x) = exp(mcmc_itp_2pal(x))
mcmc_meanvar_2pal = map(x -> mean_and_var(exp.(x)), [eachcol(opt_2pal_mcmc.value.data[:, 1:5, 1])...])

sum_2palwt = summarystats(opt_2palwt_mcmc)
mcmc_itp_2palwt = make_itpp_spline(knots_2pal, mean(opt_2palwt_mcmc).nt[:mean])
mcmc_mean_2palwt(x) = exp(mcmc_itp_2palwt(x))
mcmc_meanvar_2palwt = map(x -> mean_and_var(exp.(x)), [eachcol(opt_2palwt_mcmc.value.data[:, 1:5, 1])...])

lines!(ax_sel, xs, x -> 1; color=(:grey, 0.8))
lines!(ax_sel, xs, itp_wt.(xs); color=(colors[1], transp), linestyle=:dash)
lines!(ax_sel, xs, itp_2pal.(xs); color=(colors[2], transp), linestyle=:dash)
#lines!(ax_sel, xs, itp_2palwt.(xs); color=colors[2], linestyle=:dash)
#scatter!(ax_sel, knots_wt, itp_wt.(knots_wt); color=colors[1], linestyle=:dash)
#scatter!(ax_sel, knots_2pal, itp_2pal.(knots_2pal); color=colors[2], linestyle=:dash)

lines!(ax_sel, xs, mcmc_mean_wt.(xs); color=colors[1])
rangebars!(ax_sel, knots_wt, 
    mcmc_mean_wt.(knots_wt) .- 1.96*sqrt.(getindex.(mcmc_meanvar_wt, 2)), 
    mcmc_mean_wt.(knots_wt) .+ 1.96*sqrt.(getindex.(mcmc_meanvar_wt, 2)); 
    color=colors[1], whiskerwidth = 10)

lines!(ax_sel, xs, mcmc_mean_2pal.(xs); color=colors[2])
rangebars!(ax_sel, knots_2pal, 
    mcmc_mean_2pal.(knots_2pal) .- 1.96*sqrt.(getindex.(mcmc_meanvar_2pal, 2)), 
    mcmc_mean_2pal.(knots_2pal) .+ 1.96*sqrt.(getindex.(mcmc_meanvar_2pal, 2)); 
    color=colors[2], whiskerwidth = 10)

#lines!(ax_sel, xs, mcmc_mean_2palwt.(xs); color=colors[2])
#band!(ax_sel, knots_2pal, 
#    mcmc_mean_2palwt.(knots_2pal) .- 1.96*sqrt.(getindex.(mcmc_meanvar_2palwt, 2)), 
#    mcmc_mean_2palwt.(knots_2pal) .+ 1.96*sqrt.(getindex.(mcmc_meanvar_2palwt, 2)); 
#    color=colors[2], whiskerwidth = 10)

xlims!(ax_sel, (00.0, 120.0))

mkpath("$(plotsdir())/growth-sos/")
save("$(plotsdir())/growth-sos/selfuns.pdf", fig)


mγ_wt = AgentBasedFSP.gen_division_rate_function(haz_wt(t), paramrn)
mγ_wtitp = AgentBasedFSP.gen_division_rate_function(haz_wt(t)*itp_wt(t), paramrn)
#likelihood_hazard_pop(mγ_wt, params_wt, q0[1]; data=fluors₋Sm_test, inf=1e10, tslice=5, x0=x0) / length(fluors₋Sm_test)
#likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm_test, inf=1e10, tslice=5, x0=x0) / length(fluors₋Sm_test)
#likelihood_hazard_pop(mγ₋Sm, params₋Sm, q0[1]; data=fluors₋Sm_train, inf=1e10, tslice=5, x0=x0) / length(fluors₋Sm_train)
#likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm_train, inf=1e10, tslice=5, x0=x0) / length(fluors₋Sm_train)

println("Spline")
5*log(length(fluors_wt)) - 2*likelihood_hazard_pop(mγ_wtitp, params_wt, q0[1]; data=fluors_wt, inf=1e10, tslice=5, x0=x0)
println("No spline")
-2*likelihood_hazard_pop(mγ_wt, params_wt, q0[1]; data=fluors_wt, inf=1e10, tslice=5, x0=x0)
#-2*likelihood_hazard_pop(mγ₋Sm_, params₋Sm, q0[1]; data=fluors₋Sm_train, inf=1e10, tslice=5, x0=x0, offset=0.5)


mγ_2pal = AgentBasedFSP.gen_division_rate_function(haz_2pal(t), paramrn)
mγ_2palitp = AgentBasedFSP.gen_division_rate_function(haz_2pal(t)*itp_2pal(t), paramrn)
println("Spline")
5*log(length(fluors_2pal)) - 2*likelihood_hazard_pop(mγ_2palitp, params_2pal, q0[1]; data=fluors_2pal, inf=1e10, tslice=5, x0=x0)
println("No spline")
-2*likelihood_hazard_pop(mγ_2pal, params_2pal, q0[1]; data=fluors_2pal, inf=1e10, tslice=5, x0=x0)


#mγ₊Sm_ = AgentBasedFSP.gen_division_rate_function(haz₊Sm(t), paramrn)
#likelihood_hazard_pop(mγ₊Sm, params₊Sm, q0[1]; data=fluors₊Sm_test, inf=1e10, tslice=5, x0=x0) / length(fluors₊Sm_test)
#likelihood_hazard_pop(mγ₊Sm_, params₊Sm, q0[1]; data=fluors₊Sm_test, inf=1e10, tslice=5, x0=x0) / length(fluors₊Sm_test)
#likelihood_hazard_pop(mγ₊Sm, params₊Sm, q0[1]; data=fluors₊Sm_train, inf=1e10, tslice=5, x0=x0) / length(fluors₊Sm_train)
#likelihood_hazard_pop(mγ₊Sm_, params₊Sm, q0[1]; data=fluors₊Sm_train, inf=1e10, tslice=5, x0=x0) / length(fluors₊Sm_train)

#5*log(length(fluors₊Sm_test)) - 2*likelihood_hazard_pop(mγ₊Sm, params₊Sm, q0[1]; data=fluors₊Sm_test, inf=1e10, tslice=5, x0=x0 )
#5*log(length(fluors₊Sm_train)) - 2*likelihood_hazard_pop(mγ₊Sm, params₊Sm, q0[1]; data=fluors₊Sm_train, inf=1e10, tslice=5, x0=x0)
#-2*likelihood_hazard_pop(mγ₊Sm_, params₊Sm, q0[1]; data=fluors₊Sm_test, inf=1e10, tslice=5, x0=x0)
#-2*likelihood_hazard_pop(mγ₊Sm_, params₊Sm, q0[1]; data=fluors₊Sm_train, inf=1e10, tslice=5, x0=x0)

