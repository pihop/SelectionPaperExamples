using Glob
using DataFrames
using CSV
using JLD2

using Catalyst
using AgentBasedFSP
using KernelDensity
using Distributions
using QuadGK
using DataInterpolations
using MLJ

using StatProfilerHTML

include("$(srcdir())/fitting_helpers.jl")
include("utils.jl")

# Model.
# Register the parameters and variables used.
@parameters t k r L K δ n b δr Kr m d
@species Protein(t)

hillfs = L*((1-δ)*Protein^n / (K^n + Protein^n) + δ)
hillfr = (1-δr) * Kr^m / (Kr^m + Protein^m) + δr 

@parameters v[1:5]
#@variables Protein(t)

# Bursty model. Note that maximum burst size approximation is made. Better
# compute for the upper bound.
rxs_burst = [Reaction(k*(b-1)^(n-1)/b^(n+1), nothing, [Protein], nothing, [n]) for n in 1:20]
# We are fitting (k'b, b)
rxs_deg = [Reaction(d, [Protein], nothing)]
@named bursty_rn = ReactionSystem([rxs_burst...,], t, [Protein], [k, b, L, K, δ, n])
@named bursty_rn_base = ReactionSystem([rxs_burst...,], t, [Protein], [k, b])
@named bursty_rn_itp = ReactionSystem([rxs_burst...,], t, [Protein], [k, b, v...])

working_trn = 220
tslice = 5
nknots = 5

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

filter!(row -> 1 < row.Column1 < 50 && row.Column3 < working_trn - 1, df_div₊Sm)
filter!(row -> 1 < row.Column1 < 50 && row.Column3 < working_trn - 1, df_div₋Sm)

cycles₋Sm = rmnothing(vcat(make_cell_cycle.(lineages₋Sm)...));
fluors₋Sm = fluorescence_trajectory.(cycles₋Sm) ./ cf₋Sm
filter!(x -> 1 < length(x) < 50 && x[end] < working_trn - 1, fluors₋Sm)
fluors₋Sm_train, fluors₋Sm_test = partition(fluors₋Sm, 0.7; shuffle=true)

cycles₊Sm = rmnothing(vcat(make_cell_cycle.(lineages₊Sm)...));
fluors₊Sm = fluorescence_trajectory.(cycles₊Sm) ./ cf₋Sm
filter!(x -> 1 < length(x) < 50 && x[end] < working_trn - 1, fluors₊Sm)
fluors₊Sm_train, fluors₊Sm_test = partition(fluors₊Sm, 0.7; shuffle=true)

interdiv_times₋Sm = df_div₋Sm[!, :Column1] .- 0.5 
interdiv_times₊Sm = df_div₊Sm[!, :Column1] .- 0.5

interdiv_dist₋Sm = interdiv_dist_kde_pop(interdiv_times₋Sm; time_slice=5, bounds=(0.0, 300.0), bandwidth=5)
interdiv_dist₊Sm = interdiv_dist_kde_pop(interdiv_times₊Sm; time_slice=5, bounds=(0.0, 300.0), bandwidth=5)

# Fit the γ(τ).
# Cell inter-division time distribution. 
g₋Sm(t) = hazard(interdiv_dist₋Sm, t)
@register g₋Sm(t)
γτ₋Sm = g₋Sm(t)

g₊Sm(t) = hazard(interdiv_dist₊Sm, t)
@register g₊Sm(t)
γτ₊Sm = g₊Sm(t)

# Set-up model parameters.
function experiment_setup(;model_parameters, trn=(working_trn,), tspan=(0.0, 250.0), iters=20)
    return PopulationExperimentSetup(
        ps = model_parameters, # k,b 
        analytical_tspan = tspan,
        truncation = trn,
        iters = iters,
        )
end

x0 = fill(0.0, nknots)
quant₋Sm = quantile(last.(fluors₋Sm), [0.05, 0.95])
knots₋Sm = range(quant₋Sm[1], stop=quant₋Sm[2], length=nknots) 
itpp_spline₋Sm = Interp(x0, knots₋Sm)
Symbolics.@register_symbolic (itpp_spline₋Sm::Interp)(x, p1, p2, p3, p4, p5)

quant₊Sm = quantile(last.(fluors₊Sm), [0.05, 0.95])
knots₊Sm = range(quant₊Sm[1], stop=quant₊Sm[2], length=nknots) 
itpp_spline₊Sm = Interp(x0, knots₊Sm)
Symbolics.@register_symbolic (itpp_spline₊Sm::Interp)(x, p1, p2, p3, p4, p5)

sel_funhillfs = AgentBasedFSP.gen_division_rate_function(hillfs, bursty_rn)
