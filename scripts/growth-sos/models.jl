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

using StatProfilerHTML
include("$(srcdir())/fitting_helpers.jl")

function interdiv_dist_kde(interdiv_times; time_slice, bounds, bandwidth)
    x_ = interdiv_times .* time_slice .- 0.5 * time_slice
    return CumulativeKDE(InterpKDE(kde(x_; bandwidth=bandwidth, boundary=bounds)), tend=bounds[end])
end

# Model.
# Register the parameters and variables used.
@parameters t k b L K δ n m K2 d
@species Protein(t)

hillfs0 = Protein^n / (K^n + Protein^n)
hillfs = L*((1-δ)*Protein^m / (K2^m + Protein^m) + δ)

# Bursty model. Note that maximum burst size approximation is made. Better
# compute for the upper bound.
rxs_burst = [Reaction(k*(b-1)^(n-1)/b^(n+1), nothing, [Protein], nothing, [n]) for n in 1:50]
# We are fitting (k'b, b)
rx_deg = Reaction(d, [Protein], nothing)

@named bursty_rn = ReactionSystem([rxs_burst...,], t, [Protein], [k, b, L, K2, δ, m])
@named bursty_deg_rn = ReactionSystem([rxs_burst..., rx_deg], t, [Protein], [k, b, d, L, K2, δ, m])
@named bursty_rn_fs0 = ReactionSystem([rxs_burst...,], t, [Protein], [k, b, n, K, L, K2, δ, m])
@named bursty_rn_react = ReactionSystem([rxs_burst...,], t, [Protein], [k, b])

working_trn = 200

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

cycles_2pal = cycles_2pal ./ cf_wt
filter!(x -> 0 < length(x) < 40 && x[end] < working_trn - 1, cycles_2pal)
cycles_wt = cycles_wt ./ cf_wt
filter!(x -> 0 < length(x) < 40 && x[end] < working_trn - 1, cycles_wt)

tslice = 8

filter!(row -> 0 < row.Column1 < 40 && row.Column3 < 150, df_div_wt)
filter!(row -> 0 < row.Column1 < 40 && row.Column3 < 150, df_div_2pal)

interdiv_times_wt = df_div_wt[!, :Column1]
interdiv_times_2pal = df_div_2pal[!, :Column1]

interdiv_dist_wt = interdiv_dist_kde(interdiv_times_wt; time_slice=tslice, bandwidth=10.0, bounds=(0.0, 400.0))
interdiv_dist_2pal = interdiv_dist_kde(interdiv_times_2pal; time_slice=tslice, bandwidth=10.0, bounds=(0.0, 400.0))

# Fit the γ(τ).
# Cell inter-division time distribution. 
g_wt(t) = hazard(interdiv_dist_wt, t)
@register g_wt(t)
γτ_wt = g_wt(t)

g_2pal(t) = hazard(interdiv_dist_2pal, t) 
@register g_2pal(t)
γτ_2pal = g_2pal(t)

# Set-up model parameters.
function experiment_setup(;model_parameters, trn=working_trn)
    return PopulationExperimentSetup(
        ps = model_parameters, # k,b 
        analytical_tspan = (0.0, 340.0),
        truncation = (trn,),
        iters = 20,
        )
end
