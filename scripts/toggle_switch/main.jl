using Catalyst
using AgentBasedFSP
using FastGaussQuadrature
using OrdinaryDiffEq
using Roots
using LinearSolve
using Sundials
using Cubature
using Integrals
using Random
using FileIO
using JLD2

solver_settings = (
    atol=1e-6, 
    rtol=1e-6, 
    stol=1e-2, 
    method=Reinsert(), 
    integrator=QuadratureRule(FastGaussQuadrature.gausslegendre, n=100),
#    solver=CVODE_BDF(linear_solver = :GMRES, enable_sensitivities=false), 
    solver=CVODE_BDF(linear_solver = :GMRES), 
#    solver=QNDF(),
    rootfinder=Order2(), 
    compute_joint=true)

using Distributions
using IntervalArithmetic

include("$(srcdir())/steady_state/stochastic_dilution.jl")
include("$(srcdir())/steady_state/effective_dilution.jl")

include("$(srcdir())/plotting/plots.jl")

@variables t 
@species Protein1(t) Protein2(t) 

hillr_(p, L, K) = L / ((p/K)^2 + 1) # Repressing Hill function.
hill_(p, L, K) = L / ((K/p)^2 + 1) # Hill function.
@register hillr_(p, L, K)
@register hill_(p, L, K)

@parameters k1 k2 k3 α γ
@parameters K1 K2 K3 

rn = @reaction_network begin
    @parameters k1 K1 k2 K2 k3 K3 α γ
    @species Protein1(t) Protein2(t)
    hillr_(Protein2, k1, K1), 0 --> Protein1
    α + hillr_(Protein1, k2, K2), 0 --> Protein2
    γ, Protein1 --> 0
    γ, Protein2 --> 0
end

rn_eff_dil = @reaction_network begin
    @parameters k1 K1 k2 K2 k3 K3 α γ
    hillr_(Protein2, k1, K1), 0 --> Protein1
    α + hillr_(Protein1, k2, K2), 0 --> Protein2
    γ + hill_(Protein2, k3, K3), Protein1 --> 0
    γ + hill_(Protein2, k3, K3), Protein2 --> 0
end

γhill = hill_(Protein2, k3, K3)
γmax = k3

# Parameter definitions for the experiments. 
function experiment_setup(;model_parameters, 
        truncation=(170,80), 
        analytical_tspan=(0.0, 10.), 
        simulation_tspan=(0.0, 400.0),
        init = [0.0, 0.0],
        Ninit)

    return PopulationExperimentSetup(
        init = init,
        ps = model_parameters,
        analytical_tspan = analytical_tspan,
        iters = 50, 
        simulation_tspan = simulation_tspan,
        Δt = 0.1,
        max_pop = 1000,
        Ninit=Ninit,
        truncation = truncation)
end

#Ω = 2
#k1span = (8.0, 11.0)
#k1step = 0.5
#k1span_ = collect(k1span[1]:k1step:k1span[2])
#k1span_dense = k1span[1]:0.01:k1span[2]
#params = [last.((
#    k1 => k,      K1 => Ω*5.0, 
#    k2 => Ω*5.0,  K2 => Ω*5.0,
#    k3 => Ω*1.0,  K3 => Ω*10.0, α => Ω*0.1, γ => 0.05)) for k in k1span_]

#Ω = 3
#k1span = (13.0, 17.0)
#k1step = 0.5
#k1span_ = collect(k1span[1]:k1step:k1span[2])
#k1span_dense = k1span[1]:0.01:k1span[2]
#params = [last.((
#    k1 => k,      K1 => Ω*5.0, 
#    k2 => Ω*5.0,  K2 => Ω*5.0,
#    k3 => Ω*1.0,  K3 => Ω*10.0, α => Ω*0.1, γ => 0.1)) for k in k1span_]

Ω = 4
k1span = (16.0, 22.0)
k1step = 0.5
k1span_ = collect(k1span[1]:k1step:k1span[2])
k1span_dense = k1span[1]:0.01:k1span[2]
params = [last.((
    k1 => k,      K1 => Ω*5.0, 
    k2 => Ω*5.0,  K2 => Ω*5.0,
    k3 => Ω*1.0,  K3 => Ω*30.0, α => Ω*0.1, γ => 0.2)) for k in k1span_]
model_pop = CellPopulationModel(rn, DivisionRateBounded(γhill, γmax, rn), BinomialKernel(0.5))
model_mother = MotherCellModel(rn, DivisionRateBounded(γhill, γmax, rn), BinomialKernel(0.5))

experiment_suite_pop = [
    experiment_setup(model_parameters = x, Ninit=1) for x in params]
analyticals_pop = run_analytical(model_pop, experiment_suite_pop; solver_settings...)
failed_pop = findall(x -> x.flag==:Failed, analyticals_pop)
deleteat!(analyticals_pop, failed_pop)
deleteat!(experiment_suite_pop, failed_pop)
marginal_size_distribution!.(analyticals_pop; atol=1e-4)
division_dist_ancest!.(analyticals_pop)

experiment_suite_mother = [
    experiment_setup(model_parameters = x, 
    Ninit=100, 
    analytical_tspan=(0.0, 100.)) for x in params]
analyticals_mother = run_analytical(model_mother, experiment_suite_mother; solver_settings...)
failed_mother = findall(x -> x.flag==:Failed, analyticals_mother)
deleteat!(analyticals_mother, failed_mother)
deleteat!(experiment_suite_mother, failed_mother)
marginal_size_distribution!.(analyticals_mother; atol=1e-4)
division_dist_ancest!.(analyticals_mother)

# Dilution models
sdil_models = 
    [StochasticDilutionModel(rn_eff_dil, exp.init, exp.analytical_tspan, exp.ps) for exp in experiment_suite_pop]
stochastic_steady_state!.(
    sdil_models, 
    Ref((700, 100)), 
    Ref(CVODE_BDF(linear_solver = :GMRES)))
simulation_steady_state!.(sdil_models, Ref([50.0, 5.0]), Ref((0.0, 1.0e6)))

edil_model = EffectiveDilutionModel(rn_eff_dil, 1)
root_finding(edil_model, experiment_suite_pop[1].ps, k1span_dense; search_interval = IntervalBox(interval(0.0,200.0), 2))

params_ = last.((
    k1 => 18.9, K1 => Ω*5.0, 
    k2 => Ω*5.0, K2 => Ω*5.0,
    k3 => Ω*1.0, K3 => Ω*30.0, α => Ω*0.1, γ => 0.2)) 

experiment_ = experiment_setup(model_parameters = params_, Ninit=1000, analytical_tspan=(0.0, 100.), simulation_tspan=(0.0,400.0))
experiment_p = experiment_setup(model_parameters = params_, Ninit=1, analytical_tspan=(0.0, 100.), simulation_tspan=(0.0,400.0))
simulation_mother = run_simulation(model_mother, experiment_;)
simulation = run_simulation(model_pop, experiment_p;)

DrWatson.default_prefix(e::PopulationExperimentSetup) = "Exp_toggle_switch"
seedn = 4 # Needed to have the sample trajectories for both plots from the same simulation.

Random.seed!(seedn)
idxs = [1,5,8,13]
plot_comparison(
    [simulation, simulation_mother],
    [analyticals_pop, analyticals_mother],
    edil_model,
    sdil_models,
    experiment_suite_mother[1],
    bifurcation_parameter="α",
    ylims=(0.0, 100.0),
    tlims=(0.0, 1.0),
    tylims=(0.0, 0.6),
    divtidx=10,
    xmlims=(0,140),
    idxm=idxs,
    xticks=0:25:140,
    yticks=0:0.3:0.6,
    yticklabels=getindex(k1span_, idxs),
    offset=[0, 0.2],
    ridgelims=[140, 0.6, 0.1],
    acolors=[1,2],
    nlins=[1,1],
    trtlimx=(0.0, 400.0),
    trtlimy=(0.0, 120.0),
    sidx=1,
    heightmod = 0.5
    )
plot_comparison_extras([analyticals_pop, analyticals_mother],
    edil_model,
    sdil_models,
    experiment_suite_mother[1],
    bifurcation_parameter="α",
    ylims=(0.0, 70.0),
    tlims=(0.0, 10.0),
    tticks=0.0:1.0:10.0,
    tylims=(0.0, 0.1),
    divtidx=10,
    xmlims=(0,75),
    idxm=idxs,
    xticks=0:10:170,
    yticks=range(0, stop=0.6, length=length(idxs)),
    yticklabels=getindex(k1span_, idxs),
    offset=[0, 0.2],
    ridgelims=[75, 0.6, 0.07],
    tridgelims=[0.1, 0.5],
    acolors=[1,2],
    sidx = 1
    )


DrWatson.default_prefix(e::PopulationExperimentSetup) = "Exp_toggle_switch_2"
Random.seed!(seedn)
plot_comparison(
    [simulation, simulation_mother],
    [analyticals_pop, analyticals_mother],
    edil_model,
    sdil_models,
    experiment_suite_mother[1],
    bifurcation_parameter="α",
    ylims=(0.0, 50.0),
    tlims=(0.0, 1.0),
    tylims=(0.0, 0.6),
    divtidx=10,
    xmlims=(0,50),
    idxm=idxs,
    xticks=0:5:40,
    yticks=0:0.3:0.6,
    yticklabels=getindex(k1span_, idxs),
    offset=[0, 0.2],
    ridgelims=[50, 0.6, 0.2],
    acolors=[1,2],
    nlins=[1,1],
    trtlimx=(0.0, 400.0),
    trtlimy=(0.0, 60.0),
    sidx=2,
    heightmod = 0.5,
    linestyle = :solid
    )
plot_comparison_extras([analyticals_pop, analyticals_mother],
    edil_model,
    sdil_models,
    experiment_suite_mother[1],
    bifurcation_parameter="α",
    ylims=(0.0, 70.0),
    tlims=(0.0, 10.0),
    tticks=0.0:1.0:10.0,
    tylims=(0.0, 0.1),
    divtidx=10,
    xmlims=(0,75),
    idxm=idxs,
    xticks=0:10:170,
    yticks=range(0, stop=0.6, length=length(idxs)),
    yticklabels=getindex(k1span_, idxs),
    offset=[0, 0.2],
    ridgelims=[75, 0.6, 0.07],
    tridgelims=[0.1, 0.5],
    acolors=[1,2],
    sidx = 2
    )

params_ = last.((
    k1 => 18.9, K1 => Ω*5.0, 
    k2 => Ω*5.0, K2 => Ω*5.0,
    k3 => Ω*1.0, K3 => Ω*30.0, α => Ω*0.1, γ => 0.2)) 

# Plot distributions
experiment_ = experiment_setup(model_parameters = params_, Ninit=1, analytical_tspan=(0.0, 100.))
analytical_ = run_analytical_single(model_mother, experiment_; solver_settings...)
marginal_size_distribution!(analytical_; atol=1e-4)
division_dist_ancest!(analytical_)

analytical_pop = run_analytical_single(model_pop, experiment_; solver_settings...)
marginal_size_distribution!(analytical_pop; atol=1e-4)
division_dist_ancest!(analytical_pop)

size_inches = (17, 5.5)
size_pt = 72 .* size_inches
fig = Figure(resolution=size_pt, fontsize=16, figure_padding=1)

fig = Figure(resolution = size_pt) 
axs = [Axis(fig[1,i], xlabel="Protein x count", ylabel="Protein y count", title="Placeholder") for i in 1:3]

cmaps = [ColorScheme(range(alphacolor(colors[i], 0.2), stop=ARGB(colors[i], 1.0), length=10)) for i in 1:length(colors)]
cmapsf = [ColorScheme(range(alphacolor(colors[i], 0.0), stop=ARGB(colors[i], 1.0), length=10)) for i in 1:length(colors)]

levels_ = range(1e-8, 0.005, length=5)

contour!(axs[1], analytical_pop.results[:marginal_size]; levels=14, colormap=cmaps[1], linewidth=3)
contour!(axs[1], analytical_.results[:marginal_size]; levels=14, colormap=cmaps[2], linewidth=3)

contour!(axs[2], analytical_.results[:birth_dist]; levels=14, colormap=cmaps[2], linewidth=3, lowclip=1e-8)
contour!(axs[3], analytical_.results[:division_dist]; levels=14, colormap=cmaps[2], linewidth=3, lowclip=1e-8)

contour!(axs[2], analytical_pop.results[:birth_dist]; levels=14, colormap=cmaps[1], linewidth=3, lowclip=1e-8)
contour!(axs[3], analytical_pop.results[:division_dist]; levels=14, colormap=cmaps[1], linewidth=3, lowclip=1e-8)

axs[1].title = "Marginal snapshot distribution"
axs[2].title = "Birth protein distribution"
axs[3].title = "Division protein distribution"
ylims!(axs[1], (0, 70))
xlims!(axs[1], (0, 120))

ylims!(axs[2], (0, 50))
xlims!(axs[2], (0, 60))

ylims!(axs[3], (0, 80))
xlims!(axs[3], (0, 120))

hidedecorations!(axs[1], ticks=false, label=false, ticklabels=false)
hidedecorations!(axs[2], ticks=false, label=false, ticklabels=false)
hidedecorations!(axs[3], ticks=false, label=false, ticklabels=false)
hidespines!(axs[1], :r, :t)
hidespines!(axs[2], :r, :t)
hidespines!(axs[3], :r, :t)

save("$(plotsdir())/toggle_distribution.pdf", fig)

# Appendix simulation trajectories.
initn = [[75.0, 5.0], [0, 35]]
experiment_suite_mother_init = [
    experiment_setup(
        model_parameters = params_, 
        init = x,
        Ninit = 100,
        simulation_tspan=(0.0, 400.0)) for x in initn]

experiment_suite_init = [
    experiment_setup(
        model_parameters = params_, 
        init = x,
        Ninit=1,
        simulation_tspan=(0.0, 400.0)) for x in initn]

simulations_mother_init = run_simulation.(Ref(model_mother), experiment_suite_mother_init;)
simulations_init = run_simulation.(Ref(model_pop), experiment_suite_init;)

DrWatson.default_prefix(e::PopulationExperimentSetup) = "Exp_toggle_switch"
Random.seed!(28)
plot_simulation_traj(simulations_init, simulations_mother_init, experiment_suite_init;
    scolors = [1,2],
    nlins = [10,10],
    tlims = (0.0, 400.0),
    ylims = (0, 150),
    initn = initn,
    sidx=1,
) 

DrWatson.default_prefix(e::PopulationExperimentSetup) = "Exp_toggle_switch_2"
Random.seed!(28)
plot_simulation_traj(simulations_init, simulations_mother_init, experiment_suite_init;
    scolors = [1,2],
    nlins = [10,10],
    tlims = (0.0, 400.0),
    ylims = (0, 90),
    initn = initn,
    sidx=2,
) 

