using Catalyst
using Distributions
using AgentBasedFSP
using StatProfilerHTML
using OrdinaryDiffEq
using Roots
solver_opts = (atol=1e-6, rtol=1e-6, stol=1e-4, method=Reinsert(), solver=TRBDF2(), rootfinder=Order16(), compute_joint=true)

# Import reusable scripts. 
include("$(srcdir())/plotting/plots.jl")
include("$(srcdir())/steady_state/stochastic_dilution.jl")
#include("$(srcdir())/steady_state/.jl")

# Register the parameters and variables used.
@parameters t K Ω λ k k0
@variables Protein(t)

# Positive feedback model.
rn = @reaction_network begin
    @parameters α K Ω λ k k0
    Ω * α, 0 --> Protein
    λ, Protein --> 0
end 

# Effective diffusion model.
effective_dilution_rn = @reaction_network begin
    @parameters α K Ω λ k k0
    Ω * α, 0 --> Protein
    k0 * ((1-k) * K * Ω^2 / (K * Ω^2 + Protein^3) + k) + λ, Protein --> 0
end 

γx = k0*((1-k) * K * Ω^2  / (K * Ω^2 + Protein^3) + k)

## Parameter definitions for the experiments. 
function experiment_setup(;model_parameters, analytical_tspan, truncation, simulation_tspan, init, Ninit)
    return PopulationExperimentSetup(
        init = init,
        ps = model_parameters, # α K Ω λ
        analytical_tspan = analytical_tspan,
        iters = 60, 
        simulation_tspan = simulation_tspan,
        Δt = 1.0,
        max_pop = 100,
        Ninit=Ninit,
        truncation = truncation)
end

model_mother = MotherCellModel(rn, DivisionRateBounded(γx, k0, rn), BinomialKernel(0.5))
model_pop = CellPopulationModel(rn, DivisionRateBounded(γx, k0, rn), BinomialKernel(0.5))

# Experiment suite.
δ1 = 1.5 
δ2 = 6.0 
δspan = δ1:0.5:δ2
δspan_dense = δ1:0.001:δ2
params = [(x, 1.0, 10.0, 1.0, 0.0, 20.0) for x in δspan]

experiment_suite_mother = [
    experiment_setup(
        model_parameters = x, 
        analytical_tspan = (0.0, 400.0), 
        init = [1.0,],
        Ninit = 100,
        truncation=(150,), 
        simulation_tspan=(0.0, 100.0)) for x in params]

experiment_suite = [
    experiment_setup(
        model_parameters = x, 
        analytical_tspan = (0.0, 100.0), 
        init = [1.0,],
        Ninit = 1,
        truncation=(100,), 
        simulation_tspan=(0.0, 20.0)) for x in params]

# Run simulations & analytical.
analyticals_gr_mother = run_analytical(model_mother, experiment_suite_mother; solver_opts...)
compute_joint!.(analyticals_gr_mother)
failed_mother = findall(x -> x.flag==:Failed, analyticals_gr_mother)
deleteat!(analyticals_gr_mother, failed_mother)
deleteat!(experiment_suite_mother, failed_mother)
division_dist_ancest!.(analyticals_gr_mother)
marginal_size_distribution!.(analyticals_gr_mother; atol=1e-6, rtol=1e-6)

# Dilution models
sdilution_gr_models = 
    [StochasticDilutionModel(effective_dilution_rn, exp.init, exp.ps) for exp in experiment_suite]
birth_death_steady_state!.(sdilution_gr_models, experiment_suite[1].truncation)

edilution_gr_model = EffectiveDilutionModel(effective_dilution_rn, 1)
root_finding(edilution_gr_model, experiment_suite[1].ps, δspan_dense; search_interval = Interval(0.0, 500.0))

analyticals_gr_pop = run_analytical(model_pop, experiment_suite; solver_opts...)
compute_joint!.(analyticals_gr_pop)
failed_pop = findall(x -> x.flag==:Failed, analyticals_gr_pop)
deleteat!(analyticals_gr_pop, failed_pop)
deleteat!(experiment_suite, failed_pop)
marginal_size_distribution!.(analyticals_gr_pop; atol=1e-6, rtol=1e-6)
division_dist_ancest!.(analyticals_gr_pop)

simulation_mother = run_simulation(model_mother, experiment_suite_mother[2];)
simulation = run_simulation(model_pop, experiment_suite[2];)

DrWatson.default_prefix(e::PopulationExperimentSetup) = "Exp_gPfeedback_selection_alt"
plot_comparison(
    [simulation, simulation_mother],
    [analyticals_gr_pop, analyticals_gr_mother],
    edilution_gr_model,
    sdilution_gr_models,
    experiment_suite[1],
    bifurcation_parameter="α",
    ylims=(0.0, 70.0),
    tlims=(0.0, 1.0),
    tylims=(0.0, 0.6),
    divtidx=10,
    xmlims=(0,75),
    idxm=[1,3,5,8],
    xticks=0:25:75,
    yticks=0:0.2:0.6,
    yticklabels=getindex(δspan, [1,3,5,8]),
    offset=[0, 0.2],
    ridgelims=[75, 0.6, 0.07],
    acolors=[1,4],
    nlins=[10,10],
    trtlimx=(0.0, 100.0),
    trtlimy=(0.0, 30.0)
    )

plot_comparison_extras([analyticals_gr_pop, analyticals_gr_mother],
    edilution_gr_model,
    sdilution_gr_models,
    experiment_suite[1],
    bifurcation_parameter="α",
    ylims=(0.0, 70.0),
    tlims=(0.0, 1.0),
    tticks=0.0:0.2:1.0,
    tylims=(0.0, 0.6),
    divtidx=10,
    xmlims=(0,75),
    idxm=[1,3,5,8],
    xticks=0:25:75,
    yticks=0:0.2:0.6,
    yticklabels=getindex(δspan, [1,3,5,8]),
    offset=[0, 0.2],
    ridgelims=[75, 0.6, 0.07],
    tridgelims=[0.2, 3.0],
    acolors=[1,4]
    )


initn = [1.0, 20.0]
experiment_suite_mother_init = [
    experiment_setup(
        model_parameters = params[2], 
        analytical_tspan = (0.0, 400.0), 
        init = [x, ],
        Ninit = 100,
        truncation=(150,), 
        simulation_tspan=(0.0, 100.0)) for x in initn]

experiment_suite_init = [
    experiment_setup(
        model_parameters = params[2], 
        analytical_tspan = (0.0, 100.0), 
        init = [x, ],
        Ninit=10,
        truncation=(100,), 
        simulation_tspan=(0.0, 2.0)) for x in initn]

simulations_mother_init = run_simulation.(Ref(model_mother), experiment_suite_mother_init;)
simulations_init = run_simulation.(Ref(model_pop), experiment_suite_init;)

plot_simulation_traj(simulations_init, simulations_mother_init, experiment_suite_init;
    scolors = [1,4],
    nlins = [10,50],
    tlims = (0.0, 2.0),
    ylims = (0, 35),
    initn = initn) 

