using Catalyst
using Distributions
using AgentBasedCells
using OrdinaryDiffEq
using Roots
solver_opts = (atol=1e-8, rtol=1e-8, stol=1e-4, method=Reinsert(), solver=TRBDF2(autodiff=false), rootfinder=Order2(), compute_joint=true)

include("$(srcdir())/plotting/plots.jl")
plot_options = (bins=10000, normalization=:pdf, strokewidth=0.0, x_gap=0.0, dodge_gap=0.0)

# Register the parameters and variables used.
@parameters t δ k K Ω λ
@species Protein(t)

# Positive feedback model.
rn = @reaction_network begin
    @parameters δ k K Ω λ
    Ω * δ + Ω * k * (Protein^2 / (K * Ω^2 + Protein^2)),  0 --> Protein
end 

f(a,θ,x) = exp(logpdf(Gamma(a,θ), x) - logccdf(Gamma(a,θ), x))
@register f(a,θ,x)

# Effective dilution model.
effective_dilution_rn = @reaction_network begin
    @parameters δ k K Ω λ 
    Ω * δ + Ω * k * (Protein^2 / (K * Ω^2 + Protein^2)),  0 --> Protein
    λ, Protein --> 0 
end 

## Parameter definitions for the experiments. 
function experiment_setup(;model_parameters, analytical_tspan, simulation_tspan, truncation, init, Ninit)
    return PopulationExperimentSetup(
        init = init,
        Δt = 1.0,
        max_pop = 100,
        Ninit=Ninit,
        ps = model_parameters, # α K Ω λ
        analytical_tspan = analytical_tspan,
        simulation_tspan = simulation_tspan,
        iters = 60, 
        truncation = truncation)
end

# First experiment
# Experiment suite.
δspan = 0.00:0.05:0.7
δspan_dense = 0.01:0.001:0.7
params = [(x, 6.0, 10.0, 10.0, 1.0) for x in δspan]
γwide = f(10, 0.1*0.7177346253629298, t) 

experiment_suite_mother = [
    experiment_setup(
        model_parameters = x, 
        analytical_tspan = (0.0, 2.5), 
        truncation=(150,), 
        init=[1.0,],
        Ninit=100,
        simulation_tspan=(0.0, 100.0)) for x in params]

experiment_suite = [
    experiment_setup(
        model_parameters = x, 
        analytical_tspan = (0.0, 2.5), 
        truncation=(150,), 
        init=[1.0,],
        Ninit=1,
        simulation_tspan=(0.0, 100.0)) for x in params]

## Run simulations & analytical.
model_mother = MotherCellModel(rn, DivisionRateBounded(γwide,14.0,rn), BinomialKernel(0.5))
analyticals_mother = run_analytical(model_mother, experiment_suite_mother; solver_opts...)
compute_joint!.(analyticals_mother)
failed_mother = findall(x -> x.flag==:Failed, analyticals_mother)
deleteat!(analyticals_mother, failed_mother)
deleteat!(experiment_suite, failed_mother)
marginal_size_distribution!.(analyticals_mother;)
division_dist_ancest!.(analyticals_mother)

# Dilution models
sdilution_models = 
    [StochasticDilutionModel(effective_dilution_rn, exp.init, exp.ps) for exp in experiment_suite]
birth_death_steady_state!.(sdilution_models, experiment_suite[1].truncation)

edilution_model = EffectiveDilutionModel(effective_dilution_rn, 1)
root_finding(edilution_model, experiment_suite[1].ps, δspan_dense; search_interval = IntervalBox(Interval(0.0,100.0), 1))

# Second experiment
model_pop = CellPopulationModel(rn, DivisionRateBounded(γwide,14.0,rn), BinomialKernel(0.5))
analyticals_pop = run_analytical(model_pop, experiment_suite; solver_opts...)
compute_joint!.(analyticals_pop)
failed_pop = findall(x -> x.flag==:Failed, analyticals_pop)
deleteat!(analyticals_pop, failed_pop)
deleteat!(experiment_suite, failed_pop)
marginal_size_distribution!.(analyticals_pop;)
division_dist_ancest!.(analyticals_pop)

simulation_mother = run_simulation(model_mother, experiment_suite_mother[9];)
simulation = run_simulation(model_pop, experiment_suite[9];)

DrWatson.default_prefix(e::PopulationExperimentSetup) = "Exp_txPfeedback_test"
plot_comparison(
    [simulation, simulation_mother],
    [analyticals_pop, analyticals_mother],
    edilution_model,
    sdilution_models,
    experiment_suite[1],
    bifurcation_parameter="α",
    tlims=(0.0, 2.0),
    tylims=(0.0, 2.),
    ylims=(0, 80),
    xmlims=(0,120),
    idxm=[4, 9, 12, 15],
    xticks=0:50:100,
    yticks=0:0.1:0.3,
    yticklabels=getindex(δspan, [5,9,12,15]),
    divtidx=1,
    offset=[0, 0.1],
    ridgelims=[120, 0.3, 0.05],
    acolors=[1,4],
    nlins=[5,5],
    trtlimx=(0.0, 100.0),
    trtlimy=(0.0, 120.0),
    sidx=1
   )

plot_comparison_extras([analyticals_pop, analyticals_mother],
    edilution_model,
    sdilution_models,
    experiment_suite[1],
    bifurcation_parameter="α",
    tlims=(0.0, 2.0),
    tticks=0.0:0.5:2.0,
    tylims=(0.0, 2.),
    ylims=(0, 70),
    xmlims=(0,120),
    idxm=[4, 9, 12, 15],
    xticks=0:50:100,
    yticks=0:0.1:0.3,
    yticklabels=getindex(δspan, [5,9,12,15]),
    divtidx=1,
    offset=[0, 0.1],
    ridgelims=[120, 0.3, 0.05],
    tridgelims=[0.2, 3.0],
    acolors=[1,4]
   )

initn = [1.0, 43.0]
experiment_suite_mother_init = [
    experiment_setup(
        model_parameters = params[9], 
        analytical_tspan = (0.0, 400.0), 
        init = [x, ],
        Ninit = 10,
        truncation=(150,), 
        simulation_tspan=(0.0, 100.0)) for x in initn]

experiment_suite_init = [
    experiment_setup(
        model_parameters = params[9], 
        analytical_tspan = (0.0, 100.0), 
        init = [x, ],
        Ninit=1,
        truncation=(100,), 
        simulation_tspan=(0.0, 100.0)) for x in initn]

simulations_mother_init = run_simulation.(Ref(model_mother), experiment_suite_mother_init;)
simulations_init = run_simulation.(Ref(model_pop), experiment_suite_init;)

plot_simulation_traj(simulations_init, simulations_mother_init, experiment_suite_init;
    scolors = [1,4],
    nlins = [10,10],
    tlims = (0.0, 50.0),
    ylims = (0, 120),
    initn = initn)

