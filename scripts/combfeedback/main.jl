using Catalyst
using Distributions
using AgentBasedCells
using OrdinaryDiffEq
using Roots
using FastGaussQuadrature
using Integrals
using Sundials

solver_opts = (
    atol=1e-6, 
    rtol=1e-6, 
    stol=1e-2, 
    method=Reinsert(), 
    #solver=CVODE_BDF(linear_solver = :GMRES),
    solver=TRBDF2(autodiff=false), 
    rootfinder=Order16(), 
    compute_joint=true, 
    integrator=QuadratureRule(FastGaussQuadrature.gausslegendre, n=1000))

include("$(srcdir())/steady_state/effective_dilution.jl")
include("$(srcdir())/steady_state/stochastic_dilution.jl")
include("$(srcdir())/plotting/plots.jl")

plot_options = (bins=10000, normalization=:pdf, strokewidth=0.0, x_gap=0.0, dodge_gap=0.0)

# Register the parameters and variables used.
@parameters t δ Ω k1 K1 λ k2 K2 μ
@species Protein(t)

# Positive feedback model.
rn = @reaction_network begin
    @parameters δ Ω k1 K1 λ k2 K2 μ
    2*δ*Ω + Ω*k2 * (Protein^2 / (K2 * Ω^2 + Protein^2)),  0 --> Protein
    λ+μ, Protein --> 0
end 

γx = k1*(K1 * Ω^2 / (K1 * Ω^2 + Protein^4))

# Effective dilution model.
effective_dilution_rn = @reaction_network begin
    @parameters δ Ω k1 K1 λ k2 K2 μ
    2*δ*Ω + Ω*k2 * (Protein^2 / (K2 * Ω^2 + Protein^2)),  0 --> Protein
    λ+μ, Protein --> 0
    k1*(K1 * Ω^2 / (K1 * Ω^2 + Protein^4)), Protein --> 0 
end 

## Parameter definitions for the experiments. 
function experiment_setup(;model_parameters, analytical_tspan, simulation_tspan, truncation, Ninit, init)
    return PopulationExperimentSetup(
        init = init,
        Δt = 0.1,
        max_pop = 101,
        Ninit=Ninit,
        ps = model_parameters, # α K Ω λ
        analytical_tspan = analytical_tspan,
        simulation_tspan = simulation_tspan,
        iters = 60, 
        truncation = truncation)
end

# Experiment suite.
δspan = 4.2:0.1:5.0
δspan_dense = 4.2:0.01:5.0
#@parameters δ Ω k1 K1 λ k2 K2 μ
params = [(x, 70.0, 40.0, 15.0, 5.0, 80.0, 4.0, 20.0) for x in δspan]
#params = [(x, 45.0, 30.0, 10.0, 5.0, 80.0, 4.0, 20.0) for x in δspan]
#params = [(x, 42.0, 40.0, 5.0, 5.0, 80.0, 4.0, 20.0) for x in δspan]

experiment_suite_mother = [
    experiment_setup(
        model_parameters = x, 
        analytical_tspan = (0.0, 100.0), 
        truncation=(300,), 
        Ninit=11,
        init = [0.0,],
        simulation_tspan=(0.0, 10.0)) for x in params]

experiment_suite = [
    experiment_setup(
        model_parameters = x, 
        analytical_tspan = (0.0, 10.0), 
        truncation=(300,), 
        Ninit=1,
        init = [0.0, ],
        simulation_tspan=(0.0, 10.0)) for x in params]

model_mother = MotherCellModel(rn, DivisionRateBounded(γx,18.1,rn), BinomialKernel(0.5))
model_pop = CellPopulationModel(rn, DivisionRateBounded(γx,18.1,rn), BinomialKernel(0.5))

# Run simulations & analytical.
analyticals_mother = run_analytical(model_mother, experiment_suite_mother; solver_opts...)
compute_joint!.(analyticals_mother)
failed_mother = findall(x -> x.flag==:Failed, analyticals_mother)
deleteat!(analyticals_mother, failed_mother)
deleteat!(experiment_suite, failed_mother)
marginal_size_distribution!.(analyticals_mother;)
division_dist_ancest!.(analyticals_mother)

# Dilution models
sdilution_models = 
    [StochasticDilutionModel(effective_dilution_rn, exp.init, exp.analytical_tspan, exp.ps) for exp in experiment_suite]
birth_death_steady_state!.(sdilution_models, experiment_suite[1].truncation)

edilution_model = EffectiveDilutionModel(effective_dilution_rn, 1)
root_finding(
    edilution_model, 
    experiment_suite[1].ps, 
    δspan_dense; 
    search_interval = IntervalBox(interval(0.0,200.0), 1))

# Second experiment
analyticals_pop = run_analytical(model_pop, experiment_suite; solver_opts...)
compute_joint!.(analyticals_pop)
failed_pop = findall(x -> x.flag==:Failed, analyticals_pop)
deleteat!(analyticals_pop, failed_pop)
deleteat!(experiment_suite, failed_pop)
marginal_size_distribution!.(analyticals_pop;)
division_dist_ancest!.(analyticals_pop)

simulation_mother = run_simulation(model_mother, experiment_suite_mother[7];)
simulation = run_simulation(model_pop, experiment_suite[7];)

DrWatson.default_prefix(e::PopulationExperimentSetup) = "Exp_combfeedback"
plot_comparison(
    [simulation, simulation_mother],
    [analyticals_pop, analyticals_mother],
    edilution_model,
    sdilution_models,
    experiment_suite[1],
    bifurcation_parameter="α",
    tlims=(0.0, 2.0),
    tylims=(0.0, 2.0),
    ylims=(0,250),
    xmlims=(0,250),
    idxm=[1,3,6,9],
    xticks=0:50:250,
    yticks=0:0.1:0.3,
    yticklabels=getindex(δspan, [1,3,6,9]),
    divtidx=1,
    offset=[0, 0.1],
    ridgelims=[250, 0.3, 0.02],
    acolors=[1,4],
    nlins=[5,5],
    trtlimx=(0.0, 200.0),
    trtlimy=(0.0, 200.0),
    sidx=1
   )

plot_comparison_extras([analyticals_pop, analyticals_mother],
    edilution_model,
    sdilution_models,
    experiment_suite[1],
    bifurcation_parameter="α",
    tlims=(0.0, 0.5),
    tticks=0.0:0.1:0.5,
    tylims=(0.0, 2.),
    ylims=(0, 70),
    xmlims=(0,250),
    idxm=[1, 3, 4, 5],
    xticks=0:50:250,
    yticks=0:0.1:0.3,
    yticklabels=getindex(δspan, [1,3,6,9]),
    divtidx=1,
    offset=[0, 0.1],
    ridgelims=[250, 0.3, 0.02],
    tridgelims=[0.2, 3.0],
    acolors=[1,4])

# Appendix simulations
initn = [1.0, 37.0, 111.0]
experiment_suite_mother_init = [
    experiment_setup(
        model_parameters = params[4], 
        analytical_tspan = (0.0, 10.0), 
        init = [x, ],
        Ninit = 10,
        truncation=(150,), 
        simulation_tspan=(0.0, 10.0)) for x in initn]

experiment_suite_init = [
    experiment_setup(
        model_parameters = params[4], 
        analytical_tspan = (0.0, 10.0), 
        init = [x, ],
        Ninit=1,
        truncation=(100,), 
        simulation_tspan=(0.0, 5.0)) for x in initn]

simulations_mother_init = run_simulation.(Ref(model_mother), experiment_suite_mother_init;)
simulations_init = run_simulation.(Ref(model_pop), experiment_suite_init;)

plot_simulation_traj(simulations_init, simulations_mother_init, experiment_suite_init;
    scolors = [1,4],
    nlins = [10,10],
    tlims = (0.0, 2.0),
    ylims = (0, 200),
    initn = initn
   )
