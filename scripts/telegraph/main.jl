using ModelingToolkit
using Catalyst
using AgentBasedCells
using OrdinaryDiffEq
solver_settings = (atol=1e-6, rtol=1e-6, stol=1e-6, method=Reinsert(), solver=TRBDF2())
using FiniteStateProjection
using Distributions
using Distributions: Geometric
using FileIO

include("$(srcdir())/steady_state/effective_dilution.jl")
include("$(srcdir())/steady_state/stochastic_dilution.jl")
include("$(srcdir())/plotting/plots.jl")

# Model
@parameters t k b L K δ n
@species I(t), A(t), Protein(t)

rn = @reaction_network begin
    @parameters λ μ k b d L K δ n
    λ  * (1 - A), 0 --> A
    μ, A --> 0
    d, Protein --> 0 
end 

rxs = [Reaction(k*(b-1)^(n-1)/b^n, [A], [Protein, A], [1], [n,1]) for n in 1:20]
addreaction!.(Ref(rn), rxs)

haz(a,θ,x) = exp(logpdf(Gamma(a,θ), x) - logccdf(Gamma(a,θ), x))
@register haz(a,θ,x)

γ = haz(1/0.6, 0.6*0.859527611609782, t)
hill = L*((1-δ)*Protein^n / (K^n + Protein^n) + δ)
γmax = 2

# Parameter definitions for the experiments. 
function experiment_setup(;model_parameters, 
        truncation=(2,120), 
        analytical_tspan=(0.0, 5.0), 
        simulation_tspan=(0.0, 15.0))

    return PopulationExperimentSetup(
        init = [1.0, 0.0],
        ps = model_parameters,
        analytical_tspan = analytical_tspan,
        iters = 50, 
        simulation_tspan = simulation_tspan,
        Δt = 0.1,
        max_pop = 3000,
        truncation = truncation)
end

colors = ColorSchemes.Egypt.colors

DrWatson.default_prefix(e::PopulationExperimentSetup) = "Exp_telegraph"
params = (0.7, 0.2, 20.0, 1.2, 0.1, 1.0, 20.0, 0.0, 2.0)
exper = experiment_setup(model_parameters = params)

model = CellPopulationModel(rn, DivisionRateBounded(γ, γmax, rn), BinomialWithDuplicate(0.5, [1]))
analytical = run_analytical_single(model, exper; solver_settings...)
division_dist_ancest!(analytical)
simulation = run_simulation(model, exper;)

model_sel = CellPopulationModel(rn, DivisionRateBounded(γ*hill, γmax*2, rn), BinomialWithDuplicate(0.5, [1])) 
analytical_sel = run_analytical_single(model_sel, exper; solver_settings...)
division_dist_ancest!(analytical_sel)
simulation_sel = run_simulation(model_sel, exper;)

exper_mother = experiment_setup(;model_parameters = params, simulation_tspan=(0.0, 15000.0))
model_mother = MotherCellModel(rn, DivisionRateBounded(γ, γmax*2, rn), BinomialWithDuplicate(0.5, [1]))
analytical_mother = run_analytical_single(model_mother, exper_mother; solver_settings...)
division_dist_ancest!(analytical_mother)
simulation_mother = run_simulation(model_mother, exper_mother;)

exper_mother_sel = experiment_setup(;model_parameters = params, simulation_tspan=(0.0, 15000.0))
model_mother_sel = MotherCellModel(rn, DivisionRateBounded(γ*hill, γmax*2, rn), BinomialWithDuplicate(0.5, [1]))
analytical_mother_sel = run_analytical_single(model_mother_sel, exper_mother_sel; solver_settings...)
division_dist_ancest!(analytical_mother_sel)
simulation_mother_sel = run_simulation(model_mother_sel, exper_mother;)

truncations = 20:10:100
exper_suite_trn = [experiment_setup(model_parameters = params, truncation=(2,trn)) for trn in truncations]
analytical_trn = run_analytical(model, exper_suite_trn; solver_settings...)
analytical_sel_trn = run_analytical(model_sel, exper_suite_trn; solver_settings...)
analytical_mother_trn = run_analytical(model_mother, exper_suite_trn; solver_settings...)
analytical_mother_sel_trn = run_analytical(model_mother_sel, exper_suite_trn; solver_settings...)

exper_suite_trn_sh = [experiment_setup(model_parameters = params, truncation=(2,trn), analytical_tspan=(0.0, 2.0)) for trn in truncations]
analytical_trn_sh = run_analytical(model, exper_suite_trn_sh; solver_settings...)
analytical_sel_trn_sh = run_analytical(model_sel, exper_suite_trn_sh; solver_settings...)
analytical_mother_trn_sh = run_analytical(model_mother, exper_suite_trn_sh; solver_settings...)
analytical_mother_sel_trn_sh = run_analytical(model_mother_sel, exper_suite_trn_sh; solver_settings...)

ttruncs = 0.5:0.5:7.0 
exper_suite_trn = [
    experiment_setup(model_parameters = params,
                     analytical_tspan=(0.0, tend), truncation=(2,trn)) for (trn, tend) in Iterators.product(truncations, ttruncs)]
analytical_trn_short = run_analytical(model, exper_suite_trn; solver_settings...)
analytical_sel_trn_short = run_analytical(model_sel, exper_suite_trn; solver_settings...)
analytical_mother_trn_short = run_analytical(model_mother, exper_suite_trn; solver_settings...)
analytical_mother_sel_trn_short = run_analytical(model_mother_sel, exper_suite_trn; solver_settings...)

save("$(projectdir())/scripts/telegraph/workspace.jld2", 
    simulation, simulation_sel, simulation_mother, simulation_mother_sel,
    analytical, analytical_sel, analytical_mother, analytical_mother_sel, 
    analytical_trn, analytical_sel_trn, analytical_trn_short, analytical_sel_trn_short,
    analytical_mother_trn, analytical_mother_sel_trn, analytical_mother_trn_short, analytical_sel_trn_short)

plots_selection_compare(
    [simulation, simulation_sel], 
    [analytical, analytical_sel], 
    [analytical_mother, analytical_mother_sel], 
    exper; 
    idx=2, 
    xlimsbirth=(0.0, 30.0), 
    xlimsdiv=(0.0, 80.0), 
    ylimsbirth=(0.0, 0.1), 
    ylimsdiv=(0.0, 0.04), 
    divtbins=0.0:0.1:5.0, 
    divtrange=0.0:0.01:3.5, 
    fptlimsy=(0.0, 50), 
    fptlimsx=(0.0, 2.0),
    selection=nothing, 
    tag="selelection", 
    sumdim=1, 
    acolors=[2,1],
    acolorsm=[3,4],
    scolors=[2,1],
    traj_idxs=[])

plot_convergence_extras(
    [analytical, analytical_sel], 
    [analytical_trn, analytical_sel_trn], 
    [analytical_trn_sh, analytical_sel_trn_sh], 
    [analytical_trn_short[:, :], analytical_sel_trn_short[:, :]],
    exper; 
    simulations=[simulation, simulation_sel],
    idxs = [[1,4,11], [1,4,11]],
    acolors = [2,1],
    labels=["No selection", "Selection"],
    mlabel="Population",
    note="pop", drawgrowth=true)

plot_convergence_extras(
    [analytical_mother, analytical_mother_sel], 
    [analytical_mother_trn, analytical_mother_sel_trn], 
    [analytical_mother_trn_sh, analytical_mother_sel_trn_sh], 
    [analytical_mother_trn_short[:, :], analytical_mother_sel_trn_short[:, :]],
    exper; 
    simulations=[simulation_mother, simulation_mother_sel],
    idxs = [[1,4,11], [1,4,11]],
    acolors = [3,4],
    labels=["No selection", "Selection"],
    mlabel="Mother machine",
    note="mother",
    mother=true)

# Plot the selection function
size_inches = (3, 2.5)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize=18)
xs = 0:1:150
xs = collect(zip(ones(length(xs)), xs))
sx = AgentBasedCells.gen_division_rate_function(hill, rn)
ax = Axis(fig[1,1], xlabel="Protein counts (x)", ylabel="Selection s(x)")

lines!(ax, getindex.(xs, 2), sx.(xs, Ref(params), Ref(0.0)); color=colors[1], linewidth=3)
lines!(ax, getindex.(xs, 2), ones(length(xs)); color=colors[2], linewidth=3)

xlims!(ax, (0, 150))
ylims!(ax, (0, 1.1))

hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
hidespines!(ax, :r, :t)

CairoMakie.activate!()
mkpath("$(plotsdir())/Exp_telegraph")
save("$(plotsdir())/Exp_telegraph/selection.pdf", fig)

