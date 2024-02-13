using AgentBasedCells
using OrdinaryDiffEq
solver_opts = (stol=1e-3, atol=1e-6, rtol=1e-6, method=Reinsert(), solver=TRBDF2(autodiff=false))
using Distributions
using Interpolations

using PythonCall
@py import skopt: gp_minimize
@py import skopt: callbacks
@py import skopt.callbacks: CheckpointSaver
@py import skopt
@py import dragonfly: minimise_function

using LinearAlgebra

include("utils.jl")
include("models.jl")

bopt_params = (max_evals=15, restarts=1, tstep=5, strain=strain)
pbounds_itp = fill((-10.0, 5.0), nknots)

# BURSTY MODEL.
model_bursty₋Sm = CellPopulationModel(bursty_rn, DivisionRateBounded(γτ₋Sm,1.0,bursty_rn), BinomialKernel(0.5))
# Bayesian optimisation.
bopt_bursty₋Sm = optimise_parameters(model_bursty₋Sm; 
        ps=[], 
        pbounds=[(0.0, 1.0), (1.0, 10.0)], 
        likelihood=likelihood_joint, 
        data=df_div₋Sm, 
        name="notreatment",
        solver_opts=solver_opts,
        exp_setup=experiment_setup,
        bopt_params...
       )
bopt_bursty₋Sm = skopt.load("data/sim/checkpoint_$(strain)_notreatment.pkl")
params_bopt_bursty₋Sm = pyconvert.(Float64, bopt_bursty₋Sm.x)

# SELECTION. HILL FUNCTION.
model_bursty_mult_hill = CellPopulationModel(bursty_rn, DivisionRateBounded(γτ₋Sm*hillfs,1.0,bursty_rn), BinomialKernel(0.5))
bopt₊Sm_bursty_mult_hill = optimise_parameters(model_bursty_mult_hill; 
    ps=params_bopt_bursty₋Sm, 
    pbounds=[(0.0, 5.0), (0.0, 150.0), (0.0, 1.0), (0, 20),],
    likelihood=likelihood_joint, 
    data=df_div₊Sm, 
    name="hillsel",
    solver_opts=solver_opts,
    exp_setup=experiment_setup,
    bopt_params...
    )
bopt_bursty₊Sm_mult_hill = skopt.load("data/sim/checkpoint_$(strain)_hillsel.pkl")
params_bopt_bursty₊Sm_mult_hill = pyconvert.(Float64, bopt_bursty₊Sm_mult_hill.x)

opt₊Sm_bursty_mult_hill_adapt = optimise_parameters(model_bursty_mult_hill; 
    ps=[], 
    pbounds=[(0.0, 1.0), (1.0, 5.0), (0.0, 5.0), (0.0, 150.0), (0.0, 1.0), (0, 20), ],
    likelihood=likelihood_joint, 
    data=df_div₊Sm, 
    name="hilladapt",
    solver_opts=solver_opts,
    exp_setup=experiment_setup,
    bopt_params...
   )
bopt_bursty₊Sm_mult_adapt = skopt.load("data/sim/checkpoint_$(strain)_hilladapt.pkl")
params_bopt_bursty₊Sm_adapt = pyconvert.(Float64, bopt_bursty₊Sm_mult_adapt.x)
