using AgentBasedFSP
using OrdinaryDiffEq
solver_opts = (stol=1e-4, atol=1e-7, rtol=1e-7, method=Reinsert(), solver=TRBDF2(autodiff=false))
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

bopt_params = (max_evals=200, restarts=5, tstep=tslice, strain="sos")

# No simulation so bound not important.
model_wt = MotherCellModel(bursty_rn, DivisionRateBounded(γτ_wt, 1.0, bursty_rn), BinomialKernel(0.5))

# Bayesian optimisation.
opt_wt = optimise_parameters(model_wt; 
        ps=[],
        pbounds=[(0.0, 1.0), (1.0, 10.0)], 
        likelihood=likelihood_joint, 
        data=df_div_wt, 
        name="wt",
        exp_setup=experiment_setup,
        solver_opts=solver_opts,
        bopt_params...)
bopt_wt = skopt.load("data/sim/checkpoint_sos_wt.pkl")
params_wt = pyconvert.(Float64, bopt_wt.x)

model_2pal = MotherCellModel(bursty_rn, DivisionRateBounded(γτ_wt*hillfs, 1.0, bursty_rn), BinomialKernel(0.5))
opt_2pal_sel = optimise_parameters(model_2pal; 
    ps=params_wt, 
    pbounds=[(0.0,5.0), (0.0, 100.0), (0.0, 1.0), (0, 20)],
    likelihood=likelihood_joint, 
    data=df_div_2pal, 
    name="2palsel",
    trn=working_trn,
   exp_setup=experiment_setup,
    solver_opts=solver_opts,
    bopt_params...
    )
bopt_2palsel = skopt.load("data/sim/checkpoint_sos_2palsel.pkl")
params_2palsel = pyconvert.(Float64, bopt_2palsel.x)

model_nosel = MotherCellModel(bursty_rn, DivisionRateBounded(γτ_2pal, 1.0, bursty_rn), BinomialKernel(0.5))
opt_wt = optimise_parameters(model_wt; 
        ps=params_wt,
        pbounds=[(0.0, 1.0), (1.0, 5.0), (0, 50), (0.0, 100.0)], 
        likelihood=likelihood_joint, 
        data=df_div_2pal, 
        trn=working_trn, 
        name="growth_sos_2palnosel",
        exp_setup=experiment_setup,
        solver_opts=solver_opts,
        bopt_params...)

#ps_kde = [0.2894619313035903, 1.8880591761464778, 36, 27.80838442006454] # 6183.5850998401465
# without degredation [0.3031134152596094, 1.464297829619517] 6259.170875106497
#ps_wt_kde = [0.3031134152596094, 1.464297829619517]
#
model_2pal = MotherCellModel(bursty_rn, DivisionRateBounded(γτ_wt*hillfs0*hillfs, 1.0, bursty_rn), BinomialKernel(0.5))
opt_2pal_hill_wtparams = optimise_parameters(model_2pal; 
    ps=ps_kde, 
    pbounds=[(0.0,5.0), (0.0, 100.0), (0.0, 1.0), (0, 20)],
    likelihood=likelihood_joint, 
    data=df_div_2pal, 
    name="growth_sos_2pal",
    inf=-70000, 
    trn=working_trn,
    max_evals=1000)
# [0.9968777895952129, 43.93932683765189, 0.0, 11] 30504.558142055022
# [1.4721004302953467, 53.97132637287041, 0.0, 4] 30534.98643063118
# [1.0469129929806757, 44.678292951350954, 0.0, 11] 30369.147847181346

opt_2pal_hill_2palparams = optimise_parameters_hill(model_2pal; 
    ps=ps_2pal_kde, 
    pbounds=[(0.0, 1.0), (1.0, 20.0), (0.0,5.0), (0.0, 100.0), (0.0, 1.0), (0, 20)],
    likelihood=likelihood_joint, 
    data=df_div_2pal, 
    inf=-70000, 
    name="growth_sos_2pal_",
    trn=working_trn,
    max_evals=1000)

#model_2pal_spl = MotherCellModel(bursty_rn_react, DivisionRateBounded(γτ_wt*itp2palwt_(Protein), 1.0, bursty_rn_react), BinomialKernel(0.5))
#opt_2pal_itp = optimise_parameters(model_2pal_spl; 
#    ps=[], 
#    x0=[0.3031134152596094, 1.464297829619517],
#    pbounds=[(0.0, 1.0), (1.0, 20.0)],
#    likelihood=likelihood_joint, 
#    data=df_div_2pal, 
#    inf=-70000, 
#    name="2pal_spl",
#    trn=working_trn,
#    strain="sos",
#    tstep=tslice,
#    max_evals=1000,
#    restarts=3,
#    exp_setup=experiment_setup,
#    solver_opts=solver_opts)
# [0.3026266852208091, 3.4554597340429414] 30250.067426279616
