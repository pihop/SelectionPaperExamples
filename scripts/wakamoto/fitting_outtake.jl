#=
#model_bursty_mult = CellPopulationModel(
#    bursty_rn_itp, 
#    DivisionRateBounded(γτ₋Sm*itpp_spline₊Sm(Protein, [v...][1:nknots]...), 1.0, bursty_rn_itp), 
#    BinomialKernel(0.5))

#bopt₊Sm_bursty_mult_itp = optimise_parameters(model_bursty_mult; 
#    ps=params_bopt_bursty₋Sm, 
#    pbounds=pbounds_itp,
#    likelihood=likelihood_joint, 
#    data=df_div₊Sm, 
#    name="itpselection",
#    solver_opts=solver_opts,
#    exp_setup=experiment_setup,
#    bopt_params...
#    )
#bopt₊Sm_bursty_mult_itp = skopt.load("data/sim/checkpoint_$(strain)_itpselection.pkl")

#model_bursty₊Sm_itp = CellPopulationModel(
#    bursty_rn_itp, 
#    DivisionRateBounded(γτ₊Sm*itpp_spline₊Sm(Protein, [v...][1:nknots]...) ,1.0,bursty_rn_itp),
#    BinomialKernel(0.5))

#bopt_bursty₊Sm = optimise_parameters(model_bursty₊Sm_itp; 
#        ps=[], 
#        pbounds=[(0.0, 1.0), (1.0, 10.0), pbounds_itp...], 
#        likelihood=likelihood_joint, 
#        data=df_div₊Sm, 
#        trn=working_trn, 
#        name="treatment",
#        solver_opts=solver_opts,
#        exp_setup=experiment_setup,
#        bopt_params...)


























# Reparam reaction network.
# [1.3155011636750482, 72.17590012630849, 0.13581509283555265, 28]

#opt₊Sm_bursty_mult_hill_sel = optimise_parameters(model_bursty_mult; 
#    ps=[], 
#    pbounds=[(0.0, 1.0), (1.0, 5.0), (0.0, 5.0), (0.0, 150.0), (0.0, 1.0), (0, 20), ],
#    likelihood=likelihood_joint, 
#    data=df_div₊Sm, 
#    n_initial_points=10,
#    trn=working_trn,
#    name="hilladapt",
#    strain=strain,
#    tstep=tslice,
#    restarts=3,
##    x0= [0.6976711236526892, 2.064277138400638, 1.3771485956602671, 78.24118123524796, 0.1039282608597518, 10],
#    max_evals=250,
#    solver_opts=solver_opts,
#    exp_setup=experiment_setup)
# [0.6976711236526892, 2.064277138400638, 1.3771485956602671, 78.24118123524796, 0.1039282608597518, 10] 
# 38286.781490473775

model_bursty_mult = CellPopulationModel(bursty_rn, DivisionRateBounded(γτ₋Sm*hillfs*hillfr,1.0,bursty_rn), BinomialKernel(0.5))
#opt₊Sm_bursty_mult_hill_sel = optimise_parameters(model_bursty_mult; 
#    ps=[],
#    pbounds=[(0.0, 1.0), (1.0, 5.0), (0.0, 5.0), (0.0, 150.0), (0.0, 1.0), (0, 20), (0.0, 1.0), (0.0, 250.0), (0, 50)],
#    likelihood=likelihood_joint, 
#    data=df_div₊Sm, 
#    trn=working_trn,
#    name="hilladaptrs",
#    max_evals=1000,
#    tstep=5,
#    solver_opts=solver_opts)

# [0.7095197188938238, 2.138316409463015, 1.3228327110171074, 107.17366474408334, 0.3841850330172725, 5, 1.0, 90.30382933314671, 6]

model_bursty₊Sm = CellPopulationModel(bursty_rn, DivisionRateBounded(γτ₊Sm,1.0,bursty_rn), BinomialKernel(0.5))
#opt_bursty₊Sm = optimise_parameters(model_bursty₊Sm; 
#        ps=[], 
#        pbounds=[(0.0, 1.0), (1.0, 10.0)], 
#        likelihood=likelihood_joint, 
#        data=df_div₊Sm, 
#        trn=working_trn, 
#        name="treatmentbase",
#        max_evals=1000,
#        tstep=5,
#        solver_opts=solver_opts)
# [0.7138359186007392, 1.792217506859519]  36186.140036441466

# Spline
#model_spline₊Sm = CellPopulationModel(bursty_rn_base, DivisionRateBounded(γτ₋Sm*itp₊Sm_(Protein),1.0,bursty_rn_base), BinomialKernel(0.5))
#opt_spline₊Sm = optimise_parameters(model_spline₊Sm; 
#        ps=[], 
#        pbounds=[(0.0, 1.0), (1.0, 10.0)], 
#        likelihood=likelihood_joint, 
#        data=df_div₊Sm, 
#        trn=working_trn, 
#        name="treatmentbase",
#        strain=strain,
#        max_evals=200,
#        x0=[0.7138359186007392, 1.792217506859519],
#        tstep=5,
#        restarts=3,
#        solver_opts=solver_opts,
#        exp_setup=experiment_setup)
# [0.6467306826927604, 1.5010777373338273]  38978.37344539174
=#
