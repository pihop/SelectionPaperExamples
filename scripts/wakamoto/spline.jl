include("fitting.jl")

@named bursty_rn_itp = ReactionSystem(rxs, t, [Protein], [k, b, y1, y2, y3])

# Set up the spline.
knots = [10, 50, 90, 130]
#knots = [50, 90, 130]
itp_(k1,k2,k3,k4,v1,v2,v3,v4,x) = LinearInterpolation(knots, [y1, y2, y3, y4]; extrapolation_bc=Line())(x)
@register itp_(y1, y2, y3, y4, Protein)
f_itp = itp_(y1, y2, y3, y4, Protein) 

model_mult_itp = CellPopulationModel(bursty_rn_itp, DivisionRate(γτ₋Sm + f_itp, bursty_rn_itp), BinomialKernel(0.5))
opt₊Sm_itp_mult_hill = optimise_parameters_itp(
    model_mult_itp; 
    ps=ps_bursty₋Sm, 
#    pbounds=[[0,100],[0, 15],[0, 15], [0, 15]],
    pbounds=[[0,1],[0, 1],[0, 1],[0,1]],
    divrate=γτ₋Sm, 
    likelihood=likelihood_joint, 
    data=df_div₊Sm, 
    inf=-60000, 
    max_evals=150,
    trn=working_trn)
# LinearInterpolation
# [0.44444162, 1.7571621, 7.62078951e-04, 2.54026317e-04, 1.48376772e+00, 1.96387746e+0] 34209.4999
# [0.25842097, 2.5785703, 9.99997177e+00, 2.82251463e-05, 1.11113934e+00, 1.35805292e+0] 33670.537
# [0.25842097, 2.5785703, 1.49988569e+01, 3.81039476e-04, 9.84224966e-01, 2.04046639e+0] 33696.7219
#
# [0.25842097, 2.57857034, 3.71056241, 0.36351166, 1.11339735, 7.7709190] 34405.0846

#fit_param = [0.25842097, 2.5785703, 9.99997177e+00, 2.82251463e-05, 1.11113934e+00, 1.35805292e+0]
#
#analyticals_bursty_mult₊Sm = run_analytical_single(model_mult_itp, experiment_setup(model_parameters=fit_param, trn=working_trn);)
#division_dist!(analyticals_bursty_mult₊Sm)
#
#fig = Figure(resolution=(1000, 800))
#ax1 = Axis(fig[1,1], xlabel="Protein counts", ylabel="Probability density", title="Birth distribution")
#ax2 = Axis(fig[2,1], xlabel="Protein counts", ylabel="Probability density", title="Division distribution")
#ax3 = Axis(fig[3,1], xlabel="Time", ylabel="Probability density", title="Interdivision time")
#ax4 = Axis(fig[4,1], xlabel="Protein counts x", ylabel="Selection f(x)")
#
#xs = collect(1:working_trn)
#hist!(ax1, df_birth₋Sm[:, :Column1]; normalization=:pdf, bins=100, color=Cycled(1), label="No treatment")
#lines!(ax1, xs, analyticals_bursty₋Sm.results[:birth_dist]; color=Cycled(1), linewidth=3.0)
#
#hist!(ax1, df_birth₊Sm[:, :Column1]; normalization=:pdf, bins=100, color=Cycled(2), label="Treatment")
#lines!(ax1, xs, analyticals_bursty_mult₊Sm.results[:birth_dist]; color=Cycled(2), linewidth=3.0)
##xlims!(ax1, (0, working_trn))
#xlims!(ax1, (0, 200))
#
#hist!(ax2, df_div₋Sm[:, :Column2]; normalization=:pdf, bins=100, color=Cycled(1))
#lines!(ax2, xs, analyticals_bursty₋Sm.results[:division_dist]; color=Cycled(1), linewidth=3.0)
#
#hist!(ax2, df_div₊Sm[:, :Column2]; normalization=:pdf, bins=100, color=Cycled(2))
#lines!(ax2, xs, analyticals_bursty_mult₊Sm.results[:division_dist]; color=Cycled(2), linewidth=3.0)
#xlims!(ax2, (0, 200))
#
#ts = 0.0:0.1:150.0
#div_dist₋Sm = division_time_dist(analyticals_bursty₋Sm)
#div_dist₊Sm = division_time_dist(analyticals_bursty_mult₊Sm)
#hist!(ax3, interdiv_times₋Sm .* 5; bins=25, normalization=:pdf, color=Cycled(1), label="Data")
#hist!(ax3, interdiv_times₊Sm .* 5; bins=25, normalization=:pdf, color=Cycled(2), label="Data")
#lines!(ax3, ts, div_dist₋Sm.(ts); color=Cycled(1))
#lines!(ax3, ts, div_dist₊Sm.(ts); color=Cycled(2))
#
##lines!(ax2, ts, div_dist_hist.(ts))
##lines!(ax3, xs, fx_.(xs; ps=fit_param[2:end], n=2))
##lines!(ax3, xs, itp_fit.(xs))
#lines!(ax4, xs, itp_.(fit_param[3:end]..., xs))
#xlims!(ax4, (0, 200))
#
#Legend(fig[5, 1], ax1, orientation = :horizontal, tellwidth = false, tellheight = true)
##mkpath("$(plotsdir())/Wakamoto")
##
##save("$(plotsdir())/Wakamoto/fitting_mcmc_$strain.pdf", fig)
#
