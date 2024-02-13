using GLMakie
GLMakie.activate!()

#results = run_analytical(simulation_parameters)
#model = results.model
#solver = results.solver
#
#f = Figure(resolution = (1200, 600))
#ax1 = Axis(f[1, 1])
#ax2 = Axis(f[1, 2])
#xs = 1:results.solver.maxiters
#lines!(ax1, xs, results.convergence_monitor.growth_factor)
#
#for (i, dist) in enumerate(results.convergence_monitor.birth_dists)
#    lines!(ax2, 1:100, dist; color = (:slategray, 0.4)) 
#end


#cme = model.cme_model(results.birth_dist, model.tspan, model.parameters, solver.truncation)
#cme_solution = solve(cme, solver.solver; cb=PositiveDomain() )# isoutofdomain=(u,p,t) -> any(x -> x .< 0, u), atol=1e-8) 
#
#exp_(τ) = exp.(quadgk(s -> division_rate.(collect.(axes(model.xinit))[1] .- 1, s, model.parameters), 0.0, τ)[1])
#
#f = Figure()
#ax = Axis(f[1, 1])
#sl_t = Slider(f[2, 1], range = 0:0.0001:0.1, startvalue = 0.0)
#
#xs = 1:100
#
#ys = lift(sl_t.value) do x
#    cme_solution(x)
#end
##ys_sum = lift(ys) do x
##    sum(x)
##end
#
#ys2 = lift(sl_t.value) do x
#    cme_solution(x) .* exp_(x)
#end
##
##ys2_sum = lift(ys2) do x
##    sum(x)
##end
##
#barplot!(1:100, ys; )
#barplot!(1:100, ys2;)
#
##text!(ax, "$(ys_sum)", position=(0,0))
##Legend(f[3,1], ax, orientation=:horizontal)
#f
