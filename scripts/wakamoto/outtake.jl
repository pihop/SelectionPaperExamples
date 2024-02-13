#analyticals₋Sm = run_analytical_single(model₋Sm, experiment_setup(model_parameters=[ps[62]...], trn=150);)
#analyticals₊Sm = run_analytical_single(model₊Sm, experiment_setup(model_parameters=[ps[38]...], trn=150);)
#fl(x; L, K) = L * x^10 / (K^10 + x^10)
#
#function G(ps; γ, sol) 
#    L,K = ps
#    return t -> γ(t) / (sol.cmesol(t)[end] .* sum(fl.(collect(1:150); L=L, K=K) .* sol.cmesol(t)[1:150]))
#end

N = 50
Ls = range(0.4, stop=2.0, length=N) 
Ks = range(1.0, stop=90.0, step=10.0)
ps₋ = [[ps[62][1],L,K] for (L,K) in Iterators.product(Ls, Ks)]
ps₊ = [[ps[38][1],L,K] for (L,K) in Iterators.product(Ls, Ks)]

function likelihoodg(p; γ, sol₋Sm, data)
    g = G(ps₋[1,1]; γ=γ, sol=sol₋Sm)
    itpspan = range(0.0, stop=150.0, length=100)
    A = g.(itpspan)
    gi = linear_interpolation(A)
    Gg(t) = gi(t)

 #   model_prod₊Sm = CellPopulationModel(rn, DivisionRateProduct(γτ₋Sm, Gg, rn), BinomialKernel(0.5))
end


## From above we have protein production rates for both cases. 
## Suppose that γ(x,τ) = f(x)*g₋Sm(τ)
#f = L * Protein^4 / (K^4 + Protein^4)
#model_prod₋Sm = CellPopulationModel(rn, DivisionRateProduct(γτ₋Sm, f, rn), BinomialKernel(0.5))
#@register Gg(t)
#model_prod₊Sm = CellPopulationModel(rn, DivisionRateProduct(γτ₋Sm, Gg, rn), BinomialKernel(0.5))

#ls_prod₊Sm = likelihood.(ps₊; model=model_prod₊Sm, λ=nothing, df_data=df_birth₊Sm)
#ls_prod₋Sm = likelihood.(ps₋; model=model_prod₋Sm, λ=nothing, df_data=df_birth₋Sm)


#N = 20
#Ks = range(0.4, stop=1.5, length=N) 
#Ls = range(50.0, stop=90.0, step=10.0)
#ps_prod₊Sm = [[ps[idx₊Sm][1], L, K] for (L, K) in Iterators.product(Ks, Ls)]
#ps_prod₋Sm = [[ps[idx₋Sm][1], L, K] for (L, K) in Iterators.product(Ks, Ls)]
##ls_prod₊Sm = likelihood.(ps_prod₊Sm; model=model_prod₊Sm, λ=nothing, df_data=df_birth₊Sm)
##ls_prod₋Sm = likelihood.(ps_prod₋Sm; model=model_prod₋Sm, λ=nothing, df_data=df_birth₋Sm)
#max_prod₊Sm, idx_prod₊Sm = findmin(ls_prod₊Sm)
#max_prod₋Sm, idx_prod₋Sm = findmin(ls_prod₋Sm)
#
#idx_opt = findmin(ls_prod₊Sm .+ ls_prod₋Sm)
# Bursty.
#model_bursty = CellPopulationModel(bursty_rn, DivisionRate(γτ, bursty_rn), BinomialKernel(0.5))
#bs = range(1, stop=2, length=20)
#cs = range(0.4, stop=0.8, length=20)
#ps_bursty = [[c/b, 0, 0, b] for (b,c) in Iterators.product(bs, cs)]
#ls_bursty = likelihood.(ps_bursty; model=model_bursty, λ=λ, data=df_birth₋Sm)
#max_bursty, idx_bursty = findmin(ls_bursty)
#analyticals_bursty = run_analytical_single(model_bursty, experiment_setup(model_parameters=[ps_bursty[idx_bursty]...], trn=150);)

## Setting 2, suppose γ(x,τ) = f(x)g(τ)
## We have already fitted the marginal γ(τ) = g(τ)E[f(x)].
## Assume f(x) =  l(xⁿ / (K + xⁿ)).
#f = L * Protein^10 / (K^2 + Protein^10)
#model_prod = CellPopulationModel(rn, DivisionRateProduct(γτ, f, rn), BinomialKernel(0.5))
#model_bursty_prod = CellPopulationModel(bursty_rn, DivisionRateProduct(γτ, f, bursty_rn), BinomialKernel(0.5))
#
## Non-bursty.
#N = 10
#Ks = range(10.0, stop=100.0, length=N) 
#Ls = range(2.0, stop=10.0, step=1.0)
#ps_prod_ = [[ps[idx][1], L, K] for (L, K) in Iterators.product(Ks, Ls)]
##ls_prod = likelihood.(ps_prod; model=model_prod, λ=λ, data=df_birth₋Sm)
##max_prod, idx_prod = findmin(ls_prod)
##analyticals_prod = run_analytical_single(model_prod, experiment_setup(model_parameters=[ps_prod[idx_prod]...], trn=150);)
## Bursty
#
#ps_bursty_prod = [[ps_bursty[idx_bursty][1], L, K, ps_bursty[idx_bursty][4]] for (L, K) in Iterators.product(Ks, Ls)]
#ls_bursty_prod = likelihood.(ps_bursty_prod; model=model_bursty_prod, λ=λ)


# Compare
#analyticals₊Sm = run_analytical_single(model₊Sm, experiment_setup(model_parameters=[ps[idx₊Sm]...], trn=150);)
#analyticals₋Sm = run_analytical_single(model₋Sm, experiment_setup(model_parameters=[ps[idx₋Sm]...], trn=150);)
#analyticals_prod₊Sm = run_analytical_single(model_prod₊Sm, experiment_setup(model_parameters=[ps_prod₊Sm[7, 1]...], trn=150);)
#analyticals_prod₋Sm = run_analytical_single(model_prod₋Sm, experiment_setup(model_parameters=[ps_prod₋Sm[7, 1]...], trn=150);)
#
# Plots
fig = Figure()
#ax1 = Axis(fig[1,1], xlabel="Protein counts", title="No treatment Sm₋")
#ax2 = Axis(fig[2,1], xlabel="Protein counts", title="Treatment Sm₊")
###contourf!(ax, vec(Point2f.(getindex.(ps, 1), getindex.(ps, 4))), color=ls)
###contourf!(ax, bs, cs, ls) 
#lines!(ax2, collect(1:150), analyticals₊Sm.results[:birth_dist]; color=Cycled(2), linewidth=2.0)
#lines!(ax2, collect(1:150), analyticals_prod₊Sm.results[:birth_dist]; color=Cycled(3), linewidth=2.0)
###barplot!(ax, collect(1:150), analyticals2.results[:birth_dist])
#hist!(ax2, df_birth₊Sm[:, :Column1]; normalization=:pdf, bins=100, color=Cycled(1))
#xlims!(ax2, (0, 150))
#
#lines!(ax1, collect(1:150), analyticals₋Sm.results[:birth_dist]; color=Cycled(2), linewidth=2.0)
#lines!(ax1, collect(1:150), analyticals_prod₋Sm.results[:birth_dist]; color=Cycled(3), linewidth=2.0)
###barplot!(ax, collect(1:150), analyticals2.results[:birth_dist])
#hist!(ax1, df_birth₋Sm[:, :Column1]; normalization=:pdf, bins=100, color=Cycled(1))
#xlims!(ax1, (0, 150))
#
# Plots pareto
#ax3 = Axis(fig[3,1])
#scatter!(ax3, vec(ls_prod₊Sm), vec(ls_prod₋Sm))
#scatter!(ax3, ls_prod₊Sm[7,1], ls_prod₋Sm[7,1])
#scatter!(ax3, ls_prod₊Sm[19,1], ls_prod₋Sm[19,1])
#scatter!(ax3, -ls_prod₊Sm[12,1], -ls_prod₋Sm[12,1])

#sol₋ = analyticals₋Sm.cmesol
#sol₊ = analyticals₊Sm.cmesol
ax4 = Axis(fig[4,1])

#pl₋(u; L, K) = g₋Sm(u)*exp(sol₋(u)[end])sum(fl.(collect(1:150); L=L, K=K) .* sol₋(u)[1:end-1])
#pl₊(u; L, K) = g₋Sm(u)*exp(sol₊(u)[end])sum(fl.(collect(1:150); L=L, K=K) .* sol₊(u)[1:end-1])

tes(x;L) = L*x*exp(L*x^2)

lines!(ax4, 0.0:0.001:150.0, g.(0.0:0.001:150.0))
lines!(ax4, 0.0:0.001:150.0, tes.(0.0:0.001:150.0; L=0.001))
#lines!(ax4, 0.0:0.001:100.0, g₊Sm.(0.0:0.001:100.0))
###lines!(ax4, 0.0:0.001:150.0, pl.(0.0:0.001:150.0; L=1.0, K=50))
###lines!(ax4, 0.0:0.001:150.0, exp.(getindex.(sol(0.0:0.001:150.0).u, 151)))
#lines!(ax4, 0.0:0.001:100.0, pl₋.(0.0:0.001:100.0; L=0.733623, K=61.0))
#lines!(ax4, 0.0:0.001:100.0, pl₊.(0.0:0.001:100.0; L=0.733623, K=61.0))

