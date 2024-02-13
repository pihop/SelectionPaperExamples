using AdvancedMH
using MCMCChains

function sample_parameters(n)
    f = L * Protein^n / (K^n + Protein^n) 
    model_prod = CellPopulationModel(rn, DivisionRateProduct(γτ₋Sm, f, rn), BinomialKernel(0.5))

    insupport(θ) = 0.0 < θ[1] < 2.0 && 1.0 < θ[2] < 100.0
    f(p) = insupport(p) ? likelihood_combine([ps[31][1], p...], model=model_prod, df_data=df_div₊Sm, interdiv_times=interdiv_times₊Sm) : -Inf

    density_model = DensityModel(f)
    spl = RWMH(MvNormal([0.0, 0.0], [0.1, 10.0]))

    # Sample from the posterior.
    chain = sample(density_model, spl, 1000; init_params=[1.0, 20.0], param_names=["K", "L"], chain_type=Chains)
    return chain
end


