using Interpolations
import Distributions: pdf, logpdf, cdf, logcdf, logccdf
using Roots
using StatsBase
using LinearAlgebra
using ModelingToolkit
using AdvancedMH
using Serialization
using MCMCChains
using LogDensityProblems
using Turing
using LoopVectorization
using QuadGK
using Zygote
using ForwardDiff
using ReverseDiff

Turing.setadbackend(:forwarddiff)

using AgentBasedCellss
#using Bootstrap
#using OnlineStats

using PythonCall
@py import skopt: gp_minimize
@py import skopt: callbacks
@py import skopt.callbacks: CheckpointSaver
@py import skopt

mutable struct OnlineMean
    value::Real
    n::Int64
    function OnlineMean(value::Real, n::Int64)
        return new(value, n)
    end
end

function OnlineMean()
    return OnlineMean(0.0, 0) 
end

function update_mean!(mean::OnlineMean, x)
    if mean.n != 0
        n = mean.n
        μ = mean.value
        mean.value = n / (n + 1) * (x / n + μ)
        mean.n += 1
    else
        mean.n = 1
        mean.value = x
    end
end

function value(mean::OnlineMean)
    return mean.value
end

struct LogTargetDensity
    loglik
    dim::Int
end

LogDensityProblems.logdensity(p::LogTargetDensity, θ) = p.loglik(θ) 
LogDensityProblems.dimension(p::LogTargetDensity) = p.dim
LogDensityProblems.capabilities(::Type{LogTargetDensity}) = LogDensityProblems.LogDensityOrder{0}()

InterpKDELin(kde::UnivariateKDE) = InterpKDE(kde, BSpline(Linear()))

struct CumulativeKDE
    interp
    hazard
    kde
    age
    λ
    function CumulativeKDE(cdf_interp, hazard_interp, kde, age, λ)
        new(cdf_interp, hazard_interp, kde, age, λ)
    end
end

function CumulativeKDEPop(kde; rtol=1e-10, atol=1e-10, tend)
    # Compute the growth rate and ancestral distribution.
    f(λ) = 2 - quadgk(t -> exp(λ*t) * pdf(kde,t), 0, tend)[1]
    λ0 = find_zero(f, [0.0, 1.0]; )
    display(λ0)

    A = [quadgk(s -> 0.5*exp(λ0*s)*pdf(kde, s), 0, x)[1] for x in kde.kde.x]
#    A = [quadgk(s -> pdf(kde, s), 0, x)[1] for x in kde.kde.x]

    cdf_interp = Interpolations.scale(interpolate(A, BSpline(Linear())), kde.kde.x)
    cdf_interp = extrapolate(cdf_interp, Interpolations.Flat())

    H = [exp(log(0.5*exp(λ0*x)*pdf(kde,x)) - log(max(1e-12, 1-cdf_interp(x)))) for x in kde.kde.x]
#    H = map(x -> isinf(x) ? 10000.0 : x, H) 
#    H = [exp(log(pdf(kde,x)) - log(max(1e-16, 1-cdf_interp(x)))) for x in kde.kde.x]

    hazard_interp = Interpolations.scale(interpolate(H, BSpline(Linear())), kde.kde.x)
    hazard_interp = extrapolate(hazard_interp, Interpolations.Flat())

#pd₋Sm(t) = 2*exp(-λ₋Sm*t)*marginalγ₋Sm(t) * exp(-integ₋Sm(t))
#    pdfhaz(τ) = hazard_interp(τ)*exp(-quadgk(s -> hazard_interp(s), 0, τ)[1])
#    norm = quadgk(s -> pdfhaz(s), 0.0, tend)
#    display(norm)
#    age(τ) = exp(-quadgk(s -> hazard_interp(s), 0, τ)[1])
#    normal = quadgk(s -> age(s), 0, tend, rtol=rtol)[1]
#    display(normal)

#    Aage = age.(kde.kde.x) #./ normal
#    age_interp = Interpolations.scale(interpolate(Aage, BSpline(Linear())), kde.kde.x)
#    age_interp = extrapolate(age_interp, Interpolations.Flat())

    CumulativeKDE(cdf_interp, hazard_interp, kde.kde, t -> 0.0, λ0)
end

function CumulativeKDE(kde; rtol=1e-8, tend, λ=0.0)
    A = [quadgk(s -> pdf(kde, s), 0, x, rtol=rtol)[1] for x in kde.kde.x]
#    A[end] = 1.0

    cdf_interp = Interpolations.scale(interpolate(A, BSpline(Linear())), kde.kde.x)
    cdf_interp = extrapolate(cdf_interp, Interpolations.Flat())

    H = [exp(log(pdf(kde,x)) - log(max(1e-12, 1-cdf_interp(x)))) for x in kde.kde.x]

    hazard_interp = Interpolations.scale(interpolate(H, BSpline(Linear())), kde.kde.x)
    hazard_interp = extrapolate(hazard_interp, Interpolations.Flat())

    CumulativeKDE(cdf_interp, hazard_interp, kde.kde, t -> 0.0, λ)
end

Distributions.pdf(kde::CumulativeKDE, x) = pdf(kde.kde, x)
Distributions.logpdf(kde::CumulativeKDE, x) = log(pdf(kde.kde, x))
Distributions.cdf(kde::CumulativeKDE, x) = kde.interp(x)
Distributions.logcdf(kde::CumulativeKDE, x) = log(kde.interp(x))
Distributions.logccdf(kde::CumulativeKDE, x) = log(1-kde.interp(x))
hazard(kde::CumulativeKDE, x) = kde.hazard(x)

function likelihood_joint(θ; model, data, inf, tstep, solver_opts, exp_setup, kwargs...)
    # Run analytical.
    display(θ)
    analyticals = run_analytical_single(model, exp_setup(model_parameters=θ); solver_opts...)   
    if analyticals.flag == :Success
        joint_fpt_ancest_cdist!(analyticals; step=tstep) # Ancestral distribution.
        joint_cdist = analyticals.results[:joint_fpt_ancest_cdist]

        liks = [joint_cdist[Int64(row[:Column1])][Int(round(row[:Column3]))+1] for row in eachrow(data)]
        lik = sum(log.(liks .+ 1e-4))
        return lik 
    else 
        return inf 
    end
end

function optimise_parameters(model; 
    x0=nothing, 
    ps, 
    pbounds, 
    likelihood, 
    data, 
    inf=-100000, 
#    trn, 
    max_evals, 
    n_initial_points=10, 
    name="", 
    strain, 
    restarts, kwargs...)

    if !isnothing(x0)
        x0_ = pylist(x0)  
    else
        x0_ = pybuiltins.None
    end

    res = nothing

    for n in 1:restarts
        println("Epoch $n out of $restarts")
        checkpoint_saver = CheckpointSaver("data/sim/checkpoint_$(strain)_$name.pkl", store_objective=false)
        len = length(parameters(model.molecular_model)) - length(pbounds) - length(ps)
        #display((ps..., pyconvert.(Float64, [1.0, 1.0, 1.0, 1.0, 1.0])..., zeros(len)...))
        lik(p) = -likelihood((ps..., pyconvert.(Float64, p)..., zeros(len)...); model=model, data=data, inf=inf, kwargs...)
        res = gp_minimize(lik, pbounds, verbose=true, 
            n_calls=max_evals, n_initial_points=n_initial_points, 
            x0=x0_,
            initial_point_generator="lhs", callback=checkpoint_saver, n_jobs=Base.Threads.nthreads())

        x0_ = res.x
        val = pyconvert(Float64, res.fun)
        println("Best value $val")
    end
    return res
end

function optimise_parameters_ctn(model; ps, pbounds, likelihood, data, inf=-100000, trn, max_evals, n_initial_points=10, name="", kwargs...)
    checkpoint_saver = CheckpointSaver("data/sim/checkpoint_$name.pkl", compress=9)
    res = skopt.load("data/sim/checkpoint_$name.pkl")
    len = length(parameters(model.molecular_model)) - length(pbounds) - length(ps)
    lik(p) = -likelihood([ps..., pyconvert.(Float64, p)..., zeros(len)...]; model=model, data=data, inf=inf, trn=trn, kwargs...)
    res = gp_minimize(lik, pbounds, verbose=true, x0 = res.x_iters, y0 = res.func_vals,
        n_calls=max_evals, n_initial_points=n_initial_points, 
        initial_point_generator="lhs", callback=checkpoint_saver, n_jobs=Base.Threads.nthreads())
end

function empirical_marginal_γ(values; ts)
    marginalτ = Dict(i => OnlineMean() for i in ts)
    for val in values
        update_mean!(marginalτ[val[2]], val[1])
    end

    marginalτ[0] = OnlineMean(value(marginalτ[ts[1]]), 0)
    marginalτ[ts[1]+ts[end]] = OnlineMean(value(marginalτ[ts[end]]), 0) 
    return marginalτ
end

function add_marginal_ends(marginalτ, ts)
    marginalτ[0] = OnlineMean(value(marginalτ[ts[1]]), 0)
    marginalτ[ts[1]+ts[end]] = OnlineMean(value(marginalτ[ts[end]]), 0) 
    return marginalτ
end

function empirical_marginal_cycles_γ(γ, θ, cycles; tslice, offset=0.5)
    vals = [] 
    tss = []

    for cyc in cycles
        r = 1:length(cyc)
        ts_ = r .* tslice #.- offset*tslice 
        pairs = zip(cyc, ts_)
        push!(tss, last(ts_))
        int, last_ = trapz_saveval(x -> γ(x[1], θ, x[2]), pairs, tslice, vals;) 
    end

    maxr = 1:maximum(length.(cycles))
    maxts = maxr .* tslice #.- offset*tslice    
    γ = empirical_marginal_γ(vals; ts=maxts)
    γs = [value(γ[t]) for t in maxts]

    marginalγ = extrapolate(
        Interpolations.scale(
#            interpolate(γs, BSpline(Cubic(Line(OnGrid())))),
            interpolate(γs, BSpline(Linear())),
            maxts),
        Interpolations.Flat()
       )

    return marginalγ
end

function trapz(f, as, h::Real)
    int::Real = 0.0 
    as_ = collect(as)
    fa = f(as_[1])
    fb::Real = 0.0

    @inbounds @fastmath for (a, b) in zip(as_, @view as_[2:end])
        fb = f(b)
        int += 0.5 * h * (fa + fb)    
        fa = fb
    end
    return int, (fb, as_[end][2])
end

function trapz_marginal(f, as, h::Real)
    int::Real = 0.0 
    as_ = collect(as)
    fa = f(as_[1])
    fb::Real = 0.0

    @inbounds @fastmath for (a, b) in zip(as_, @view as_[2:end])
        fb = f(b)
        int += 0.5 * h * (fa + fb)    
        fa = fb
    end
    return int
end


function trapz_saveval(f, as, h::Real, save)
    int::Real = 0.0 
    as_ = collect(as)
    fa = f(as_[1])
    fb::Real = 0.0

    @inbounds @fastmath for (a, b) in zip(as_, @view as_[2:end])
        fb = f(b)
        int += 0.5 * h * (fa + fb)    
        push!(save, (fa, a[2]))
        fa = fb
    end

    push!(save, (fb, as_[end][2]))
    return int, fb
end

function trapz_updatemarginal(f, as, h::Real, marginalτ)
    int::Real = 0.0 
    as_ = collect(as)
    fa = f(as_[1])
    fb::Real = 0.0

    @inbounds @fastmath for (a, b) in zip(as_, @view as_[2:end])
        fb = f(b)
        int += 0.5 * h * (fa + fb)    
        update_mean!(marginalτ[a[2]], fa)
        fa = fb
    end

    update_mean!(marginalτ[as_[end][2]], fb)
    return int, (fb, as_[end][2])
end

D(f) = x -> ForwardDiff.derivative(f,float(x))

function likelihood_hazard_pop(γ::Function, θ, q; data, inf, tslice, x0, offset=0.5)
    maxlen = maximum(length.(data))
    maxr = 0:maxlen-1
    maxts = maxr .* tslice 
    marginalτ = Dict(i => OnlineMean() for i in maxts)

    out::Real = 0.0
    tss = []

    for cyc in data
        ts_ = maxts[1:length(cyc)] 
        pairs = zip(cyc, ts_)
        int, last_ = trapz_updatemarginal(x -> γ(x[1], θ, x[2]), pairs, tslice, marginalτ) 
        out += log(last_[1])
        push!(tss, last_[2])
        out -= int 
    end

    γ = add_marginal_ends(marginalτ, maxts)

    times = sort(collect(keys(γ)))
    pdfγ(t) = value(γ[t])*exp(-trapz_marginal(t -> value(γ[t]), times[times .<= t], tslice))
    fx(λ) = 2 - trapz_marginal(t -> exp(ForwardDiff.value(λ)*t) * pdfγ(t), times, tslice)

    λ0 = find_zero(fx, zero(eltype(θ)))
    out -= sum(λ0 .* tss)

    if !isnothing(x0)
        out -= q*norm(x0 - θ)^2
    end

    if !isnan(out) && !isinf(out) 
        return out
    else
        return inf
    end
end

function growth_rate(hazard, θ, paramrn; data, tslice, offset=0.5)
    γ = AgentBasedCellss.gen_division_rate_function(hazard, paramrn)
    
    vals = [] 
    tss = []

    for cyc in data
        r = 0:length(cyc)-1
        ts_ = r .* tslice 
        pairs = zip(cyc, ts_)
        push!(tss, last(ts_))
        int, last_ = trapz_saveval(x -> γ(x[1], θ, x[2]), pairs, tslice, vals;) 
    end

    maxr = 0:maximum(length.(data))-1
    maxts = maxr .* tslice# .- offset*tslice    
    γ = empirical_marginal_γ(vals; ts=maxts)
    times = sort(collect(keys(γ)))

    pdfγ(t) = value(γ[t])*exp(-trapz_marginal(t -> value(γ[t]), times[times .<= t], tslice))
    f(λ) = 2 - trapz_marginal(t -> exp(λ*t) * pdfγ(t), times, tslice)
    λ0 = find_zero(f, 0;)

    return λ0
end

function likelihood_hazard_mother(γ::Function, θ, q; data, inf, tslice, x0)
    out = 0.0

    for cyc in data
        r = 0:length(cyc)-1
        ts_ = r .* tslice 

        pairs = zip(cyc, ts_)
        int, last_ = trapz(x -> γ(x[1], θ, x[2]), pairs, tslice) 
        out += log(last_[1])
        out -= int
    end
   
    if !isnothing(x0)
        out -= q*norm(x0 - θ)^2
    end

    if !isnan(out) && !isinf(out) 
        return out
    else
        return inf
    end
end

function optimisation(f; 
    x0, 
    rn, 
    pbounds, 
    data, 
    likelihood, 
    max_evals, 
    n_initial_points=10, 
    name="", 
    strain, 
    restarts, 
    tslice, 
    kwargs...)

    println("Saving checkpoint at data/sim/checkpoint_hazard_$(strain)_$name.pkl")
    checkpoint_saver = CheckpointSaver("data/sim/checkpoint_hazard_$(strain)_$name.pkl", store_objective=false)

    opt = nothing
    x0_ = x0 
    display(x0_)
    
    for n in 1:restarts
        println("Epoch $n out of $restarts")
        opt = gp_minimize(likelihood, pbounds, verbose=true, 
            n_calls=max_evals, n_initial_points=n_initial_points, 
            x0=x0_,
            initial_point_generator="lhs", n_jobs=Base.Threads.nthreads(), callback=checkpoint_saver)
        x0_ = opt.x
        val = pyconvert(Float64, opt.fun)
        println("Best value $val")
    end

    println("Saved the final checkpoint at data/sim/checkpoint_hazard_$(strain)_$name.pkl")
    println("Optimisation finished")
    return opt
end

function fit_spline(hazard, paramrn; 
    data, 
    likelihood, 
    name, 
    strain, 
    tslice, 
    x0=nothing, 
    q0=(nothing, true),
    qbounds=nothing, 
    pbounds, 
    inf=-100000, 
    max_evals,
    restarts,
    kwargs...)

    if !isnothing(x0)
        x0_ = pylist(x0)  
    else
        x0_ = pybuiltins.None
    end

    nparams = length(pbounds)
    γ = AgentBasedCellss.gen_division_rate_function(hazard, paramrn)
    lik(p, q) = -likelihood(γ, pyconvert.(Float64, p), q; 
        data=data, inf=inf, tslice=tslice, x0=x0)

    if q0[2]
        pbounds_ = (pbounds..., qbounds...)
        liks_ = p -> lik(p[1:end-1], p[end])
    else
        pbounds_ = (pbounds...,)
        liks_ = p -> lik(p, q0[1])
    end

    if !isnothing(q0[1]) && q0[2]
        x0_ = pylist([x0_..., q0[1]])
    end


    opt = optimisation(
        hazard;
        rn=paramrn,
        pbounds = pbounds_,
        x0 = x0_, 
        data = data,
        likelihood=liks_,
        n_initial_points=10,
        tslice=tslice,
        name="$(name)_$nparams",
        strain=strain,
        restarts=restarts,
        max_evals=max_evals,
        kwargs...)
    return opt
end

function fit_spline_mcmc(hazard, paramrn; priorσ, sampler, data, likelihood, name, strain, tslice, x0=nothing, nsamples)
    γ = AgentBasedCellss.gen_division_rate_function(hazard, paramrn)
    loglik(p) = likelihood(γ, p, 0.0; data=data, inf=-Inf, tslice=tslice, x0=nothing)

    @model function model()
        prior ~ MvNormal(x0, priorσ)
        logl = loglik(prior)
        Turing.@addlogprob!(logl)
        return logl
    end

    model = model()

    opt = sample(model, sampler, nsamples; chain_type=Chains)
    serialize("data/sim/chainfile_$(strain)_$(name)_$(nsamples).jls", opt)
    println("Saved the chain at data/sim/chainfile_$(strain)_$(name)_$(nsamples).jls")
    return opt
end

function std_range(data, nstd::Float64; length::Int64=5) 
    # range containing nstd to either side of the mean
    m = mean(data)
    s = std(data)
    return range(m-nstd*s, stop=m+nstd*s, length=length)
end

function std_range(data, nstd::Tuple{Float64, Float64}; length::Int64=5)
    # range containing nstd to either side of the mean
    m = mean(data)
    s = std(data)
    return range(m-nstd[1]*s, stop=m+nstd[2]*s, length=length)
end

function compute_λ(dist, tspan)
    f(λ) = 2 - quadgk(t -> exp(λ*t) * pdf(dist,t), tspan...)[1]
    λ0 = find_zero(f, 0)
end

function itpp_spline(x, knots::AbstractRange{Float64}, itpvals::Vector{Float64})
    exp(extrapolate(
        Interpolations.scale(
            interpolate(itpvals, BSpline(Cubic(Line(OnGrid())))),
            #interpolate(itpvals, BSpline(Cubic())),
            knots),
        Interpolations.Flat()
       )(x))
end

function make_itpp_spline(knots::AbstractRange{Float64}, itpvals)
    extrapolate(
        Interpolations.scale(
#            interpolate([itpvals...], BSpline(Cubic(Line(OnGrid())))),
            interpolate([itpvals...], BSpline(Linear())),
            knots),
        Interpolations.Flat()
       )
end

mutable struct Interp 
    knots
    values
    iterp

    function Interp(values, knots)
        interp = make_itpp_spline(knots, values)
        return new(knots, values, interp) 
    end
end

function (itp::Interp)(x, values...)
    if itp.values == values 
        return exp(itp.iterp(x))
    else
        itp.iterp = make_itpp_spline(itp.knots, values)
        itp.values = values
        return exp(itp.iterp(x))
    end
end

function equantiles(x; quantile=0.95)
    eff = ecdf(x)
    display(eff(37.5))
    flow(q) = 0.025 - eff(q)
    fhigh(q) = 0.975 - eff(q)
    
    qlow = find_zero(flow, 100.0)
    display(qlow)
    qhigh = find_zero(fhigh, 100.0)

    return qlow, qhigh
end

function interdiv_dist_kde_pop(interdiv_times; time_slice, bounds, bandwidth)
    x_ = interdiv_times .* time_slice 
    return CumulativeKDEPop(InterpKDELin(kde(x_; bandwidth=bandwidth, boundary=bounds)), tend=bounds[end])
end

function interdiv_dist_kde_mother(interdiv_times; time_slice, bounds, bandwidth)
    x_ = interdiv_times .* time_slice 
    return CumulativeKDE(InterpKDELin(kde(x_; bandwidth=bandwidth, boundary=bounds)), tend=bounds[end])
end
