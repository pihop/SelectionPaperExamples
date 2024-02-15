using Catalyst
using FiniteStateProjection
using DifferentialEquations.EnsembleAnalysis
using SparseArrays
using AgentBasedFSP
using JumpProcesses
using ImageFiltering

mutable struct StochasticDilutionModel
    rn::ReactionSystem
    init
    tspan
    ps
    steady_state

    function StochasticDilutionModel(rn, init, tspan, ps)
        new(rn, init, tspan, ps, [])
    end
end

#_StochasticDilutionModel(rn,  init, ps) = StochasticDilutionModel(rn, init, ps)
#Broadcast.broadcasted(::typeof(StochasticDilutionModel), rn, init, ps) = 
#    broadcast(_StochasticDilutionModel, Ref(rn), Ref(init), ps)

struct BirthDeathException <: Exception end
Base.showerror(io::IO, e::BirthDeathException) = print("Not a birth death process.")

function is_birth_death(rn::ReactionSystem)
    # Do we care if either birth or death process doesn't exist?
    return length(Catalyst.species(rn)) == 1
end

function birth_death_steady_state!(model::StochasticDilutionModel, truncation)
    # Before we do anything check that the rn is a birth death process.

    if !is_birth_death(model.rn)
        throw(BirthDeathException())
    end

    reacts = reactions(model.rn)
    rates = [r.rate for r in reacts]
    stoich = Catalyst.netstoichmat(model.rn)

    ratef = AgentBasedFSP.gen_division_rate_function.(rates, model.rn)
    eval_rates = [[f(state, model.ps, 0.0) for f in ratef] for state in 1:truncation[1]] 
    # currently evaluated for t = 0.0. Ok if rates are not time dependent.
    sumλ = sum.(map(x -> x[vec(stoich .== 1)], eval_rates))
    sumμ = sum.(map(x -> x[vec(stoich .== -1)], eval_rates)) .* (collect(1:truncation[1]) .- 1)

    cprods = cumprod([λ/sumμ[i+1] for (i, λ) in enumerate(sumλ[1:end-1])])
    
    π₀ = 1 / (1 + sum(cprods))
    πₖ = π₀ .* cprods 
    model.steady_state = vcat(π₀, πₖ)
end

function mean_steady_state(model::StochasticDilutionModel)
    return sum(model.steady_state .* (collect(1:length(model.steady_state)) .- 1))
end

_birth_death_steady_state!(model, truncation) = birth_death_steady_state!(model, truncation)
Broadcast.broadcasted(::typeof(birth_death_steady_state!), model, truncation) = 
    broadcast(_birth_death_steady_state!, model, Ref(truncation))

function stochastic_steady_state!(model::StochasticDilutionModel, truncation, solver)
    println("Computing steady state of stochastic dilution ...")
    fsp = FSPSystem(model.rn)  
    A = convert(SparseMatrixCSC, fsp, truncation, model.ps, 0)
    bndA = abs.(vec(sum(A, dims=1))) 
    bndA = reshape(bndA, truncation...)
    u0 = zeros(truncation...)
    u0[1] = 1.0
    sz = ones(size(u0)) 
    sz = sz ./ sum(sz)

    function fu!(du, u, p, t)
        du .= reshape(A*vec(u), size(u))
        du .+= sum(u .* bndA) .* sz 
    end

    prob = ODEProblem(fu!, u0, model.tspan, model.ps) 
    @time sol = solve(prob, solver)
    model.steady_state = sol.u[end]  
    println("Done")
end

function simulation_steady_state!(model::StochasticDilutionModel, init, tspan)
    println("Computing steady state of stochastic dilution ...")
    @time begin 
        dprob = DiscreteProblem(model.rn, init, tspan, model.ps)
        jprob = JumpProblem(model.rn, dprob, Direct(); save_positions=(false,false))
        jsol = solve(jprob, SSAStepper(); saveat=tspan[1]:1.:tspan[end])
        n = length(species(model.rn)) 
        Nmax = [maximum(getindex.(jsol.u, i)) + 1 for i in 1:n]
        resmat = zeros(Int64.(Nmax)...)
        for res in jsol.u
            resmat[Int64.(res .+ 1)...] += 1
        end
        resmat = resmat ./ sum(resmat)
        resmat = imfilter(resmat, Kernel.gaussian(1.5))
        model.steady_state = resmat
    end
end
