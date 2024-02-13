using IntervalArithmetic
using IntervalRootFinding
using Catalyst
using StaticArrays
import SciMLBase: ODEFunction

mutable struct EffectiveDilutionModel
    rn::ReactionSystem
    roots_function::Function
    bif_idx::Int64
    ps
    roots::Vector{Any}

    function EffectiveDilutionModel(rn, bif_idx)
        # 1D model
        odefun = ODEFunction(convert(ODESystem,rn),jac=true)
        n = length(species(rn))
        rootsf(ps) = 
            x -> begin
                odefun.f(SVector(x...), ps, 0.0)
            end
        new(rn, rootsf, bif_idx, nothing, [])
    end
end

function root_finding(model::EffectiveDilutionModel, params, bif_range; search_interval) 
    model.roots = []
    for p in bif_range
        _ps = [params...]
        _ps[model.bif_idx] = p 
        rts = roots(model.roots_function(_ps), search_interval)
        midpoints = mid.(getfield.(rts, :interval))
        pts = [(p, m) for m in midpoints]
        model.ps = params
        push!(model.roots, pts...)
    end
end
