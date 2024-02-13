using StatsBase
using Distributions
using GLM
using Intervals
import Distributions: cdf, logccdf, logpdf
using KernelDensity
using QuadGK
using DataFrames

struct LineageData
    lineages::Vector{DataFrame}
    mother_cells::DataFrame
    slice::Tuple{Int64, Int64} 
end

struct CellCycle
    cell::DataFrame
    cycle::DataFrame
end

function divided_cells(lineages, df)
    # Count cells at each time slice. 
    slices = union(getindex.(lineages, :, :Slice)...)
     
    count_cells = [] 
    mother_cells = DataFrame()
    for s in slices 
        frames = subset.(lineages, :Slice => x -> isequal.(x,s))
        push!(count_cells, countmap(vcat(getindex.(frames, :, :1)...)))
    end
    reverse!(count_cells)

    for i in 1:(length(count_cells) - 1)
        if length(count_cells[i]) < length(count_cells[i+1])
            for mother_idx in keys(count_cells[i])
                daughtercells = subset(df, :PreviousCell => x -> isequal.(x, mother_idx))
                if nrow(daughtercells) > 1
                    push!(mother_cells, subset(df, :1 => x -> isequal.(x, mother_idx))[1, :])
                end
            end
        end
    end

    return mother_cells
end

function compute_lineages(df; slice = [41, 121])
    lineages = []
    cells_at_end = subset(df, :LastIndex => x -> x .== 1) 

    for idx in cells_at_end[:, 1]
        # choose index to start with.
        df_ = DataFrame()
        while idx != 0
            push!(df_, df[idx, :])
            idx = df[idx, :PreviousCell]
        end
        push!(lineages, unique(subset(df_, :Slice => x -> slice[1] .≤ x .≤ slice[end])))
    end
    filter!(x -> !isempty(x), lineages)
    mother_cells = divided_cells(lineages, df)
    return LineageData(sort!.(lineages, :Slice), mother_cells, (slice[1], slice[end]))
end

function make_mother_time_pairs(df)
    return [[r[:1], r[:Slice]] for r in eachrow(df)]
end

function interdivision_times(lineages)
    if nrow(lineages.mother_cells) > 0
        mother_idxs = lineages.mother_cells[:, :1]
        df_ = unique(subset.(lineages.lineages, :1 => x -> in.(x, Ref(mother_idxs))))
        pairs_containers = make_mother_time_pairs.(df_)
        res = []
        for cont in pairs_containers
            push!(res, [[cont[i][1], cont[i][2] - cont[i-1][2]] for i in 2:length(cont)]...)
        end
        unique!(res)
        return getindex.(res, 2)
    end
end

function joint_fluor_times(lineages)
    if nrow(lineages.mother_cells) > 0
        mother_idxs = lineages.mother_cells[:, :1]
        df_ = unique(subset.(lineages.lineages, :1 => x -> in.(x, Ref(mother_idxs))))
        pairs_containers = make_mother_time_pairs.(df_)
        res = []
        for cont in pairs_containers
            push!(res, [[cont[i][1], cont[i][2] - cont[i-1][2]] for i in 2:length(cont)]...)
        end
        unique!(res)
        fluors = [fluorescence_intensity.(eachrow(filter(x -> x[:1] == r[1], lineages.mother_cells)))[1] for r in res]
        return collect.(collect(zip(getindex.(res, 2), fluors, fluors)))
    end
end

function get_slice(row)
    return row[:Slice]
end

function make_cell_cycle(lineages)
    if nrow(lineages.mother_cells) > 0 
        mother_idxs = lineages.mother_cells[:, :1]
        df_ = subset.(lineages.lineages, :1 => x -> in.(x, Ref(mother_idxs)))
        slices = [get_slice.(eachrow(d)) for d in df_] 
        cell_cycles = []
        for (slice, d) in zip(slices, lineages.lineages)
            for i in 2:length(slice)
                push!(cell_cycles, subset(d, :Slice => x -> slice[i-1] .< x .≤ slice[i]))
            end
        end
        return unique(cell_cycles)
    end
end

function fluor_at_age(cycle; n, cf)
    res = []
    if nrow(cycle) ≥ n
        return [fluorescence_intensity(cycle[n,:]) / cf,  nrow(cycle) == n]
    end
end

function fluorescence_intensity(row)
    return 0.067 * row[:Area] * ((row[:Mean] - row[:Background]) / 500)
end

function fluorescence_tuples(lineages)
    if nrow(lineages.mother_cells) > 0
        mother_idxs = lineages.mother_cells[:, :1]
        df_daughter = unique(vcat(subset.(lineages.lineages, :PreviousCell => x -> in.(x, Ref(mother_idxs)))...))

        tuples = []
        
        for mother in eachrow(lineages.mother_cells)
            mother_fluorescence = fluorescence_intensity(mother)
            daughters = subset(df_daughter, :PreviousCell => x -> isequal.(x, mother[:1]))
            daughter_fluorescence = fluorescence_intensity.(eachrow(daughters))

            push!(tuples, vcat(mother_fluorescence, daughter_fluorescence))
        end
        
        return filter(x -> length(x) == 3, tuples)
    end
end

function fluorescence_trajectory(cell_cycle)
    return fluorescence_intensity.(eachrow(cell_cycle)) 
end

function div_fluorescence(lineages)
    if nrow(lineages.mother_cells) > 0
        intensities = fluorescence_intensity.(eachrow(lineages.mother_cells))
        return intensities 
    end
end

function div_joint(lineages)
    if nrow(lineages.mother_cells) > 0
        intensities = fluorescence_intensity.(eachrow(lineages.mother_cells))
        return intensities 
    end
end

function birth_fluorescence(lineages)
    if nrow(lineages.mother_cells) > 0
        mother_idxs = lineages.mother_cells[:, :1]
        df_daughter = unique(vcat(subset.(lineages.lineages, :PreviousCell => x -> in.(x, Ref(mother_idxs)))...))
        intensities = fluorescence_intensity.(eachrow(df_daughter))
        return intensities 
    end
end

function interdiv_dist_gamma(interdiv_times; time_slice) 
    # Fit interdivision time distribution Gamma(a, b) via moment matching.
    
    # Convert time-slices to time in minutes. Each slice is 5 minutes.
    μ = mean(interdiv_times .* time_slice)
    σ = var(interdiv_times .* time_slice)

    # Mean and variance parametrisation of Gamma distribution.
    opt_init = [μ^2 / σ, σ / μ]

    return Gamma(opt_init...)
end
#
#function interdiv_dist_kde(interdiv_times; time_slice, bandwidth)
#    x_ = interdiv_times .* time_slice
#    fit1 = kde(x_, bandwidth=6) 
#    return InterpKDE(fit1) 
#end
#
#function logpdf(kde::InterpKDE, x)
#    return log(pdf(kde, x))
#end
#
#function logccdf(kde::InterpKDE, x)
#    integral, err = quadgk(s -> pdf(kde, s), 0, x, rtol=1e-4)
#    return log(1 - integral)
#end

function aggr_pdf(dist; rnge)
    res = [] 
    for i in 2:length(rnge)
        push!(res, cdf(dist, rnge[i]) - cdf(dist, rnge[i-1]))
    end
    return res
end

function bootstrap_sample_var(values; N=5000, K=1000)
    len = length(values)
    samples = sample(values, (K, len))
    stats = map(r -> var(r), eachrow(samples))
    return mean(stats), var(stats) 
end

function conversion_fact_regression(mdtuples; nsplits=6, fltr=(10, 30))
    # Set up the partitions.
    minm = minimum(filter(x -> fltr[1] < x < fltr[2], first.(mdtuples)))
    maxm = maximum(filter(x -> fltr[1] < x < fltr[2], first.(mdtuples)))
    diff = maxm - minm
    lrange = range(minm, stop=maxm - diff/nsplits, length=nsplits)
    urange = range(minm + diff/nsplits, stop=maxm, length=nsplits)
    intervals = Intervals.Interval{Closed, Closed}.(lrange, urange)

    bins = Dict(x => [] for x in midpoint.(intervals)) 

    # Partition. 
    for tuple in mdtuples
        # Find the interval
        interval = find_interval(tuple[1], intervals)
        # find the center
        if interval != nothing
            push!(bins[midpoint(interval)], [tuple[2], tuple[3]]...)
        end
    end
    
    boots = bootstrap_sample_var.(values(bins))
    df = DataFrame(X=collect(keys(bins)), Y=var.(values(bins)), Z = getindex.(boots, 1), W = getindex.(boots,2), S = length.(values(bins)))
    ols = lm(@formula(Y ~ 0 + X), df, wts=df.S) 

    return df, ols, 4*coef(ols)[1]
end

function rmnothing(array)
    return filter(x -> !isnothing(x), array)
end

function rmnan(array)
    return filter(x -> !isnan(x), array)
end

# Interval
midpoint(i::Intervals.Interval) = (i.last + i.first)/2

function find_interval(a, intervals)
    for interval in intervals
        if a in interval
            return interval
        end
    end
end
