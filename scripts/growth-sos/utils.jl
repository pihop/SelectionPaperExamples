using Clustering
using GLM
using PaddedViews
using StatsBase
using Intervals
using RollingFunctions
using KernelDensity
using QuadGK
import Distributions: cdf, logpdf, logccdf

function make_padded_matrix(df)
    mlen = maximum(length.(df))
    pdv = PaddedView.(0.0, df, Ref((1,mlen)))
    return reduce(vcat, pdv)
end

function make_gfp_matrix(mat_area, mat_gfp)
    return mat_area .* mat_gfp
end

function remove_consec(div_idxs)
    res = [div_idxs[1],]
    for i in 2:length(div_idxs)
        if div_idxs[i]-1 != div_idxs[i-1]
            push!(res, div_idxs[i])
        end
    end
    return res
end

function rmnan(array)
    return filter(x -> !isnan(x), array)
end

function div_time_analysis(mat; window=5, threshold=(-1.0, 0.0))
    differences = diff(mat; dims=2) ./ mat[:,1:end-1]
    sz_diff = size(differences)
    sz_mat = size(mat)
    idxs = findall(x -> threshold[1] < x < threshold[2], differences)
    
    idxs_to_remove = []
    for i in idxs 
        wd = CartesianIndices((i.I[1]:i.I[1], i.I[2]:min(i.I[2] + window, sz_diff[2])))
        if length(wd) > 1
            rnmean = reduce(hcat, runmean(vec(differences[wd]), length(vec(differences[wd]))))'
            midx = map(x -> rnmean[x[1]] < 0.0 ? x[2] : nothing, collect(enumerate(collect(wd))))
            push!(idxs_to_remove, filter(x -> !isnothing(x) && x != i, midx)...)
        end
    end
    mothers = filter(x -> !in(x, idxs_to_remove), idxs)

    md_pairs = []
    for i in mothers
        rowm_idxs = filter(x -> x.I[1] == i.I[1] && x > i, mothers)
        if !isempty(rowm_idxs)
            nextm_idx = rowm_idxs[1]

            wd = CartesianIndices((i.I[1]:i.I[1], i.I[2]:nextm_idx.I[2]))
            lens =  mat[wd]
            if length(wd) > 1
                bidx = findmin(vec(lens))
                bidx = wd[bidx[2]] 
                push!(md_pairs, (i, bidx))
            end
        end
    end

    dm_pairs = []
    mothers = idxs 
    daughters = getindex.(md_pairs, 2)
    for d in daughters
        idxm_ = findfirst(x -> x.I[1] == d.I[1] && x > d, mothers)
        if idxm_ != nothing
            push!(dm_pairs, (d, mothers[idxm_]))
        end
    end

    return md_pairs, dm_pairs
#    return getindex.(res, 3), getindex.(res,1), getindex.(res,2)#differences[idxs], idxs, idxs_b
end


function division_times(mat)
    md_idxs_ = []
    birth_idxs_ = []
    div_times_ = []
    div_idxs_ = []
    differences = diff(mat; dims=2) ./ mat[:,1:end-1]

    for row in eachrow(differences)
        # Cluster
        r_ = rmnan(row)
        kms = kmeans(reshape(r_, length(r_), 1)', 2)
        cidx = findmin(vec(kms.centers))
        div_idxs = map(x -> x[2] == cidx[2] ? x[1] : nothing, enumerate(kms.assignments)) 
        filter!(x -> x != nothing, div_idxs)
        div_idxs = remove_consec(div_idxs)
        push!(div_times_, diff(div_idxs))
        push!(div_idxs_, div_idxs)

        md_idxs = [] 
        birth_idxs = [] 
        for idx in div_idxs
            bidx = findnext(map(x->x>0.0, r_[idx+1:end]), 1)
            if bidx != nothing
                push!(md_idxs, (idx, bidx + idx))
                push!(birth_idxs, bidx + idx)
            end
        end
        push!(md_idxs_, md_idxs)
        push!(birth_idxs_, Int64.(birth_idxs))
    end

    return div_times_, div_idxs_, md_idxs_, birth_idxs_
end

function fluorecence_at_division(div_times, gfps)
    fluors = []
    for (div, gfp) in zip(div_times, eachrow(gfps))
        push!(fluors, gfp[div]...)
    end
    return fluors
end

function fluorescence_tuples(md_times, gfps)
    fluors = []
    for (age, div, gfp) in zip(getindex(md_times,1), getindex(md_times,3), eachrow(gfps))
        for (a,d) in zip(age, div)
            push!(fluors, (gfp[d[1]], gfp[d[2]], a))
        end
    end
    return fluors
end

function md_fluorescence(md_times, gfps; ftime)
    fluors = []
    for md in md_times
        if ftime[1] < interdiv(md) < ftime[2] 
            push!(fluors, (gfps[md[1]], gfps[md[2]]))
        end
    end
    return fluors
end

function interdiv(dm)
    return dm[2].I[2] - dm[1].I[2]
end

function interdivs(dm_pairs)
    res = []
    for dm in dm_pairs
        push!(res, dm[2].I[2] - dm[1].I[2])
    end
    return res
end

function birthdiv_fluorescence_time(dm_times, gfps; ftime)
    fluors = []
    for dm in dm_times
        if ftime[1] < interdiv(dm) < ftime[2] 
            push!(fluors, (interdiv(dm), gfps[dm[1]], gfps[dm[2]]))
        end
    end
    return fluors
end

function bootstrap_sample_var(values; N=5000, K=1000)
    len = length(values)
    samples = sample(values, (K, len))
    stats = map(r -> var(r), eachrow(samples))
    return mean(stats), var(stats) 
end

function conversion_fact_regression(mdtuples; nsplits=6, fltr=(0.0, 3.0))
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
            #push!(bins[midpoint(interval)], abs(tuple[1] - 2*tuple[2])) #[tuple[2], tuple[1] - tuple[2]]...)
            #push!(bins[midpoint(interval)], tuple[2])#[tuple[2], tuple[1] - tuple[2]]...)
            push!(bins[midpoint(interval)], [tuple[2], tuple[1] - tuple[2]]...)
        end
    end
    
    boots = bootstrap_sample_var.(values(bins))
    df = DataFrame(X=collect(keys(bins)), Y=var.(values(bins)), Z = getindex.(boots, 1), W = getindex.(boots,2), S = length.(values(bins)))
    ols = lm(@formula(Y ~ 0 + X), df) 
    return df, ols, coef(ols)[1] * 4
end

function make_cell_cycles(gfp, dm_divtimes)
    out = [] 
    for divt in dm_divtimes
        push!(out, vec(gfp[range(divt...)]))
    end
    return out
end

function interdiv_dist_gamma(interdiv_times; time_slice) 
    # Fit interdivision time distribution Gamma(a, b) via moment matching.
    
    # Convert time-slices to time in minutes.
    μ = mean(interdiv_times .* time_slice)
    σ = var(interdiv_times .* time_slice)

    # Mean and variance parametrisation of Gamma distribution.
    opt_init = [μ^2 / σ, σ / μ]

    return Gamma(opt_init...)
end

function interdiv_dist_kde(interdiv_times; time_slice)
    x_ = interdiv_times .* time_slice
#    b = KernelDensity.default_bandwidth(x_)
    fit1 = kde(x_) 
#    fit2 = MixtureModel(Normal.(x_, b))
    return InterpKDE(fit1) 
end

function cdf(kde::InterpKDE, x)
    integral, err = quadgk(s -> pdf(kde, s), 0, x, rtol=1e-8)
    return integral
end

function logpdf(kde::InterpKDE, x)
    return log(pdf(kde, x))
end

function logccdf(kde::InterpKDE, x)
    integral, err = quadgk(s -> pdf(kde, s), 0, x, rtol=1e-4)
    return log(1 - integral)
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
