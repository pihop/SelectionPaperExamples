using AgentBasedCells
using CairoMakie
import CairoMakie: Axis
using StatsBase
using LinearAlgebra
using Peaks
using PaddedViews

using Colors
using ColorSchemes
colors = ColorSchemes.Egypt.colors
colors = vcat(colors, RGB(81/255,179/255,179/255))

transp = 0.2

CairoMakie.activate!(type = "svg")

function mylt((x1,y1), (x2,y2))
    y1>y2
end

function make_normalised_bar_data(array)
    counts = countmap(array)
    vals = collect(values(counts))
    vals = vals ./ sum(vals)
    return collect(zip(collect(Float64.(keys(counts))), vals))
end

function plot_selection(simulation, analytical, experiment;
    idx=1, 
    xlimsbirth, 
    xlimsdiv, 
    ylimsbirth, 
    ylimsdiv, 
    divtbins, 
    divtrange,
    selection=nothing, 
    tag="selection", 
    sumdim, 
    traj_idxs=[])

    CairoMakie.activate!()
    # Comparison of simulation and analytical results.
    tspan_simulation = simulation.problem.tspan
    tspan_analytical = analytical.problem.tspan

    size_inches = (17, 7)
    size_pt = 72 .* size_inches
    figure = Figure(resolution=size_pt, fontsize=19)
    ax_divt = Axis(figure[1,1], title="Interdivision time distribution") 
    #ax_birth = Axis(figure[1,2], title="Birth distribution")
    ax_div = Axis(figure[1,2], title="Division distribution")
    ax_joint = Axis(figure[1,3], title="First passage time distribution")

    divn_times = normalize(fit(Histogram, simulation.results[:division_times_all], divtbins);
        mode=:pdf)
    stairs!(ax_divt, collect(midpoints(divn_times.edges[1])), divn_times.weights; 
        color=colors[1], step=:center, linewidth=1)
    barplot!(ax_divt, collect(midpoints(divn_times.edges[1])), divn_times.weights;  
        color=(colors[1], transp), gap=0.0, dodge_gap=0.0, 
        linecolor=:match)

    lines!(ax_divt, 
        collect(tspan_analytical[1]:0.001:divtbins[end]), 
        t -> division_time_dist(analytical)(t); 
        color=colors[1], label="Analytical (no selection)", linewidth=2, overdraw=true) 

    xlims!(ax_divt, (tspan_analytical[1], divtbins[end]))
    ylims!(ax_divt, low=0.0)
    ax_divt.xlabel = "Time"
    ax_divt.ylabel = "Probability density"

    bar_division = make_normalised_bar_data(
        getindex.(simulation.results[:molecules_at_division], idx))
    sort!(bar_division; by=x->x[1])
    barplot!(ax_div, getindex.(bar_division, 1), getindex.(bar_division, 2); 
        color=(colors[1], transp), strokewidth = 0.0, gap=0.0, dodge_gap=0.0, strokecolor = (colors[1], transp),
        linecolor=(colors[1], transp))
    stairs!(ax_div, getindex.(bar_division, 1), getindex.(bar_division, 2); 
        color=colors[1], step=:center, linewidth=1)

    dd = sum(analytical.results[:division_dist_ancest]; dims=sumdim)
    lines!(ax_div, collect(0:length(dd)-1), vec(dd); color=colors[1], linewidth=2, overdraw=true) 
#    dd = sum(analytical.results[:division_dist_ancest]; dims=sumdim)
#    lines!(ax_div, collect(0:length(dd)-1), vec(dd); color=colors[1], linewidth=2, overdraw=true) 

    xlims!(ax_div, xlimsdiv)
    ylims!(ax_div, ylimsdiv)

    ax_div.xlabel = "Protein counts"
    ax_div.ylabel = "Probability density"
    
    zs1 = hcat(sum.(fpt_dist_ancest(analytical).(divtrange), dims=sumdim)...)
    cmap1 = range(alphacolor(colors[2], 0.0), stop=ARGB(colors[2], 0.8), length=15)
    contourf!(ax_joint, divtrange, 1:length(analytical.results[:birth_dist]), zs1'; colormap=cmap1, linewidth=2)
    xlims!(ax_joint, (divtrange[1], divtrange[end]))
    ylims!(ax_joint, xlimsdiv)

    ax_joint.xlabel = "Time"
    ax_joint.ylabel = "Protein counts"

    for i in traj_idxs
        divtime = cell_division_age(simulation.results[:all_population][i])
        traj = simulation.results[:all_population][i].sim(0.0:0.01:divtime)
        lines!(ax_joint, 0.0:0.01:divtime, getindex.(traj.u, idx); color=(colors[2], 0.7), linewidth=2)
        scatter!(ax_joint, traj.t[end], traj.u[end][idx]; color=(colors[2], 0.8), strokewidth=2, marker=:xcross)
    end

    sim_colors = [PolyElement(color=colors[1]), PolyElement(color=(colors[2]))]
    analytical_colors = [LineElement(color=colors[1]), LineElement(color=(colors[2]))]

#    Legend(figure[2,1:3], [sim_colors, analytical_colors], [["No selection", "Selection"], ["No selection", "Selection"]], ["Simulation", "Analytical"]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true)

    prefix = DrWatson.default_prefix(experiment)
    sname = savename(experiment; allowedtypes = (Vector, String, Int, Float64))
    sname = sname * tag
    path = mkpath("$(plotsdir())/$prefix/")
    save("$path/$(sname)_plot.pdf", figure, pt_per_unit = 1)
    save("$path/$(sname)_plot.svg", figure, pt_per_unit = 1)
    return figure
end


function plots_selection_compare(
    simulation, analytical, analytical_mother, experiment; 
    idx=1, 
    xlimsbirth, 
    xlimsdiv, 
    ylimsbirth, 
    ylimsdiv, 
    divtbins, 
    divtrange,
    fptlimsx,
    fptlimsy,
    selection=nothing, 
    tag="selection", 
    sumdim, 
    acolors,
    acolorsm,
    scolors,
    traj_idxs=[])

    CairoMakie.activate!()
    # Comparison of simulation and analytical results.
    tspan_simulation = simulation[1].problem.tspan
    tspan_analytical = analytical[1].problem.tspan

    size_inches = (14, 7)
    size_pt = 72 .* size_inches
    figure = Figure(resolution=size_pt, fontsize=19)
    #ax_divt = Axis(figure[1,1], title="Interdivision time distribution") 
    #ax_birth = Axis(figure[1,2], title="Birth distribution")
    ax_div = Axis(figure[1,1], title="Division protein distribution")
    ax_joint = Axis(figure[1,2], title="Division distribution")

    for (sim, scol) in zip(simulation, scolors)
        divn_times = normalize(fit(Histogram, sim.results[:division_times_ancest], divtbins);
            mode=:pdf)
#        stairs!(ax_divt, collect(midpoints(divn_times.edges[1])), divn_times.weights; 
#            color=colors[i], step=:center, linewidth=1)
#        barplot!(ax_divt, collect(midpoints(divn_times.edges[1])), divn_times.weights;  
#            color=(colors[i], transp), gap=0.0, dodge_gap=0.0, 
#            linecolor=:match)

        bar_division = make_normalised_bar_data(
            getindex.(sim.results[:molecules_at_division], idx))
        sort!(bar_division; by=x->x[1])
        barplot!(ax_div, getindex.(bar_division, 1), getindex.(bar_division, 2); 
            color=(colors[scol], transp), strokewidth = 0.0, gap=0.0, dodge_gap=0.0, strokecolor = (colors[scol], transp),
            linecolor=(colors[scol], transp))
        stairs!(ax_div, getindex.(bar_division, 1), getindex.(bar_division, 2); 
            color=colors[scol], step=:center, linewidth=1)
    end

    for (analy, acol) in zip(analytical, acolors)
            linestyle = :solid
#        lines!(ax_divt, 
#            collect(tspan_analytical[1]:0.001:divtbins[end]), 
#            t -> division_time_ancest(analy)(t); 
#            color=colors[i], label="Analytical (no selection)", linewidth=3, overdraw=true, linestyle=linestyle) 

        dd = sum(analy.results[:division_dist_ancest]; dims=sumdim)
        lines!(ax_div, collect(0:length(dd)-1), vec(dd); color=colors[acol], linewidth=3, overdraw=true, linestyle = linestyle) 
    end

    for (analy, acol) in zip(analytical_mother, acolorsm)
        linestyle = :dash
        dd = sum(analy.results[:division_dist_ancest]; dims=sumdim)
        lines!(ax_div, collect(0:length(dd)-1), vec(dd); color=colors[acol], linewidth=3, overdraw=true, linestyle = linestyle) 
    end

#    xlims!(ax_divt, tspan_analytical[1], divtrange[end])
#    ylims!(ax_divt, low=0.0)
#    ax_divt.xlabel = "Interdivision time"
#    ax_divt.ylabel = "Probability density"

    xlims!(ax_div, xlimsdiv)
    ylims!(ax_div, ylimsdiv)

    ax_div.xlabel = "Protein counts"
    ax_div.ylabel = "Probability density"
    tspan = divtrange[1]:0.001:divtrange[end]
    
    zs1_ = sum.(fpt_dist_ancest(analytical[1]).(tspan), dims=sumdim)
    zs1 = reshape(hcat(zs1_...), (length(zs1_[1]), length(tspan)))
    cmap1 = range(alphacolor(colors[acolors[1]], 0.5), stop=ARGB(colors[acolors[1]], 1.0), length=5)
    contour!(ax_joint, tspan, 1:length(zs1_[1]), zs1'; colormap=cmap1, linewidth=3, levels=range(0., 0.04, length=5))

    zs2_ = sum.(fpt_dist_ancest(analytical[2]).(tspan), dims=sumdim)
    zs2 = reshape(hcat(zs2_...), (length(zs2_[1]), length(tspan)))
    cmap2 = range(alphacolor(colors[acolors[2]], 0.0), stop=ARGB(colors[acolors[2]], 1.0), length=5)
    contour!(ax_joint, tspan, 1:length(zs2_[1]), zs2'; colormap=cmap2, linewidth=3, levels=range(0., 0.04, length=5))
    xlims!(ax_joint, fptlimsx) 
    ylims!(ax_joint, fptlimsy) 

    ax_joint.xlabel = "Interdivision time"
    ax_joint.ylabel = "Protein counts"

    for i in traj_idxs
        divtime = cell_division_age(simulation[1].results[:all_population][i])
        traj = simulation[1].results[:all_population][i].sim(0.0:0.01:divtime)
        lines!(ax_joint, 0.0:0.01:divtime, getindex.(traj.u, idx); color=(:gray, 0.8), linewidth=2)
        scatter!(ax_joint, traj.t[end], traj.u[end][idx]; color=(:gray, 0.8), strokewidth=2, marker=:xcross)
    end

    sim_colors = [PolyElement(color=colors[scol]) for scol in scolors]
    analytical_colors = [LineElement(color=colors[acol], linewidth=3) for acol in acolors]

    analytical_mother_colors = [LineElement(color=(colors[acol]), linestyle=:dash, linewidth=3) for acol in acolorsm]

    Legend(figure[2,1:2], 
        [sim_colors, analytical_colors, analytical_mother_colors], 
        [["No selection", "Selection"], 
         ["No selection", "Selection"], ["No selection", "Selection"]], 
        ["Population simulation", "Analytical (population)", "Analytical (mother machine)"]; 
        orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true)

    prefix = DrWatson.default_prefix(experiment)
    sname = savename(experiment; allowedtypes = (Vector, String, Int, Float64))
    sname = sname * tag
    path = mkpath("$(plotsdir())/$prefix/")
    save("$path/$sname.pdf", figure, pt_per_unit = 1)
    save("$path/$sname.svg", figure, pt_per_unit = 1)
    return figure
end

function plots(simulation, analytical, experiment; idx=1, xlimsbirth, xlimsdiv, ylimsbirth, ylimsdiv, divtbins, selection=nothing, tag="", sumdim)
    CairoMakie.activate!()
    # Comparison of simulation and analytical results.
    tspan_simulation = simulation.problem.tspan
    tspan_analytical = analytical.problem.tspan

    figure = Figure(resolution = (1100, 500), fontsize = 16)
    ax_divt = Axis(figure[1,1]) 
    ax_birth = Axis(figure[1,2])
    ax_div = Axis(figure[1,3])

    divn_times_all = normalize(fit(Histogram, simulation.results[:division_times_all], divtbins);
        mode=:pdf)
    stairs!(ax_divt, collect(midpoints(divn_times_all.edges[1])), divn_times_all.weights; 
        color=colors[1], step=:center, linewidth=1)
    barplot!(ax_divt, collect(midpoints(divn_times_all.edges[1])), divn_times_all.weights;  
        color=(colors[1], transp), gap=0.0, dodge_gap=0.0, 
        linecolor=:match)

#    divn_times_ancest = normalize(fit(Histogram, simulation.results[:division_times_ancest], divtbins);
#        mode=:pdf)
#    stairs!(ax_divt, collect(midpoints(divn_times_ancest.edges[1])), divn_times_ancest.weights; 
#        color=colors[2], step=:center, linewidth=1)
#    barplot!(ax_divt, collect(midpoints(divn_times_ancest.edges[1])), divn_times_ancest.weights;  
#        color=(colors[2], transp), gap=0.0, dodge_gap=0.0, 
#        linecolor=:match, strokewidth=0.0)

    lines!(ax_divt, 
        collect(tspan_analytical[1]:0.001:divtbins[end]), 
        t -> division_time_dist(analytical)(t); 
        color=colors[1], label="Analytical", linewidth=2, overdraw=true) 
#    lines!(ax_divt, 
#        collect(tspan_analytical[1]:0.001:divtbins[end]), 
#        t -> 2*division_time_dist(analytical)(t)*exp(-analytical.results[:growth_factor] * t); 
#        color=colors[2], label="Analytical", linewidth=2, overdraw=true) 

    xlims!(ax_divt, tspan_analytical[1], divtbins[end])
    ylims!(ax_divt, low=0.0)
    ax_divt.xlabel = "Time at division"
    ax_divt.ylabel = "Probability density"

    bar_birth_all = make_normalised_bar_data(getindex.(simulation.results[:molecules_at_birth_all], idx))
    sort!(bar_birth_all; by=x->x[1])
    if selection != nothing
        ax_sel = Axis(figure[1, 1],
            width=Relative(0.4),
            height=Relative(0.3),
            halign=0.85,
            valign=0.9,
            backgroundcolor=:white,
            xlabelsize=14,
            xticksize=14,
            ylabelsize=14,
            yticksize=14,
           )
        translate!(ax_sel.blockscene, 0, 0, 1000)
        lines!(ax_sel, getindex.(bar_birth_all, 1), selection.(getindex.(bar_birth_all, 1)); 
            linewidth=2, overdraw=true, color=colors[3])
        ax_sel.xlabel = "Molecule numbers"
        ax_sel.ylabel = "Selection strength"
        translate!(ax_sel.scene, 0, 0, 10)
        xlims!(ax_sel; low=0.0)  
        ylims!(ax_sel; low=0.0)  
    end

    bar_birth_all = make_normalised_bar_data(getindex.(simulation.results[:molecules_at_birth_all], idx))
    sort!(bar_birth_all; by=x->x[1])
    barplot!(ax_birth, getindex.(bar_birth_all, 1), getindex.(bar_birth_all, 2); 
        width=1.0, color=(colors[1], transp), strokewidth = 0.0, gap=0.0, dodge_gap=0.0,
        linecolor=(colors[1], transp))
    stairs!(ax_birth, getindex.(bar_birth_all, 1), getindex.(bar_birth_all, 2); 
        color=colors[1], step=:center, linewidth=1)

#    bar_birth = make_normalised_bar_data(getindex.(simulation.results[:molecules_at_birth], idx))
#    sort!(bar_birth; by=x->x[1])
#    barplot!(ax_birth, getindex.(bar_birth, 1), getindex.(bar_birth, 2);
#        width=1.0, color=(colors[2], transp), strokewidth = 0.0, gap=0.0, dodge_gap=0.0, 
#        linecolor=(colors[2], transp))
#    stairs!(ax_birth, getindex.(bar_birth, 1), getindex.(bar_birth, 2); 
#        color=colors[2], step=:center, linewidth=1)

    bd = sum(analytical.results[:birth_dist]; dims=sumdim)
    lines!(ax_birth, collect(0:length(bd)-1), vec(bd); color=colors[1], linewidth=2, overdraw=true)

#    bd_ancest = sum(analytical.results[:birth_dist_ancest]; dims=sumdim)
#    lines!(ax_birth, collect(0:length(bd_ancest)-1), vec(bd_ancest); color=colors[2], linewidth=2, overdraw=true) 

    ax_birth.xlabel = "Molecule numbers at birth"
    ax_birth.ylabel = "Probability density"
#    ax_birth.title = "Distribution of molecules at birth"
    xlims!(ax_birth, xlimsbirth)
    ylims!(ax_birth, ylimsbirth)

    bar_division_all = make_normalised_bar_data(
        getindex.(simulation.results[:molecules_at_division_all], idx))
    sort!(bar_division_all; by=x->x[1])
    barplot!(ax_div, getindex.(bar_division_all, 1), getindex.(bar_division_all, 2); 
        color=(colors[1], transp), strokewidth = 0.0, gap=0.0, dodge_gap=0.0, strokecolor = (colors[1], transp),
        linecolor=(colors[1], transp))
    stairs!(ax_div, getindex.(bar_division_all, 1), getindex.(bar_division_all, 2); 
        color=colors[1], step=:center, linewidth=1)

#    bar_division = make_normalised_bar_data(
#        getindex.(simulation.results[:molecules_at_division], idx))
#    sort!(bar_division; by=x->x[1])
#    barplot!(ax_div, getindex.(bar_division, 1), getindex.(bar_division, 2); 
#        color=(colors[2], transp), strokewidth = 0.0, gap=0.0, dodge_gap=0.0, strokecolor = (colors[2], transp),
#        linecolor=(colors[2], transp))
#    stairs!(ax_div, getindex.(bar_division, 1), getindex.(bar_division, 2); 
#        color=colors[2], step=:center, linewidth=1)
    dd = sum(analytical.results[:division_dist]; dims=sumdim)
    lines!(ax_div, collect(0:length(dd)-1), vec(dd); color=colors[1], linewidth=2, overdraw=true) 

#    dd_ancest = sum(analytical.results[:division_dist_ancest]; dims=sumdim)
#    lines!(ax_div, collect(0:length(dd_ancest)-1), vec(dd_ancest); color=colors[2], linewidth=2, overdraw=true) 

    xlims!(ax_div, xlimsdiv)
    ylims!(ax_div, ylimsdiv)

    ax_div.xlabel = "Protein counts at division"
    ax_div.ylabel = "Probability density"
#    ax_div.title = "Distribution of molecules at division"

#    sim_colors = [PolyElement(color=colors[1]), PolyElement(color=(colors[2]))]
    sim_colors = [PolyElement(color=colors[1]), ]
    Legend(figure[2,2], ax_divt, orientation=:vertical, framevisible = false, tellwidth = false, tellheight = true)
    #Legend(figure[2,1:2], sim_colors, ["Simulation", "Simulation ancestral"], orientation=:horizontal, framevisible = false, tellwidth = false, tellheight = true)
    Legend(figure[3,2], sim_colors, ["Simulation", ], orientation=:vertical, framevisible = false, tellwidth = false, tellheight = true)

    prefix = DrWatson.default_prefix(experiment)
    sname = savename(experiment; allowedtypes = (Vector, String, Int, Float64))
    sname = sname * tag
    path = mkpath("$(plotsdir())/$prefix/")
    save("$path/$sname.pdf", figure, pt_per_unit = 1)
    save("$path/$sname.svg", figure, pt_per_unit = 1)
    return figure
end

function filter_modes(modes)
    res = []
    for m in modes
        flt = filter(x->x[2]>1e-6, m)
        if !isempty(flt)
            push!(res, flt)
        else
            push!(res, [m[1],])
        end
    end
    return res
end

function process_modes(modes, pspans)
    # Find the min index. This is where we switch the major modes. 
    modes = filter_modes(modes)
    maxm = maximum(length.(modes))

    paths = [] 
    for i in 1:maxm
        p = []  
        for m in modes
            push!(p, Float64.([m[max(1, length(m) - i + 1)]...]))
        end
        push!(paths, p)
    end

    verts = []

    if length(paths) > 1
        for i in 1:length(modes)  
            p_ = getindex.(paths, i)
            if !allequal(p_)
                push!(verts, [(pspans[i], x...) for x in unique(p_)])
            end
        end
    end

    return [collect.(zip(pspans, p)) for p in paths], verts
end

function plot_comparison(
        simulation::Vector{CellSimulationResults},
        agentbased::Vector{Vector{AnalyticalResults}},
        ddilution::EffectiveDilutionModel,
        sdilution::Vector{StochasticDilutionModel},
        experiment;
        bifurcation_parameter="α", title="", xlims=300, ylims, tlims, tylims, divtidx, xmlims, idxm, minprom=1e-5,
        xticks,
        yticks,
        yticklabels,
        offset,
        ridgelims, 
        acolors,
        nlins,
        trtlimx,
        trtlimy,
        sidx,
        heightmod = 1.0,
        linestyle = :solid
    )

    # Plotting.
    size_inches = (17, 6 * heightmod)
    size_pt = 72 .* size_inches

    fig = Figure(resolution=size_pt, fontsize=16, figure_padding=1)
    sub_top = GridLayout(fig[1,:])
    axsim = Axis(sub_top[1,1]; xlabel="Time", ylabel="Protein count", title="Simulated lineages")
    axmode = Axis(sub_top[1,2]; xlabel="Parameter α", ylabel="Protein count", title="Marginal protein distribution modes")
    axmarginals = Axis3(sub_top[1,3]; xlabel="Protein count", title="Marginal protein distributions", 
        azimuth=1.6pi, 
        elevation=0.08pi,
        xspinecolor_2=:transparent,
        xspinecolor_3=:transparent,
        yspinecolor_1=:transparent,
        yspinecolor_3=:transparent,
        zspinecolor_1=:transparent,
        zspinecolor_2=:transparent,
        aspect=(1,1,1/3),
        xticks=xticks,
        yticks=(yticks, ["α = $a" for a in yticklabels]),
        xlabelrotation=0.0,
        ylabel="",
        protrusions=20,
        zlabel="Probability density (arbitrary units)",
        yticklabelpad=-45,
        viewmode=:stretch,
       )

    for (sim, colr, n) in zip(simulation, acolors, nlins)
        idxs = rand(collect(1:length(sim.results[:final_population_traj])), n) 
        cidxs = findall(
            x -> x.idx in cell_index.(sim.results[:final_population_traj][idxs]), 
            sim.results[:trajectories])
        lins_ = lineage.(Ref(sim.results[:trajectories]), cidxs)
        lins = [state_time.(l) for l in lins_]

        for lin in lins
            ts = vcat([traj.t for traj in reverse(lin)]...)
            ys = vcat([getindex.(traj.u, sidx) for traj in reverse(lin)]...)
            lines!(axsim, ts, ys; color=(colors[colr], 0.2), linewidth=1.0)
        end

        linlast = last(lins)
        ts = vcat([traj.t for traj in reverse(linlast)]...)
        ys = vcat([getindex.(traj.u, sidx) for traj in reverse(linlast)]...)
        lines!(axsim, ts, ys; color=colors[colr], linewidth=1.0, linestyle=linestyle)

        divs = [(traj.t[end], traj.u[end][sidx])  for traj in reverse(linlast)]
        scatter!(axsim, getindex.(divs, 1), getindex.(divs, 2); color=(colors[colr], 0.6), markersize=10)
    end

    xlims!(axsim, trtlimx)
    ylims!(axsim, trtlimy)

    # Deterministic dilution bifurcations.
    roots_ = [(x,y[sidx]) for (x,y) in vcat(ddilution.roots...)]
    lines!(axmode, sort(roots_, lt=mylt); 
        color=(:black, 0.8), 
        label="Effective dilution (deterministic)", 
        overdraw=true, linewidth=3, linestyle=:dash)

    cmaps = [ColorScheme(range(alphacolor(colors[i], 0.2), stop=ARGB(colors[i], 0.3), length=length(idxm)+1)) for i in 1:length(colors)]

    for (agent, colr) in zip(agentbased, acolors)
        pspan = [analytical.problem.ps[ddilution.bif_idx] for analytical in agent]
        # Modes
        dists = [vcat(sum(a.results[:marginal_size]', dims=sidx)...) for a in agent]
        mode_idxs = [findmaxima(dist)[1] for dist in dists]
        mode_idxs = [peakproms(idx, dist; minprom=minprom)[1] for (dist, idx) in zip(dists, mode_idxs)]
    
        mode_idxs = map(
            x -> 
            !isempty(x[1]) ? x[1] : findmax(x[2])[2], 
            zip(mode_idxs, dists)) 
    
        modes = [collect(zip(idxs, dist[idxs])) for (idxs, dist) in zip(mode_idxs, dists)]
        modes, vert = process_modes(modes, pspan)

        for mode in modes
            scatter!(axmode, Point2f.([[m[1], m[2][1]] for m in mode]); color=colors[colr], label="Agent based", markersize=10)
            lines!(axmode, Point2f.([[m[1], m[2][1]] for m in mode]); 
                color=colors[colr], 
                colormap=range(to_color((colors[colr], 0.3)), stop=to_color(colors[colr]), length=5),
                overdraw=true, linewidth=3, linestyle=linestyle)

        end

    end
    xoff_ = offset[1] * (length(idxm)+1)
    yoff_ = offset[2] * (length(idxm)+1)
    
    for (i, idx) in enumerate(reverse(idxm))
        xoff = offset[1] * (length(idxm)-i)
        yoff = offset[2] * (length(idxm)-i)
    
        for (agent, colr) in zip(agentbased, acolors)
            dist = sum(agent[idx].results[:marginal_size]', dims=sidx)
            len_ = length(dist)
            lenn = length(dist[1:min(len_, xmlims[end])])
            valsn = dist[1:min(len_, xmlims[end])]
            lowern = Point3f.(collect(1:lenn)[1:min(len_, xmlims[end])] .- 1, fill(yoff, lenn), 0.0)
            uppern = Point3f.(collect(1:lenn)[1:min(len_, xmlims[end])] .- 1, fill(yoff, lenn), valsn)
    
            l1 = lines!(axmarginals, uppern; color=colors[colr], linewidth=3, linestyle=linestyle)
            if i != 1
                lines!(axmarginals, lowern; color=(:black, 0.8), linewidth=1, linestyle=linestyle)
            end
    
            bnd1 = band!(axmarginals,
                lowern,
                uppern; color=valsn, colormap=cmaps[colr], rasterize = 5)
        end

        # Stochastic dilution
        sdil_ss = sum(sdilution[idx].steady_state', dims=sidx) 
        len_ = length(sdil_ss)
        lenn = length(sdil_ss[1:min(len_, xmlims[end])])
        valsn = sdil_ss[1:min(len_, xmlims[end])]

        lowern = Point3f.(collect(1:lenn)[1:min(len_, xmlims[end])] .- 1, fill(yoff, lenn), 0.0)
        uppern = Point3f.(collect(1:lenn)[1:min(len_, xmlims[end])] .- 1, fill(yoff, lenn), valsn)
    
        l1 = lines!(axmarginals, uppern; color=colors[2], linewidth=3, linestyle=linestyle)
        if i != length(idxm)
            lines!(axmarginals, lowern; color=(:black, 0.8), linewidth=1, linestyle=linestyle)
        end
    
        bnd1 = band!(axmarginals,
            lowern,
            uppern; color=valsn, colormap=cmaps[2], rasterize = 5)
    end

    # Stochastic dilution.
    pspansdil = getindex.(getfield.(sdilution, :ps), 1)
    
    # Modes
    sddist = [vec(sum(sd.steady_state', dims=sidx)) for sd in sdilution]
    mode_idxs_sd = [findmaxima(sd)[1] for sd in sddist]
    mode_idxs_sd = [peakproms(idx, dist; minprom=minprom)[1] for (dist, idx) in zip(sddist, mode_idxs_sd)]
    mode_idxs_sd = map(
        x -> !isempty(x[1]) ? x[1] : findmax(x[2])[2], 
        zip(mode_idxs_sd, sddist)) 

    modes_sd = [collect(zip(idxs, sd[idxs])) for (idxs, sd) in zip(mode_idxs_sd, sddist)]
    modes_sd, vert_sd = process_modes(modes_sd, pspansdil)

    for mode in modes_sd
        scatter!(axmode, Point2f.([[m[1], m[2][1]] for m in mode]); color=colors[2], label="Stochastic dilution", markersize=10)
        lines!(axmode, Point2f.([[m[1], m[2][1]] for m in mode]); 
               color=colors[2], 
               linewidth=3)
    end

    if !isnothing(vert_sd)
        for v in vert_sd
            pts = [Point2f(v_[1], v_[2]) for v_ in v]
            cls = getindex.(v, 3)

            lines!(axmode, pts; 
                linestyle=:dot,
                linewidth=3,
                color=(colors[2], 0.3))
        end
    end

    # Heatmap
    mdists = hcat(paddedviews(0.0, [vec(sum(d.steady_state', dims=sidx)) for d in sdilution]...)...)'
    len = length(mdists[1,:])
    wd = diff(pspansdil) ./ 2
        
    for (i, p) in enumerate(pspansdil)
        xs = range(p-wd[1], stop=p+wd[1], length=9)
        ys = collect(0:len-1)
        zs = [mdists[i,j] for j in 1:length(ys), k in xs[2:end-1]]
        heatmap!(axmode, xs[2:end-1], ys, zs'; colormap=range(to_color((colors[2], 0.0)), stop=to_color((colors[2], 0.5)), length=20))
    end
    ylims!(axmode, ylims) 

    ylims!(axmarginals, low=0.0)
    zlims!(axmarginals, (0.0, ridgelims[3]))
    xlims!(axmarginals, (0.0, ridgelims[1]))
    ylims!(axmarginals, (0.0, ridgelims[2]))

    hidedecorations!(axmode, ticks=false, label=false, ticklabels=false)
    hidespines!(axmode, :r, :t)

    hideydecorations!(axmarginals, ticklabels=false, label=false)
    hidexdecorations!(axmarginals, ticklabels=false, ticks=false, label=false)
    hidezdecorations!(axmarginals)

    hidedecorations!(axsim, ticks=false, label=false, ticklabels=false)
    hidespines!(axsim, :r, :t)


    effdil_colors = [LineElement(color=(:black, 0.8), linestyle=:dash, linewidth=3), 
                     LineElement(color=(colors[2]), linewidth=3)]
    agent_colors = [LineElement(color=colors[i], linewidth=3) for i in acolors]
    
    dist_colors = [LineElement(color=colors[2], linewidth=3),
                   LineElement(color=colors[1], linewidth=3), 
                   LineElement(color=colors[4], linewidth=3)]

    axislegend(axmode, [agent_colors, effdil_colors], [["Population modes", "Mother machine modes" ], ["Deterministic", "Stochastic modes"]], ["Agent based", "Effective dilution"]; orientation=:vertical, framevisible=false, tellwidth=false, tellheight=true, groupgap=10, position=:lt, nbanks=1, halign=:left,titlehalign=:left, gridshalign=:left)

    axislegend(axmarginals, dist_colors, ["Effective dilution", "Population modes", "Mother machine modes" ]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true, groupgap=10, position=:rt, nbanks=3)

    axislegend(axsim, agent_colors, ["Population histories", "Mother machine lineages" ]; orientation=:vertical, framevisible=false, tellwidth=false, tellheight=true, groupgap=10, position=:lt, nbanks=1, halign=:left,titlehalign=:left, gridshalign=:left)


    sname = savename(experiment; allowedtypes = (Tuple, Vector, String, Int, Float64), ignores=["iters", "max_pop", "simulation_tspan", "init", "jitt", "Δt"])
    prefix = DrWatson.default_prefix(experiment)
    mkpath("$(plotsdir())/$prefix")
    save("$(plotsdir())/$prefix/$(sname)_effective_dilution.pdf", fig)
    save("$(plotsdir())/$prefix/$(sname)_effective_dilution.svg", fig)
end

function plot_comparison_extras(
        agentbased::Vector{Vector{AnalyticalResults}},
        ddilution::EffectiveDilutionModel,
        sdilution::Vector{StochasticDilutionModel},
        experiment;
        bifurcation_parameter="α", title="", xlims=300, ylims, tlims, tylims, divtidx, xmlims, idxm, minprom=1e-4,
        tticks,
        xticks,
        yticks,
        yticklabels,
        offset,
        ridgelims, 
        tridgelims, 
        acolors,
        sidx = 1
        )

    # Plotting.
    size_inches = (17, 12)
    size_pt = 72 .* size_inches

    fig = Figure(resolution=size_pt, fontsize=16, figure_padding=1)
    axinterdiv = Axis3(fig[1,1]; xlabel="Interdivision time", title="Interdivision time distributions", 
        azimuth=1.6pi, 
        elevation=0.08pi,
        xspinecolor_2=:transparent,
        xspinecolor_3=:transparent,
        yspinecolor_1=:transparent,
        yspinecolor_3=:transparent,
        zspinecolor_1=:transparent,
        zspinecolor_2=:transparent,
        aspect=(1,1,1/3),
        xticks=tticks,
        yticks=(yticks, ["α = $a" for a in yticklabels]),
        xlabelrotation=0.0,
        ylabel="",
        protrusions=35,
        zlabel="Probability density (arbitrary units)",
        yticklabelpad=-45,
        viewmode=:stretch,
       )

    axbirth = Axis3(fig[1,2]; xlabel="Protein count", title="Birth protein distributions", 
        azimuth=1.6pi, 
        elevation=0.08pi,
        xspinecolor_2=:transparent,
        xspinecolor_3=:transparent,
        yspinecolor_1=:transparent,
        yspinecolor_3=:transparent,
        zspinecolor_1=:transparent,
        zspinecolor_2=:transparent,
        aspect=(1,1,1/3),
        xticks=xticks,
        yticks=(yticks, ["α = $a" for a in yticklabels]),
        xlabelrotation=0.0,
        ylabel="",
        protrusions=35,
        zlabel="Probability density (arbitrary units)",
        yticklabelpad=-45,
        viewmode=:stretch,
       )

    axdiv = Axis3(fig[2,1]; xlabel="Protein count", title="Division protein distributions", 
        azimuth=1.6pi, 
        elevation=0.08pi,
        xspinecolor_2=:transparent,
        xspinecolor_3=:transparent,
        yspinecolor_1=:transparent,
        yspinecolor_3=:transparent,
        zspinecolor_1=:transparent,
        zspinecolor_2=:transparent,
        aspect=(1,1,1/3),
        xticks=xticks,
        yticks=(yticks, ["α = $a" for a in yticklabels]),
        xlabelrotation=0.0,
        ylabel="",
        protrusions=35,
        zlabel="Probability density (arbitrary units)",
        yticklabelpad=-45,
        viewmode=:stretch,
       )

    cmaps = [ColorScheme(range(alphacolor(colors[i], 0.1), stop=ARGB(colors[i], 0.1), length=length(idxm)+1)) for i in 1:4]

    for (i, idx) in enumerate(reverse(idxm))
        xoff = offset[1] * (length(idxm)-i)
        yoff = offset[2] * (length(idxm)-i)
        
        for (agent, colr) in zip(agentbased, acolors)
            # Birth distribution
            bdist = sum(agent[idx].results[:birth_dist]', dims=sidx)

            #lenb = length(agent[idx].results[:birth_dist][1:xmlims[end]])
            lenb = length(bdist[1:xmlims[end]])
#            valsb = agent[idx].results[:birth_dist][1:xmlims[end]]
            valsb = bdist[1:xmlims[end]]
            lowerb = Point3f.(collect(1:lenb)[1:xmlims[end]] .- 1, fill(yoff, lenb), 0.0)
            upperb = Point3f.(collect(1:lenb)[1:xmlims[end]] .- 1, fill(yoff, lenb), valsb)
    
            lines!(axbirth, upperb; color=colors[colr], linewidth=3)
            if i != length(idxm)
                lines!(axbirth, lowerb; color=:black, linewidth=1)
            end
    
            band!(axbirth,
                lowerb,
                upperb; color=valsb, colormap=cmaps[colr], rasterize = 5)

            # Division distribution
            ddist = sum(agent[idx].results[:division_dist_ancest]', dims=sidx)
#            lend = length(agent[idx].results[:division_dist_ancest][1:xmlims[end]])
            lend = length(ddist[1:xmlims[end]])

            valsd = ddist[1:xmlims[end]]
#            valsd = agent[idx].results[:division_dist_ancest][1:xmlims[end]]
            lowerd = Point3f.(collect(1:lend)[1:xmlims[end]] .- 1, fill(yoff, lend), 0.0)
            upperd = Point3f.(collect(1:lend)[1:xmlims[end]] .- 1, fill(yoff, lend), valsd)
    
            lines!(axdiv, upperd; color=colors[colr], linewidth=3)
            if i != length(idxm)
                lines!(axdiv, lowerd; color=:black, linewidth=1)
            end
    
            band!(axdiv,
                lowerd,
                upperd; color=valsd, colormap=cmaps[colr], rasterize = 5)
 
            # Interdivision time distribution
            ts = collect(tlims[1]:0.001:tlims[2]) 
            lent = length(ts)
            valst = division_time_ancest(agent[idx]).(ts)
            lowert = Point3f.(ts, fill(yoff, lent), 0.0)
            uppert = Point3f.(ts, fill(yoff, lent), valst)

            band!(axinterdiv,
                lowert,
                uppert; color=valst, colormap=cmaps[colr], rasterize = 5)

            lines!(axinterdiv, uppert; color=colors[colr], linewidth=3)

            if i != length(idxm)
                lines!(axinterdiv, lowert; color=:black, linewidth=1)
            end
        end

        # Stochastic dilution
        sdist = sum(sdilution[idx].steady_state', dims=sidx)
        lenn = length(sdist[1:xmlims[end]])
#        lenn = length(sdilution[idx].steady_state[1:xmlims[end]])
#        valsn = sdilution[idx].steady_state[1:xmlims[end]]
        valsn = sdist[1:xmlims[end]]

        lowern = Point3f.(collect(1:lenn)[1:xmlims[end]] .- 1, fill(yoff, lenn), 0.0)
        uppern = Point3f.(collect(1:lenn)[1:xmlims[end]] .- 1, fill(yoff, lenn), valsn)
    
        lines!(axbirth, uppern; color=colors[2], linewidth=3)
        lines!(axdiv, uppern; color=colors[2], linewidth=3)
        if i != 1
            lines!(axbirth, lowern; color=:black, linewidth=1)
        end
    
        band!(axbirth,
            lowern,
            uppern; color=valsn, colormap=cmaps[2], rasterize = 5)
        band!(axdiv,
            lowern,
            uppern; color=valsn, colormap=cmaps[2], rasterize = 5)
    end

    zlims!(axbirth, (0.0, ridgelims[3]))
    xlims!(axbirth, (0.0, ridgelims[1]))
    ylims!(axbirth, (0.0, ridgelims[2]))

    hideydecorations!(axbirth, ticklabels=false, label=false)
    hidexdecorations!(axbirth, ticklabels=false, ticks=false, label=false)
    hidezdecorations!(axbirth)

    hideydecorations!(axdiv, ticklabels=false, label=false)
    hidexdecorations!(axdiv, ticklabels=false, ticks=false, label=false)
    hidezdecorations!(axdiv)
    zlims!(axdiv, (0.0, ridgelims[3]))
    xlims!(axdiv, (0.0, ridgelims[1]))
    ylims!(axdiv, (0.0, ridgelims[2]))

    ylims!(axinterdiv, low=0.0)
    xlims!(axinterdiv, tlims)
    zlims!(axinterdiv, (0.0, tridgelims[2]))

    hideydecorations!(axinterdiv, ticklabels=false, label=false)
    hidexdecorations!(axinterdiv, ticklabels=false, ticks=false, label=false)
    hidezdecorations!(axinterdiv)

    effdil_colors = [LineElement(color=colors[4], linestyle=:dot, linewidth=3), 
                     LineElement(color=(colors[2]), linewidth=3)]
    agent_colors = [LineElement(color=colors[1], linewidth=3), 
                    LineElement(color=colors[4], linewidth=3)]
    
    dist_colors = [LineElement(color=colors[2], linewidth=3),
                   LineElement(color=colors[1], linewidth=3), 
                   LineElement(color=colors[4], linewidth=3)]

    axislegend(axdiv, dist_colors, ["Effective dilution", "Population modes", "Mother machine modes" ]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true, groupgap=10, position=:rt, nbanks=3)

    axislegend(axbirth, dist_colors, ["Effective dilution", "Population modes", "Mother machine modes" ]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true, groupgap=10, position=:rt, nbanks=3)

    axislegend(axinterdiv, dist_colors[2:end], ["Population modes", "Mother machine modes" ]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true, groupgap=10, position=:rt, nbanks=3)

    sname = savename(experiment; allowedtypes = (Tuple, Vector, String, Int, Float64), ignores=["iters", "max_pop", "simulation_tspan", "init", "jitt", "Δt"])
    prefix = DrWatson.default_prefix(experiment)
    mkpath("$(plotsdir())/$prefix")
    save("$(plotsdir())/$prefix/$(sname)_effective_dilution_extras.pdf", fig)
    save("$(plotsdir())/$prefix/$(sname)_effective_dilution_extras.svg", fig)
#    save("$(plotsdir())/$prefix/$(sname)_effective_dilution_extras.eps", fig)
end


function plot_convergence(analytical, truncation_error, time_trunc_error, experiment; idxs, truncations, label)
    size_inches = (17, 7)
    size_pt = 72 .* size_inches

    fig = Figure(resolution=size_pt, fontsize=19)
    ax_birth = Axis(fig[1,1]; xlabel="Protein counts", title="Birth distribution convergence -- $label", 
                    yticks = ([0.0, 0.08], ["Selection", "No selection"]))

    ax_growth = Axis(fig[1,1]; 
                        width=Relative(0.25),
                        height=Relative(0.25),
                        halign=0.97,
                        valign=0.97,
                        xlabelsize=16,
                        xticklabelsize=16,
                        xticksize=2,
                        ylabelsize=16,
                        yticklabelsize=16,
                        yticksize=2,
                        xlabel="Iteration", 
                        ylabel="Growth rate λ")
    
    ax_error3d = Axis3(fig[1,2]; 
                        xlabel=L"\tau_{\text{max}}", 
                        ylabel="Truncation size", 
                        zlabel="Exit probability", 
                        title="Exit probability -- $label", 
                        viewmode=:stretch, 
                        azimuth=2.2pi, 
                        elevation=0.08pi, 
                        protrusions=50)


    for (i,a) in enumerate(analytical)
        len = length(a.convergence_monitor.growth_factor)
        
        lines!(ax_growth, collect(2:len) .- 1, a.convergence_monitor.growth_factor[2:end]; color=colors[i])

        toplot = a.convergence_monitor.birth_dists[idxs[i]]
        cmap = range(alphacolor(colors[i], 0.1), stop=ARGB(colors[i], 0.6), length=length(toplot))

        for (j,dist) in enumerate(toplot)
            len_birth = length(sum(dist; dims=1))
            band!(ax_birth, 
                collect(0:len_birth-1), 
                fill(0.0, len_birth) .+ (length(analytical)-i) * 0.08,
                vec(sum(dist; dims=1)) .+ (length(analytical)-i) * 0.08; 
                color = cmap[j])

            lines!(ax_birth, collect(0:length(sum(dist; dims=1))-1), vec(sum(dist; dims=1)) .+ (length(analytical)-i) * 0.08; color=cmap[j])
        end
    end
    ylims!(ax_growth, (0.5, 1.2))
    xlims!(ax_growth, (1, 10))
    hidedecorations!(ax_growth, ticks=false, label=false, ticklabels=false)
    hidespines!(ax_growth, :r, :t)

    hidedecorations!(ax_birth, ticks=false, label=false, ticklabels=false)
    hidespines!(ax_birth, :r, :t)
    xlims!(ax_birth, (0, 50))
    ylims!(ax_birth, low=0.0)


    hidedecorations!(ax_error3d, ticks=false, label=false, ticklabels=false)

    legend_colors = [LineElement(color=colors[1], linewidth=3), 
                     LineElement(color=colors[2], linewidth=3)]
    
    axislegend(ax_error3d, legend_colors, ["$label - no selection", "$label - selection"]; framevisible=false)

    for (i, errors) in enumerate(time_trunc_error)
        cmap = ColorScheme(range(alphacolor(colors[i], 0.1), stop=ARGB(colors[i], 1.0), length=100))
        times = [e.problem.tspan[end] for e in time_trunc_error[i]][1,:]
        truncs = [e.problem.approx.truncation[end]  for e in time_trunc_error[i]][:,1]
        errs = [last(e.results[:errors])  for e in time_trunc_error[i]]
        surface!(ax_error3d, times, Float64.(truncs), getindex.(errs, 1)'; colormap=cmap, rasterize=10)
    end

    xlims!(ax_error3d, (0.0, 7.))
    ylims!(ax_error3d, (20, 100))
    zlims!(ax_error3d, low=0.0)

    prefix = DrWatson.default_prefix(experiment)
    mkpath("$(plotsdir())/$prefix")
    save("$(plotsdir())/$prefix/convergence.pdf", fig)
    save("$(plotsdir())/$prefix/convergence.svg", fig)
end

function plot_convergence_extras(analytical, truncation_error, truncation_error_sh, time_trunc_error, experiment; simulations, idxs, acolors, note="", drawgrowth=false, mlabel, labels, mother=false)
    size_inches = (17, 21)
    size_pt = 72 .* size_inches

    fig = Figure(resolution=size_pt, fontsize=19)
    sub_birth = GridLayout(fig[1,1]; title="Birth distribution convergence")
    ax_birth = [Axis(sub_birth[i,1]; xlabel="Protein count at birth", ylabel="Probability density") for (i,a) in  enumerate(analytical)]

    sub_div = GridLayout(fig[2,1]; )
    ax_div = [Axis(sub_div[i,1]; xlabel="Protein count at division", ylabel="Probability density") for (i,a) in  enumerate(analytical)]

    ax_growth = Axis(fig[1,1]; 
                        width=Relative(0.25),
                        height=Relative(0.25),
                        halign=0.97,
                        valign=0.97,
                        xlabelsize=16,
                        xticklabelsize=16,
                        xticksize=2,
                        ylabelsize=16,
                        yticklabelsize=16,
                        yticksize=2,
                        xlabel="Iteration", 
                        ylabel="Growth rate λ")
    ylims!(ax_growth, (0.5, 1.2))
    xlims!(ax_growth, (1, 10))
    hidedecorations!(ax_growth, ticks=false, label=false, ticklabels=false)
    hidespines!(ax_growth, :r, :t)

    ax_error3d = Axis3(fig[1,2]; 
                        xlabel=L"\tau_{\text{max}}", 
                        ylabel="Truncation size", 
                        zlabel="Exit probability", 
                        title="Exit probability — $mlabel", 
                        viewmode=:stretch, 
                        azimuth=2.2pi, 
                        elevation=0.08pi, 
                        protrusions=50)


    ax_error = Axis(fig[3,1]; xlabel="Truncation size", ylabel="Exit probability", title="Exit probability — $mlabel")
    ax_error2 = Axis(fig[3,2]; xlabel="Truncation size", ylabel="Exit probability", title="Exit probability — $mlabel")

    i = 1
    for (a, colr) in zip(analytical, acolors)
        len = length(a.convergence_monitor.monitor[:growth_rate])
        
        lines!(ax_growth, collect(2:len) .- 1, a.convergence_monitor.monitor[:growth_rate][2:end]; color=colors[colr])

        toplot_birth = a.convergence_monitor.monitor[:birth_dist][idxs[i]]
        toplot_div = a.convergence_monitor.monitor[:division_dist][idxs[i]]
        cmap = range(alphacolor(colors[colr], 0.1), stop=ARGB(colors[colr], 0.6), length=length(toplot_birth))

        for (j,dist) in enumerate(toplot_birth)
            len_birth = length(sum(dist; dims=1))
            band!(ax_birth[i], 
                collect(0:len_birth-1), 
                fill(0.0, len_birth),
                vec(sum(dist; dims=1)); 
                color = cmap[j])

            lines!(ax_birth[i], collect(0:length(sum(dist; dims=1))-1), vec(sum(dist; dims=1)); color=cmap[j])
        end

        for (j,dist) in enumerate(toplot_div)
            len_div = length(sum(dist; dims=1))
            band!(ax_div[i], 
                collect(0:len_div-1), 
                fill(0.0, len_div),
                vec(sum(dist; dims=1)); 
                color = cmap[j])

            lines!(ax_div[i], collect(0:length(sum(dist; dims=1))-1), vec(sum(dist; dims=1)); color=cmap[j])
        end

        bar_birth_all = make_normalised_bar_data(getindex.(simulations[i].results[:molecules_at_birth_all], 2))
        sort!(bar_birth_all; by=x->x[1])
        stairs!(ax_birth[i], getindex.(bar_birth_all, 1), getindex.(bar_birth_all, 2); 
            color=colors[colr], step=:center, linewidth=5, linestyle=:dot)

        bar_div_all = make_normalised_bar_data(getindex.(simulations[i].results[:molecules_at_division_all], 2))
        sort!(bar_div_all; by=x->x[1])
        stairs!(ax_div[i], getindex.(bar_div_all, 1), getindex.(bar_div_all, 2); 
            color=colors[colr], step=:center, linewidth=5, linestyle=:dot)
        display("hi")

        i += 1 
    end

    for ax in ax_birth
        hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
        hidespines!(ax, :r, :t)
        xlims!(ax, (0, 70))
        ylims!(ax, low=0.0, high=0.08)
    end
    hidexdecorations!(ax_birth[1])
    rowgap!(sub_birth, 5)
    ax_birth[1].title = "Birth protein distribution convergence"

    for ax in ax_div
        hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
        hidespines!(ax, :r, :t)
        xlims!(ax, (0, 70))
        ylims!(ax, low=0.0, high=0.04)
    end
    hidexdecorations!(ax_div[1])
    rowgap!(sub_div, 5)
    ax_div[1].title = "Division protein distribution convergence"

    hidedecorations!(ax_error3d, ticks=false, label=false, ticklabels=false)
    legend_colors = [PolyElement(color=colors[i], linewidth=3) for i in acolors]
    axislegend(ax_error3d, legend_colors, labels; framevisible=true)

    analytical_colors = [PolyElement(color=colors[colr]) for colr in acolors]
    analytical_colors = [analytical_colors..., LineElement(color=:grey, linestyle=:dot, linewidth=5)]
    axislegend(ax_div[2], analytical_colors, ["No selection", "Selection", "Simulation"]; orientation=:vertical, framevisible=false, tellwidth=false, tellheight=true)
    axislegend(ax_birth[2], analytical_colors, ["No selection", "Selection", "Simulation"]; orientation=:vertical, framevisible=false, tellwidth=false, tellheight=true)


    for (error, colr) in zip(truncation_error, acolors)
        errors = [last(e.results[:errors]) for e in error] 
        lines!(ax_error, truncations, getindex.(errors, 2); color=colors[colr], linewidth=5)
    end

    for (error, colr) in zip(truncation_error_sh, acolors)
        errors = [last(e.results[:errors]) for e in error] 
        lines!(ax_error2, truncations, getindex.(errors, 2); color=colors[colr], linewidth=3)
    end

    i = 1
    for (errors, colr) in zip(time_trunc_error, [acolors..., acolors...])
        cmap = ColorScheme(range(alphacolor(colors[colr], 0.1), stop=ARGB(colors[colr], 1.0), length=100))
        times = [e.problem.tspan[end] for e in time_trunc_error[i]][1,:]
        truncs = [e.problem.approx.truncation[end]  for e in time_trunc_error[i]][:,1]
        errs = [last(e.results[:errors]) for e in time_trunc_error[i]]
        surface!(ax_error3d, times, Float64.(truncs), getindex.(errs, 1)'; colormap=cmap, rasterize=10)
        i += 1 
    end

    xlims!(ax_error3d, (0.0, 7.))
    ylims!(ax_error3d, (20, 100))
    zlims!(ax_error3d, low=0.0)

    ylims!(ax_error, (0.0, 1.0))
    xlims!(ax_error, (20, 100))
    hidedecorations!(ax_error, ticks=false, label=false, ticklabels=false)
    hidespines!(ax_error, :r, :t)

    ylims!(ax_error2, (0.0, 1.0))
    xlims!(ax_error2, (20, 100))
    hidedecorations!(ax_error2, ticks=false, label=false, ticklabels=false)
    hidespines!(ax_error2, :r, :t)

    analytical_colors = [PolyElement(color=colors[colr]) for colr in acolors]
    axislegend(ax_error, analytical_colors, ["No selection", "Selection"]; orientation=:vertical, framevisible=false, tellwidth=false, tellheight=true)
    axislegend(ax_error2, analytical_colors, ["No selection", "Selection"]; orientation=:vertical, framevisible=false, tellwidth=false, tellheight=true)
    colsize!(fig.layout, 1, Relative(1.3/3))

    prefix = DrWatson.default_prefix(experiment)
    mkpath("$(plotsdir())/$prefix")
    save("$(plotsdir())/$prefix/convergence_$note.pdf", fig)
    save("$(plotsdir())/$prefix/convergence_$note.svg", fig)
end

function plot_simulation_traj(
        simulations::Vector{CellSimulationResults}, 
        simulations_mother::Vector{CellSimulationResults}, 
        experiment;
        scolors,
        nlins,
        tlims,
        ylims,
        initn,
        sidx=1
    )

    size_inches = (5.5*length(simulations), 6)
    size_pt = 72 .* size_inches

    fig = Figure(resolution=size_pt, fontsize=16, figure_padding=(1,15,1,1))
    axsim = [Axis(fig[1,i]; xlabel="Time", ylabel="Protein count", title="Simulated lineages")
             for i in 1:length(simulations)] 

    # Population
    for (sim, ax) in zip(simulations, axsim)
        idxs = rand(collect(1:length(sim.results[:final_population_traj])), nlins[1]) 
        cidxs = findall(
            x -> x.idx in cell_index.(sim.results[:final_population_traj][idxs]), 
            sim.results[:trajectories])
        lins_ = lineage.(Ref(sim.results[:trajectories]), cidxs)
        lins = [state_time.(l) for l in lins_]

        for lin in lins
            ts = vcat([traj.t for traj in reverse(lin)]...)
            ys = vcat([getindex.(traj.u, sidx) for traj in reverse(lin)]...)
            lines!(ax, ts, ys; color=(colors[scolors[1]], 0.2), linewidth=1.0)
        end

        linlast = last(lins)
        ts = vcat([traj.t for traj in reverse(linlast)]...)
        ys = vcat([getindex.(traj.u, sidx) for traj in reverse(linlast)]...)

        lines!(ax, ts, ys; color=colors[scolors[1]], linewidth=1.0)

        divs = [(traj.t[end], traj.u[end][sidx])  for traj in reverse(linlast)]
        scatter!(ax, getindex.(divs, 1), getindex.(divs, 2); color=(colors[scolors[1]], 0.6), markersize=10)
    end

    # Mother machine
    for (sim, ax) in zip(simulations_mother, axsim)
        idxs = rand(collect(1:length(sim.results[:final_population_traj])), nlins[2]) 
        cidxs = findall(
            x -> x.idx in cell_index.(sim.results[:final_population_traj][idxs]), 
            sim.results[:trajectories])
        lins_ = lineage.(Ref(sim.results[:trajectories]), cidxs)
        lins = [state_time.(l) for l in lins_]

        for lin in lins

            ts = vcat([traj.t for traj in reverse(lin)]...)
            ys = vcat([getindex.(traj.u, sidx) for traj in reverse(lin)]...)
            lines!(ax, ts, ys; color=(colors[scolors[2]], 0.2), linewidth=1.0)
        end

        linlast = last(lins)
        ts = vcat([traj.t for traj in reverse(linlast)]...)
        ys = vcat([getindex.(traj.u, sidx) for traj in reverse(linlast)]...)
        lines!(ax, ts, ys; color=colors[scolors[2]], linewidth=1.0)

        divs = [(traj.t[end], traj.u[end][sidx])  for traj in reverse(linlast)]
        scatter!(ax, getindex.(divs, 1), getindex.(divs, 2); color=(colors[scolors[2]], 0.6), markersize=10)
    end

    for (i, ax) in zip(initn, axsim)
        scatter!(ax, 0.0, getindex(i, sidx); color=:black, markersize=15)
    end

    xlims!.(axsim, Ref(tlims))
    ylims!.(axsim, Ref(ylims))

    hidedecorations!.(axsim, ticks=false, label=false, ticklabels=false)
    hidespines!.(axsim, :r, :t)

    sim_colors = [PolyElement(color=colors[1]), PolyElement(color=(colors[4]))]
    Legend(fig[2,1:length(simulations)], sim_colors, ["Lineage tree histories", "Mother machine lineages" ]; orientation=:horizontal, framevisible=false, tellwidth=false, tellheight=true)


    prefix = DrWatson.default_prefix(experiment[1])
    mkpath("$(plotsdir())/$prefix")
    save("$(plotsdir())/$prefix/trajectories.pdf", fig)
#    save("$(plotsdir())/$prefix/trajectories.svg", fig)
end
