using CairoMakie
using JLD2

const DATA_DIR = joinpath(@__DIR__, "../../../data/jacobi_scan")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan/analysis")
mkpath(RESULTS_DIR)

loop_data = load(joinpath(DATA_DIR, "loop_before_r02.jld2"))
energy_data = load(joinpath(DATA_DIR, "kepler_energy_h100.jld2"))

planar_loop = loop_data["planar_summary"]
spatial_loop = loop_data["spatial_summary"]
planar_energy = energy_data["planar_summary"]
spatial_energy = energy_data["spatial_summary"]

function plot_loop_summary()
    fig = Figure(size = (1200, 900))

    ax_fraction = Axis(fig[1, 1],
                       xlabel = "C",
                       ylabel = "loop fraction",
                       title = "Looped before r_M = 0.2")
    lines!(ax_fraction, planar_loop[:, 1], planar_loop[:, 5]; label = "planar net", linewidth = 2)
    scatter!(ax_fraction, planar_loop[:, 1], planar_loop[:, 5]; markersize = 7)
    lines!(ax_fraction, spatial_loop[:, 1], spatial_loop[:, 5]; label = "spatial net", linewidth = 2)
    scatter!(ax_fraction, spatial_loop[:, 1], spatial_loop[:, 5]; markersize = 7)
    lines!(ax_fraction, planar_loop[:, 1], planar_loop[:, 7]; label = "planar total", linewidth = 2, linestyle = :dash)
    lines!(ax_fraction, spatial_loop[:, 1], spatial_loop[:, 7]; label = "spatial total", linewidth = 2, linestyle = :dash)
    axislegend(ax_fraction, position = :lt)

    ax_count = Axis(fig[1, 2],
                    xlabel = "C",
                    ylabel = "looped orbits",
                    title = "Looped orbit counts")
    lines!(ax_count, planar_loop[:, 1], planar_loop[:, 4]; label = "planar net", linewidth = 2)
    scatter!(ax_count, planar_loop[:, 1], planar_loop[:, 4]; markersize = 7)
    lines!(ax_count, spatial_loop[:, 1], spatial_loop[:, 4]; label = "spatial net", linewidth = 2)
    scatter!(ax_count, spatial_loop[:, 1], spatial_loop[:, 4]; markersize = 7)

    ax_mean = Axis(fig[2, 1],
                   xlabel = "C",
                   ylabel = "mean turns",
                   title = "Mean xy turns before r_M = 0.2")
    lines!(ax_mean, planar_loop[:, 1], planar_loop[:, 8]; label = "planar net", linewidth = 2)
    lines!(ax_mean, spatial_loop[:, 1], spatial_loop[:, 8]; label = "spatial net", linewidth = 2)
    lines!(ax_mean, planar_loop[:, 1], planar_loop[:, 10]; label = "planar total", linewidth = 2, linestyle = :dash)
    lines!(ax_mean, spatial_loop[:, 1], spatial_loop[:, 10]; label = "spatial total", linewidth = 2, linestyle = :dash)
    axislegend(ax_mean, position = :lt)

    ax_max = Axis(fig[2, 2],
                  xlabel = "C",
                  ylabel = "max turns",
                  title = "Maximum xy turns before r_M = 0.2")
    lines!(ax_max, planar_loop[:, 1], planar_loop[:, 9]; label = "planar net", linewidth = 2)
    lines!(ax_max, spatial_loop[:, 1], spatial_loop[:, 9]; label = "spatial net", linewidth = 2)
    lines!(ax_max, planar_loop[:, 1], planar_loop[:, 11]; label = "planar total", linewidth = 2, linestyle = :dash)
    lines!(ax_max, spatial_loop[:, 1], spatial_loop[:, 11]; label = "spatial total", linewidth = 2, linestyle = :dash)

    return fig
end

function plot_energy_summary()
    fig = Figure(size = (1200, 900))

    ax_fraction = Axis(fig[1, 1],
                       xlabel = "C",
                       ylabel = "bound fraction",
                       title = "Negative Moon Kepler energy at h = 100 km")
    lines!(ax_fraction, planar_energy[:, 1], planar_energy[:, 5]; label = "planar", linewidth = 2)
    scatter!(ax_fraction, planar_energy[:, 1], planar_energy[:, 5]; markersize = 7)
    lines!(ax_fraction, spatial_energy[:, 1], spatial_energy[:, 5]; label = "spatial", linewidth = 2)
    scatter!(ax_fraction, spatial_energy[:, 1], spatial_energy[:, 5]; markersize = 7)
    axislegend(ax_fraction, position = :lt)

    ax_count = Axis(fig[1, 2],
                    xlabel = "C",
                    ylabel = "bound orbits",
                    title = "Bound orbit counts at h = 100 km")
    lines!(ax_count, planar_energy[:, 1], planar_energy[:, 4]; label = "planar", linewidth = 2)
    scatter!(ax_count, planar_energy[:, 1], planar_energy[:, 4]; markersize = 7)
    lines!(ax_count, spatial_energy[:, 1], spatial_energy[:, 4]; label = "spatial", linewidth = 2)
    scatter!(ax_count, spatial_energy[:, 1], spatial_energy[:, 4]; markersize = 7)

    ax_mean = Axis(fig[2, 1],
                   xlabel = "C",
                   ylabel = "mean energy",
                   title = "Mean Moon Kepler energy at h = 100 km")
    hlines!(ax_mean, [0.0]; color = :black, linestyle = :dot)
    lines!(ax_mean, planar_energy[:, 1], planar_energy[:, 6]; label = "planar", linewidth = 2)
    lines!(ax_mean, spatial_energy[:, 1], spatial_energy[:, 6]; label = "spatial", linewidth = 2)
    axislegend(ax_mean, position = :rt)

    ax_range = Axis(fig[2, 2],
                    xlabel = "C",
                    ylabel = "energy",
                    title = "Moon Kepler energy range at h = 100 km")
    hlines!(ax_range, [0.0]; color = :black, linestyle = :dot)
    band!(ax_range, planar_energy[:, 1], planar_energy[:, 7], planar_energy[:, 8]; color = (:dodgerblue, 0.2))
    lines!(ax_range, planar_energy[:, 1], planar_energy[:, 6]; label = "planar mean", linewidth = 2)
    band!(ax_range, spatial_energy[:, 1], spatial_energy[:, 7], spatial_energy[:, 8]; color = (:orange, 0.2))
    lines!(ax_range, spatial_energy[:, 1], spatial_energy[:, 6]; label = "spatial mean", linewidth = 2)
    axislegend(ax_range, position = :rt)

    return fig
end

function plot_combined_fraction()
    fig = Figure(size = (1000, 700))
    ax = Axis(fig[1, 1],
              xlabel = "C",
              ylabel = "fraction",
              title = "Loop and bound fractions")
    lines!(ax, planar_loop[:, 1], planar_loop[:, 5]; label = "planar loop", linewidth = 2)
    scatter!(ax, planar_loop[:, 1], planar_loop[:, 5]; markersize = 7)
    lines!(ax, spatial_loop[:, 1], spatial_loop[:, 5]; label = "spatial loop", linewidth = 2)
    scatter!(ax, spatial_loop[:, 1], spatial_loop[:, 5]; markersize = 7)
    lines!(ax, planar_energy[:, 1], planar_energy[:, 5]; label = "planar bound", linewidth = 2, linestyle = :dash)
    scatter!(ax, planar_energy[:, 1], planar_energy[:, 5]; markersize = 7)
    lines!(ax, spatial_energy[:, 1], spatial_energy[:, 5]; label = "spatial bound", linewidth = 2, linestyle = :dash)
    scatter!(ax, spatial_energy[:, 1], spatial_energy[:, 5]; markersize = 7)
    axislegend(ax, position = :lt)
    return fig
end

save(joinpath(RESULTS_DIR, "loop_before_r02_summary.pdf"), plot_loop_summary())
save(joinpath(RESULTS_DIR, "kepler_energy_h100_summary.pdf"), plot_energy_summary())
save(joinpath(RESULTS_DIR, "loop_bound_fraction_summary.pdf"), plot_combined_fraction())

println("Saved to $(relpath(RESULTS_DIR, joinpath(@__DIR__, "../../..")))")
