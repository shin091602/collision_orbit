using CairoMakie
using JLD2

const DATA_DIR = joinpath(@__DIR__, "../../../data/jacobi_scan")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan/analysis")
mkpath(RESULTS_DIR)

data = load(joinpath(DATA_DIR, "loop_classification.jld2"))
planar = data["planar_summary"]
spatial = data["spatial_summary"]
r_start = data["r_start"]
r_end = data["r_end"]

fig = Figure(size = (1200, 420))

ax1 = Axis(fig[1, 1], xlabel = "C", ylabel = "fraction",
           title = "Planar loop fractions")
lines!(ax1, planar[:, 1], planar[:, 4]; label = "azimuthal", linewidth = 2)
lines!(ax1, planar[:, 1], planar[:, 5]; label = "radial", linewidth = 2)
lines!(ax1, planar[:, 1], planar[:, 6]; label = "combined", linewidth = 2, linestyle = :dash)
axislegend(ax1, position = :lt)

ax2 = Axis(fig[1, 2], xlabel = "C", ylabel = "fraction",
           title = "Spatial loop fractions")
lines!(ax2, spatial[:, 1], spatial[:, 4]; label = "azimuthal", linewidth = 2)
lines!(ax2, spatial[:, 1], spatial[:, 5]; label = "radial", linewidth = 2)
lines!(ax2, spatial[:, 1], spatial[:, 6]; label = "combined", linewidth = 2, linestyle = :dash)
axislegend(ax2, position = :lt)

Label(fig[0, :], "Loop classification from r_M = $(r_start) to r_M = $(r_end)", fontsize = 22)

out = joinpath(RESULTS_DIR, "loop_classification_summary.pdf")
save(out, fig)

println("Saved to $(relpath(out, joinpath(@__DIR__, "../../..")))")
