using CairoMakie
using JLD2

const DATA_DIR = joinpath(@__DIR__, "../../../data/jacobi_scan")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan/analysis")
mkpath(RESULTS_DIR)

data = load(joinpath(DATA_DIR, "kepler_bound_altitudes.jld2"))
alts = Int.(data["altitudes_km"])
planar = data["planar_summary"]
spatial = data["spatial_summary"]

fig = Figure(size = (1500, 620), fontsize = 28, figure_padding = (30, 30, 30, 60))

for (col, summary, title) in [(1, planar, "Planar"), (2, spatial, "Spatial")]
    ax = Axis(fig[1, col];
              xlabel = "C",
              ylabel = "bound fraction",
              title = title,
              titlesize = 54,
              xlabelsize = 42,
              ylabelsize = 42,
              xticklabelsize = 28,
              yticklabelsize = 28)
    for j in eachindex(alts)
        lines!(ax, summary[:, 1], summary[:, j + 1]; label = "$(alts[j]) km", linewidth = 3)
    end
    axislegend(ax, position = :lt, labelsize = 24)
end

Label(fig[0, :], "Negative Moon Kepler energy by altitude", fontsize = 42)

out = joinpath(RESULTS_DIR, "kepler_bound_altitudes_summary.pdf")
save(out, fig)

println("Saved to $(relpath(out, joinpath(@__DIR__, "../../..")))")
