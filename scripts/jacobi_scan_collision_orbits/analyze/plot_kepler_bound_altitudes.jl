using CairoMakie
using JLD2

const DATA_DIR = joinpath(@__DIR__, "../../../data/jacobi_scan")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan/analysis")
mkpath(RESULTS_DIR)

data = load(joinpath(DATA_DIR, "kepler_bound_altitudes.jld2"))
alts = Int.(data["altitudes_km"])
planar = data["planar_summary"]
spatial = data["spatial_summary"]

fig = Figure(size = (1200, 420))

for (col, summary, title) in [(1, planar, "Planar"), (2, spatial, "Spatial")]
    ax = Axis(fig[1, col], xlabel = "C", ylabel = "bound fraction", title = title)
    for j in eachindex(alts)
        lines!(ax, summary[:, 1], summary[:, j + 1]; label = "$(alts[j]) km", linewidth = 2)
    end
    axislegend(ax, position = :lt)
end

Label(fig[0, :], "Negative Moon Kepler energy by altitude", fontsize = 22)

out = joinpath(RESULTS_DIR, "kepler_bound_altitudes_summary.pdf")
save(out, fig)

println("Saved to $(relpath(out, joinpath(@__DIR__, "../../..")))")
