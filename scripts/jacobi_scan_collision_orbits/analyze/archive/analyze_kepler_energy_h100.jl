using JLD2
using Printf

const DATA_DIR = joinpath(@__DIR__, "../../../data/jacobi_scan")
const EARTH_MOON_DISTANCE_KM = 384400.0
const MOON_RADIUS_KM = 1737.4
const ALTITUDE_KM = 100.0
const R_EVAL = (MOON_RADIUS_KM + ALTITUDE_KM) / EARTH_MOON_DISTANCE_KM

const METRIC_NAMES = [
    "orbit_index",
    "reached_h100",
    "crossing_index_h100",
    "kepler_energy_h100",
    "bound_h100",
]

const SUMMARY_NAMES = [
    "C",
    "n_orbits",
    "n_reached_h100",
    "n_bound_h100",
    "bound_fraction_h100",
    "mean_energy_h100",
    "min_energy_h100",
    "max_energy_h100",
]

function moon_radius(traj, mu; spatial)
    x_moon = traj[1, :] .- (1 - mu)
    y = traj[2, :]
    if spatial
        return sqrt.(x_moon .^ 2 .+ y .^ 2 .+ traj[3, :] .^ 2)
    end
    return sqrt.(x_moon .^ 2 .+ y .^ 2)
end

function crossing_state(traj, r)
    i = findfirst(>=(R_EVAL), r)
    isnothing(i) && return nothing, 0
    i == 1 && return traj[:, i], i

    a = (R_EVAL - r[i - 1]) / (r[i] - r[i - 1])
    return (1 - a) .* traj[:, i - 1] .+ a .* traj[:, i], i
end

function kepler_energy_planar(z, mu)
    r = sqrt((z[1] - 1 + mu)^2 + z[2]^2)
    v2 = z[3]^2 + (z[4] - (1 - mu))^2
    return 0.5 * v2 - mu / r
end

function kepler_energy_spatial(z, mu)
    r = sqrt((z[1] - 1 + mu)^2 + z[2]^2 + z[3]^2)
    v2 = z[4]^2 + (z[5] - (1 - mu))^2 + z[6]^2
    return 0.5 * v2 - mu / r
end

function analyze_case(path; spatial)
    data = load(path)
    all_cart = data["cart"]
    C_values = data["C_values"]
    mu = data["mu"]

    metrics = Dict{Float64, Matrix{Float64}}()
    summary = zeros(length(C_values), length(SUMMARY_NAMES))

    for (j, C) in enumerate(C_values)
        trajs = all_cart[C]
        rows = zeros(length(METRIC_NAMES), length(trajs))

        for (i, traj) in enumerate(trajs)
            r = moon_radius(traj, mu; spatial)
            z, crossing_index = crossing_state(traj, r)
            energy = isnothing(z) ? NaN : (spatial ? kepler_energy_spatial(z, mu) : kepler_energy_planar(z, mu))
            rows[:, i] .= [
                i,
                isnothing(z) ? 0.0 : 1.0,
                crossing_index,
                energy,
                isfinite(energy) && energy < 0 ? 1.0 : 0.0,
            ]
        end

        metrics[C] = rows
        energies = rows[4, isfinite.(rows[4, :])]
        n_orbits = length(trajs)
        n_bound = sum(rows[5, :])
        summary[j, :] .= [
            C,
            n_orbits,
            sum(rows[2, :]),
            n_bound,
            n_bound / n_orbits,
            sum(energies) / length(energies),
            minimum(energies),
            maximum(energies),
        ]

        @printf("C = %.2f: %.0f/%d bound at h = %.0f km\n",
                C, n_bound, n_orbits, ALTITUDE_KM)
    end

    return metrics, summary, C_values, mu
end

function write_summary_tsv(path, summary)
    open(path, "w") do io
        println(io, join(SUMMARY_NAMES, '\t'))
        for i in axes(summary, 1)
            println(io, join(summary[i, :], '\t'))
        end
    end
end

println("Planar Kepler-energy analysis")
planar_metrics, planar_summary, planar_C_values, mu = analyze_case(
    joinpath(DATA_DIR, "large_jacobi_scan_planar.jld2");
    spatial = false,
)

println("Spatial Kepler-energy analysis")
spatial_metrics, spatial_summary, spatial_C_values, _ = analyze_case(
    joinpath(DATA_DIR, "large_jacobi_scan_spatial.jld2");
    spatial = true,
)

out_path = joinpath(DATA_DIR, "kepler_energy_h100.jld2")
save(out_path,
     "altitude_km", ALTITUDE_KM,
     "moon_radius_km", MOON_RADIUS_KM,
     "earth_moon_distance_km", EARTH_MOON_DISTANCE_KM,
     "eval_radius", R_EVAL,
     "metric_names", METRIC_NAMES,
     "summary_names", SUMMARY_NAMES,
     "planar_metrics", planar_metrics,
     "planar_summary", planar_summary,
     "planar_C_values", planar_C_values,
     "spatial_metrics", spatial_metrics,
     "spatial_summary", spatial_summary,
     "spatial_C_values", spatial_C_values,
     "mu", mu)

write_summary_tsv(joinpath(DATA_DIR, "kepler_energy_h100_planar.tsv"), planar_summary)
write_summary_tsv(joinpath(DATA_DIR, "kepler_energy_h100_spatial.tsv"), spatial_summary)

println("Saved to $(relpath(out_path, joinpath(@__DIR__, "../../..")))")
