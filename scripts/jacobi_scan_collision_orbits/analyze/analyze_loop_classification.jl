using JLD2
using Printf

const DATA_DIR = joinpath(@__DIR__, "../../../data/jacobi_scan")
const R_START = 0.02
const R_END = 0.2

const METRIC_NAMES = [
    "measured",
    "azimuthal_loop",
    "radial_loop",
    "combined_loop",
]

const SUMMARY_NAMES = [
    "C",
    "n_orbits",
    "n_measured",
    "azimuthal_fraction",
    "radial_fraction",
    "combined_fraction",
]

function unwrap_angles(theta)
    out = copy(theta)
    offset = 0.0
    for i in 2:length(theta)
        dtheta = theta[i] - theta[i - 1]
        dtheta > π && (offset -= 2π)
        dtheta < -π && (offset += 2π)
        out[i] += offset
    end
    return out
end

function radii(traj, mu; spatial)
    xM = traj[1, :] .- (1 - mu)
    xE = traj[1, :] .+ mu
    y = traj[2, :]
    if spatial
        z = traj[3, :]
        return sqrt.(xM .^ 2 .+ y .^ 2 .+ z .^ 2),
               sqrt.(xE .^ 2 .+ y .^ 2 .+ z .^ 2)
    end
    return sqrt.(xM .^ 2 .+ y .^ 2), sqrt.(xE .^ 2 .+ y .^ 2)
end

function n_apogees(r)
    count = 0
    for i in 2:length(r)-1
        (r[i - 1] < r[i] && r[i] > r[i + 1]) && (count += 1)
    end
    return count
end

function classify(traj, mu; spatial)
    rM, rE = radii(traj, mu; spatial)

    i0 = findfirst(>=(R_START), rM)
    isnothing(i0) && return [0.0, 0.0, 0.0, 0.0]

    i1 = findfirst(i -> i > i0 && rM[i] >= R_END, eachindex(rM))
    isnothing(i1) && (i1 = length(rM))
    inds = i0:i1

    xM = traj[1, inds] .- (1 - mu)
    y = traj[2, inds]
    theta = unwrap_angles(atan.(y, xM))
    az_turns = maximum(abs.(theta .- theta[1])) / 2π
    apogees = n_apogees(rE[inds])

    az_loop = az_turns >= 1
    radial_loop = apogees >= 2
    combined_loop = az_loop || radial_loop

    return [
        1.0,
        az_loop ? 1.0 : 0.0,
        radial_loop ? 1.0 : 0.0,
        combined_loop ? 1.0 : 0.0,
    ]
end

function analyze(path; spatial)
    data = load(path)
    C_values = data["C_values"]
    mu = data["mu"]
    all_cart = data["cart"]

    metrics = Dict{Float64, Matrix{Float64}}()
    summary = zeros(length(C_values), length(SUMMARY_NAMES))

    for (j, C) in enumerate(C_values)
        trajs = all_cart[C]
        rows = zeros(length(METRIC_NAMES), length(trajs))

        for i in eachindex(trajs)
            rows[:, i] = classify(trajs[i], mu; spatial)
        end

        n = length(trajs)
        n_measured = sum(rows[1, :])

        summary[j, :] .= [
            C,
            n,
            n_measured,
            sum(rows[2, :]) / n,
            sum(rows[3, :]) / n,
            sum(rows[4, :]) / n,
        ]
        metrics[C] = rows

        @printf("C = %.2f: az %.3f, radial %.3f, combined %.3f\n",
                C, summary[j, 4], summary[j, 5], summary[j, 6])
    end

    return metrics, summary, C_values, mu
end

function write_tsv(path, summary)
    open(path, "w") do io
        println(io, join(SUMMARY_NAMES, '\t'))
        for i in axes(summary, 1)
            println(io, join(summary[i, :], '\t'))
        end
    end
end

println("Planar")
planar_metrics, planar_summary, planar_C, mu = analyze(
    joinpath(DATA_DIR, "large_jacobi_scan_planar.jld2");
    spatial = false,
)

println("Spatial")
spatial_metrics, spatial_summary, spatial_C, _ = analyze(
    joinpath(DATA_DIR, "large_jacobi_scan_spatial.jld2");
    spatial = true,
)

out = joinpath(DATA_DIR, "loop_classification.jld2")
save(out,
     "r_start", R_START,
     "r_end", R_END,
     "metric_names", METRIC_NAMES,
     "summary_names", SUMMARY_NAMES,
     "planar_metrics", planar_metrics,
     "spatial_metrics", spatial_metrics,
     "planar_summary", planar_summary,
     "spatial_summary", spatial_summary,
     "planar_C_values", planar_C,
     "spatial_C_values", spatial_C,
     "mu", mu)

write_tsv(joinpath(DATA_DIR, "loop_classification_planar.tsv"), planar_summary)
write_tsv(joinpath(DATA_DIR, "loop_classification_spatial.tsv"), spatial_summary)

println("Saved to $(relpath(out, joinpath(@__DIR__, "../../..")))")
