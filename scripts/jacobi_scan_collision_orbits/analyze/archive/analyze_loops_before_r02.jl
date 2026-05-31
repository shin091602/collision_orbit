using JLD2
using Printf

const DATA_DIR = joinpath(@__DIR__, "../../../data/jacobi_scan")
const R_LOOP = 0.2
const LOOP_THRESHOLD = 1.0

const METRIC_NAMES = [
    "orbit_index",
    "reached_r02",
    "exit_index_r02",
    "n_points_used",
    "net_turns_xy",
    "total_turns_xy",
    "looped_net_xy",
    "looped_total_xy",
]

const SUMMARY_NAMES = [
    "C",
    "n_orbits",
    "n_reached_r02",
    "n_looped_net_xy",
    "loop_fraction_net_xy",
    "n_looped_total_xy",
    "loop_fraction_total_xy",
    "mean_net_turns_xy",
    "max_net_turns_xy",
    "mean_total_turns_xy",
    "max_total_turns_xy",
]

function unwrap_angles(theta)
    out = copy(theta)
    offset = 0.0
    for i in 2:length(theta)
        dtheta = theta[i] - theta[i - 1]
        if dtheta > π
            offset -= 2π
        elseif dtheta < -π
            offset += 2π
        end
        out[i] += offset
    end
    return out
end

function moon_radius(traj, mu; spatial)
    x_moon = traj[1, :] .- (1 - mu)
    y = traj[2, :]
    if spatial
        return sqrt.(x_moon .^ 2 .+ y .^ 2 .+ traj[3, :] .^ 2)
    end
    return sqrt.(x_moon .^ 2 .+ y .^ 2)
end

function loop_metrics(traj, mu; spatial)
    x_moon = traj[1, :] .- (1 - mu)
    y = traj[2, :]
    r = moon_radius(traj, mu; spatial)

    exit_index = findfirst(>=(R_LOOP), r)
    last_index = isnothing(exit_index) ? length(r) : exit_index
    reached = !isnothing(exit_index)

    valid = findall(i -> isfinite(x_moon[i]) && isfinite(y[i]) && r[i] > 0, 1:last_index)
    if length(valid) < 2
        return reached, isnothing(exit_index) ? 0 : exit_index, length(valid), 0.0, 0.0
    end

    theta = unwrap_angles(atan.(y[valid], x_moon[valid]))
    net_turns = abs(theta[end] - theta[1]) / 2π
    total_turns = sum(abs.(diff(theta))) / 2π
    return reached, isnothing(exit_index) ? 0 : exit_index, length(valid), net_turns, total_turns
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
            reached, exit_index, n_points, net_turns, total_turns = loop_metrics(traj, mu; spatial)
            rows[:, i] .= [
                i,
                reached ? 1.0 : 0.0,
                exit_index,
                n_points,
                net_turns,
                total_turns,
                net_turns >= LOOP_THRESHOLD ? 1.0 : 0.0,
                total_turns >= LOOP_THRESHOLD ? 1.0 : 0.0,
            ]
        end

        metrics[C] = rows
        n_orbits = length(trajs)
        n_looped_net = sum(rows[7, :])
        n_looped_total = sum(rows[8, :])
        summary[j, :] .= [
            C,
            n_orbits,
            sum(rows[2, :]),
            n_looped_net,
            n_looped_net / n_orbits,
            n_looped_total,
            n_looped_total / n_orbits,
            sum(rows[5, :]) / n_orbits,
            maximum(rows[5, :]),
            sum(rows[6, :]) / n_orbits,
            maximum(rows[6, :]),
        ]

        @printf("C = %.2f: net %.0f/%d, total %.0f/%d\n",
                C, n_looped_net, n_orbits, n_looped_total, n_orbits)
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

println("Planar loop analysis")
planar_metrics, planar_summary, planar_C_values, mu = analyze_case(
    joinpath(DATA_DIR, "large_jacobi_scan_planar.jld2");
    spatial = false,
)

println("Spatial loop analysis")
spatial_metrics, spatial_summary, spatial_C_values, _ = analyze_case(
    joinpath(DATA_DIR, "large_jacobi_scan_spatial.jld2");
    spatial = true,
)

out_path = joinpath(DATA_DIR, "loop_before_r02.jld2")
save(out_path,
     "radius", R_LOOP,
     "loop_threshold", LOOP_THRESHOLD,
     "metric_names", METRIC_NAMES,
     "summary_names", SUMMARY_NAMES,
     "planar_metrics", planar_metrics,
     "planar_summary", planar_summary,
     "planar_C_values", planar_C_values,
     "spatial_metrics", spatial_metrics,
     "spatial_summary", spatial_summary,
     "spatial_C_values", spatial_C_values,
     "mu", mu)

write_summary_tsv(joinpath(DATA_DIR, "loop_before_r02_planar.tsv"), planar_summary)
write_summary_tsv(joinpath(DATA_DIR, "loop_before_r02_spatial.tsv"), spatial_summary)

println("Saved to $(relpath(out_path, joinpath(@__DIR__, "../../..")))")
