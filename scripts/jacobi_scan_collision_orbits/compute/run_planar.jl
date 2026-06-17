using CairoMakie
using DifferentialEquations
using JLD2
using Printf

include("../../../src/lc_canonical_cr3bp.jl")
include("../../../src/utils_regularization.jl")

const mu = 0.01215058560962404

# Collision constraint: u = 0, |w|² = 8μ (from regularized energy at r_M = 0)
const W_COLL = sqrt(8 * mu)

C_values = [2.94; 3.16]
N_angles = 18
tspan  = 10pi   # fictitious time for backward propagation

cb_escape = ContinuousCallback((u, _, _) -> u[1]^2 + u[2]^2 - 5, terminate!)

buf_cart = Vector{Vector{Matrix{Float64}}}(undef, length(C_values))
buf_lc   = Vector{Vector{Matrix{Float64}}}(undef, length(C_values))

Threads.@threads for j in eachindex(C_values)
    C = C_values[j]
    trajs_cart = Vector{Matrix{Float64}}()
    trajs_lc   = Vector{Matrix{Float64}}()

    for k in 0:N_angles - 1
        alpha = 2π * k / N_angles
        zeta0 = [0.0, 0.0, W_COLL * cos(alpha), W_COLL * sin(alpha)]
        prob = ODEProblem(lc_canonical_cr3bp, zeta0, (0.0, -tspan), [mu, C])
        sol  = solve(prob, Vern9(); abstol = 1e-13, reltol = 1e-13, maxiters = Int(1e6), callback = cb_escape)
        num_points = length(sol.t)
        traj_cart = zeros(4, num_points)
        traj_lc = zeros(4, num_points)
        for i in 1:num_points
            traj_cart[:, i] = lc2cart(sol.u[i], mu)
            traj_lc[:, i] = sol.u[i]
        end
        push!(trajs_cart, traj_cart)
        push!(trajs_lc,   traj_lc)
    end

    buf_cart[j] = trajs_cart
    buf_lc[j]   = trajs_lc
    println("C = $C: $(length(trajs_cart)) orbits")
end

all_cart = Dict(C_values[j] => buf_cart[j] for j in eachindex(C_values))
all_lc   = Dict(C_values[j] => buf_lc[j]   for j in eachindex(C_values))

const DATA_DIR = joinpath(@__DIR__, "../../../data")
const DATA_PATH = joinpath(DATA_DIR, "jacobi_scan_planar.jld2")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan/planar_poster")
const MOON_POS = (1 - mu, 0.0)
const ZOOM_XLIMS = (0.8, 1.2)
const ZOOM_YLIMS = (-0.2, 0.2)
const SINGULARITY_MARKERSIZE = 22
mkpath(DATA_DIR)
save(DATA_PATH,
     "cart",     all_cart,
     "lc",       all_lc,
     "C_values", C_values,
     "mu",       mu)
println("Saved to $(relpath(DATA_PATH, joinpath(@__DIR__, "../../..")))")

fmask(traj, rows) = vec(all(isfinite.(traj[rows, :]), dims = 1))
trajectory_colors(n) = cgrad(:turbo, n, categorical = true)

function zoom_prefix_mask(traj)
    finite = fmask(traj, [1, 2])
    inside = finite .&
             (ZOOM_XLIMS[1] .<= traj[1, :] .<= ZOOM_XLIMS[2]) .&
             (ZOOM_YLIMS[1] .<= traj[2, :] .<= ZOOM_YLIMS[2])
    outside_idx = findfirst(!, inside)
    if isnothing(outside_idx)
        return inside
    end
    m = falses(length(inside))
    outside_idx > 1 && (m[1:outside_idx - 1] .= true)
    return m
end

function plot_cart_xy_poster!(pos, trajs, colors, C; zoom = false)
    ttl = "Planar collision orbits\nC = $(@sprintf("%.2f", C))"
    zoom && (ttl *= " (zoom)")
    ax = Axis(pos;
              xlabel = "x",
              ylabel = "y",
              title = ttl,
              titlesize = 87,
              titlelineheight = 0.9,
              xlabelsize = 36,
              ylabelsize = 36,
              xticklabelsize = 24,
              yticklabelsize = 24,
              aspect = DataAspect())
    if zoom
        CairoMakie.xlims!(ax, ZOOM_XLIMS...)
        CairoMakie.ylims!(ax, ZOOM_YLIMS...)
    end
    scatter!(ax,
             [MOON_POS[1]], [MOON_POS[2]];
             marker = :circle,
             markersize = SINGULARITY_MARKERSIZE,
             color = :black)
    for (i, traj) in enumerate(trajs)
        m = zoom ? zoom_prefix_mask(traj) : fmask(traj, [1, 2])
        any(m) || continue
        lines!(ax, traj[1, m], traj[2, m]; color = colors[i], linewidth = 1.2)
    end
    return ax
end

function save_planar_poster_plots(C)
    c_str = @sprintf("C%.2f", C)
    c_dir = joinpath(RESULTS_DIR, c_str)
    colors = trajectory_colors(length(all_cart[C]))
    mkpath(c_dir)

    fig = Figure(size = (1400, 1200), fontsize = 28)
    plot_cart_xy_poster!(fig[1, 1], all_cart[C], colors, C)
    save(joinpath(c_dir, "cart_xy_poster.pdf"), fig)

    fig_zoom = Figure(size = (1400, 1200), fontsize = 28)
    plot_cart_xy_poster!(fig_zoom[1, 1], all_cart[C], colors, C; zoom = true)
    save(joinpath(c_dir, "cart_xy_poster_zoom.pdf"), fig_zoom)

    println("C = $C poster plots: $(relpath(c_dir, joinpath(@__DIR__, "../../..")))")
end

mkpath(RESULTS_DIR)
for C in C_values
    save_planar_poster_plots(C)
end

println("Saved poster plots to $(relpath(RESULTS_DIR, joinpath(@__DIR__, "../../..")))")
