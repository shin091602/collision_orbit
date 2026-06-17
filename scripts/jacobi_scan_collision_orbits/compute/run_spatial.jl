using CairoMakie
using DifferentialEquations
using JLD2
using Printf

include("../../../src/ks_canonical_cr3bp.jl")
include("../../../src/utils_regularization.jl")

const mu = 0.01215058560962404

# Collision constraint: u = 0, |w|^2 = 8μ (KS analogue of the LC collision sphere)
const W_COLL = sqrt(8 * mu)

C_values = [2.94; 3.16]
N_ETA = 2
N_PHI = 3
N_PSI = 3
tspan  = 10pi   # fictitious time for backward propagation

cb_escape = ContinuousCallback((u, _, _) -> u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2 - 5, terminate!)

function collision_directions(n_eta, n_phi, n_psi)
    dirs = Vector{Vector{Float64}}()
    sizehint!(dirs, n_eta * n_phi * n_psi)

    for i in 1:n_eta
        # For S^3 in Hopf coordinates, ρ₁² is uniform on [0, 1].
        rho1 = sqrt((i - 0.5) / n_eta)
        rho2 = sqrt(1 - rho1^2)
        for j in 0:n_phi - 1
            phi = 2π * j / n_phi
            for k in 0:n_psi - 1
                psi = 2π * k / n_psi
                push!(dirs, [
                    rho1 * cos(phi),
                    rho1 * sin(phi),
                    rho2 * cos(psi),
                    rho2 * sin(psi),
                ])
            end
        end
    end

    return dirs
end

directions = collision_directions(N_ETA, N_PHI, N_PSI)

buf_cart = Vector{Vector{Matrix{Float64}}}(undef, length(C_values))
buf_ks   = Vector{Vector{Matrix{Float64}}}(undef, length(C_values))

Threads.@threads for j in eachindex(C_values)
    C = C_values[j]
    trajs_cart = Vector{Matrix{Float64}}()
    trajs_ks   = Vector{Matrix{Float64}}()

    for dir in directions
        zeta0 = vcat(zeros(4), W_COLL .* dir)
        prob = ODEProblem(ks_canonical_cr3bp, zeta0, (0.0, -tspan), [mu, C])
        sol  = solve(prob, Vern9(); abstol = 1e-13, reltol = 1e-13, maxiters = Int(1e6), callback = cb_escape)
        num_points = length(sol.t)
        traj_cart = zeros(6, num_points)
        traj_ks = zeros(8, num_points)
        for i in 1:num_points
            traj_cart[:, i] = ks2cart(sol.u[i], mu)
            traj_ks[:, i] = sol.u[i]
        end
        push!(trajs_cart, traj_cart)
        push!(trajs_ks,   traj_ks)
    end

    buf_cart[j] = trajs_cart
    buf_ks[j]   = trajs_ks
    println("C = $C: $(length(trajs_cart)) orbits")
end

all_cart = Dict(C_values[j] => buf_cart[j] for j in eachindex(C_values))
all_ks   = Dict(C_values[j] => buf_ks[j]   for j in eachindex(C_values))

const DATA_DIR = joinpath(@__DIR__, "../../../data")
const DATA_PATH = joinpath(DATA_DIR, "jacobi_scan_spatial.jld2")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan/spatial_3d")
const MOON_POS = (1 - mu, 0.0, 0.0)
const ZOOM_XLIMS = (0.8, 1.2)
const ZOOM_YLIMS = (-0.2, 0.2)
const ZOOM_ZLIMS = (-0.2, 0.2)
const ZOOM_CLIP_RADIUS = 0.2
const SINGULARITY_RADIUS = 0.012
const SINGULARITY_RADIUS_OVERALL = 0.008
mkpath(DATA_DIR)
save(DATA_PATH,
     "cart",      all_cart,
     "ks",        all_ks,
     "C_values",  C_values,
     "mu",        mu,
     "directions", directions,
     "N_ETA",     N_ETA,
     "N_PHI",     N_PHI,
     "N_PSI",     N_PSI)
println("Saved to $(relpath(DATA_PATH, joinpath(@__DIR__, "../../..")))")

fmask(traj, rows) = vec(all(isfinite.(traj[rows, :]), dims = 1))
trajectory_colors(n) = cgrad(:turbo, n, categorical = true)

function zoom_prefix_mask(traj)
    finite = fmask(traj, [1, 2, 3])
    r_moon = sqrt.((traj[1, :] .- MOON_POS[1]).^2 .+
                   (traj[2, :] .- MOON_POS[2]).^2 .+
                   (traj[3, :] .- MOON_POS[3]).^2)
    inside = finite .&
             (ZOOM_XLIMS[1] .<= traj[1, :] .<= ZOOM_XLIMS[2]) .&
             (ZOOM_YLIMS[1] .<= traj[2, :] .<= ZOOM_YLIMS[2]) .&
             (ZOOM_ZLIMS[1] .<= traj[3, :] .<= ZOOM_ZLIMS[2]) .&
             (r_moon .<= ZOOM_CLIP_RADIUS)
    outside_idx = findfirst(!, inside)
    if isnothing(outside_idx)
        return inside
    end
    m = falses(length(inside))
    outside_idx > 1 && (m[1:outside_idx - 1] .= true)
    return m
end

function position_widths(trajs)
    mins = fill(Inf, 3)
    maxs = fill(-Inf, 3)
    for traj in trajs
        m = fmask(traj, [1, 2, 3])
        any(m) || continue
        for row in 1:3
            vals = traj[row, m]
            mins[row] = min(mins[row], minimum(vals))
            maxs[row] = max(maxs[row], maximum(vals))
        end
    end
    widths = maxs .- mins
    return max.(widths, eps(Float64))
end

function overall_singularity_markersize(trajs)
    widths = position_widths(trajs)
    return Vec3f(SINGULARITY_RADIUS_OVERALL .* widths ./ maximum(widths))
end

function plot_cart_xyz_poster!(pos, trajs, colors, C; zoom = false)
    ttl = "Spatial collision orbits\nC = $(@sprintf("%.2f", C))"
    zoom && (ttl *= " (zoom)")
    ax = Axis3(pos;
               xlabel = "x",
               ylabel = "y",
               zlabel = "z",
               title = ttl,
               titlesize = 87,
               titlegap = 24,
               xlabelsize = 36,
               ylabelsize = 36,
               zlabelsize = 36,
               xticklabelsize = 24,
               yticklabelsize = 24,
               zticklabelsize = 24,
               aspect = zoom ? :data : (1, 1, 1))
    if zoom
        CairoMakie.xlims!(ax, ZOOM_XLIMS...)
        CairoMakie.ylims!(ax, ZOOM_YLIMS...)
        CairoMakie.zlims!(ax, ZOOM_ZLIMS...)
    end
    meshscatter!(ax,
                 [MOON_POS[1]], [MOON_POS[2]], [MOON_POS[3]];
                 markersize = zoom ? SINGULARITY_RADIUS : overall_singularity_markersize(trajs),
                 color = :black,
                 shininess = 32)
    for (i, traj) in enumerate(trajs)
        m = zoom ? zoom_prefix_mask(traj) : fmask(traj, [1, 2, 3])
        any(m) || continue
        lines!(ax, traj[1, m], traj[2, m], traj[3, m]; color = colors[i], linewidth = 1.2)
    end
    return ax
end

function save_spatial_3d_plots(C)
    c_str = @sprintf("C%.2f", C)
    c_dir = joinpath(RESULTS_DIR, c_str)
    colors = trajectory_colors(length(all_cart[C]))
    mkpath(c_dir)

    fig = Figure(size = (1400, 1450), fontsize = 28, figure_padding = (30, 30, 30, 190))
    plot_cart_xyz_poster!(fig[1, 1], all_cart[C], colors, C)
    save(joinpath(c_dir, "cart_xyz_3d.pdf"), fig)

    fig_zoom = Figure(size = (1400, 1450), fontsize = 28, figure_padding = (30, 30, 30, 190))
    plot_cart_xyz_poster!(fig_zoom[1, 1], all_cart[C], colors, C; zoom = true)
    save(joinpath(c_dir, "cart_xyz_3d_zoom.pdf"), fig_zoom)

    println("C = $C 3D plots: $(relpath(c_dir, joinpath(@__DIR__, "../../..")))")
end

mkpath(RESULTS_DIR)
for C in C_values
    save_spatial_3d_plots(C)
end

println("Saved 3D plots to $(relpath(RESULTS_DIR, joinpath(@__DIR__, "../../..")))")
