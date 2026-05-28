using CairoMakie
using JLD2
using Printf

const DATA_DIR    = joinpath(@__DIR__, "../../../data")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan/spatial")
mkpath(RESULTS_DIR)

data     = load(joinpath(DATA_DIR, "jacobi_scan_spatial.jld2"))
all_cart = data["cart"]
all_ks   = data["ks"]
C_values = data["C_values"]

fmask(traj, rows) = vec(all(isfinite.(traj[rows, :]), dims = 1))
trajectory_colors(n) = cgrad(:turbo, n, categorical = true)

const planes = [
    (1, 2, "x", "y", "xy"),
    (2, 3, "y", "z", "yz"),
    (3, 1, "z", "x", "zx"),
]

function plot_cart_plane!(pos, trajs, colors, C, ri, rj, xl, yl, fname; zoom = false)
    ttl = "C = $(round(C, digits=2)), $(fname)"
    zoom && (ttl *= ", zoom")
    ax = Axis(pos, xlabel = xl, ylabel = yl, title = ttl)
    if zoom
        CairoMakie.xlims!(ax, ri == 1 ? (0.8, 1.2) : (-0.2, 0.2))
        CairoMakie.ylims!(ax, rj == 1 ? (0.8, 1.2) : (-0.2, 0.2))
    end
    for (i, traj) in enumerate(trajs)
        m = fmask(traj, [ri, rj])
        lines!(ax, traj[ri, m], traj[rj, m]; color = colors[i], linewidth = 0.5)
    end
    return ax
end

function plot_cart_xyz!(pos, trajs, colors, C; zoom = false)
    ttl = "C = $(round(C, digits=2)), xyz"
    zoom && (ttl *= ", zoom")
    ax = Axis3(pos, xlabel = "x", ylabel = "y", zlabel = "z", title = ttl)
    if zoom
        CairoMakie.xlims!(ax, 0.8, 1.2)
        CairoMakie.ylims!(ax, -0.2, 0.2)
        CairoMakie.zlims!(ax, -0.2, 0.2)
    end
    for (i, traj) in enumerate(trajs)
        m = fmask(traj, [1, 2, 3])
        lines!(ax, traj[1, m], traj[2, m], traj[3, m]; color = colors[i], linewidth = 0.5)
    end
    return ax
end

function plot_ks_component!(pos, trajs, colors, C, row, name)
    ax = Axis(pos, xlabel = "step", ylabel = name, title = "C = $(round(C, digits=2)), $(name)")
    for (i, traj) in enumerate(trajs)
        lines!(ax, axes(traj, 2), traj[row, :]; color = colors[i], linewidth = 0.5)
    end
    return ax
end

function plot_C(C)
    c_str = @sprintf("C%.2f", C)
    c_dir = joinpath(RESULTS_DIR, c_str)
    colors = trajectory_colors(length(all_ks[C]))
    mkpath(c_dir)

    # ── 2D projection planes ─────────────────────────────────────
    for (ri, rj, xl, yl, fname) in planes
        fig = Figure(size = (900, 700))
        plot_cart_plane!(fig[1, 1], all_cart[C], colors, C, ri, rj, xl, yl, fname)
        save(joinpath(c_dir, "cart_$(fname).pdf"), fig)

        fig_zoom = Figure(size = (900, 700))
        plot_cart_plane!(fig_zoom[1, 1], all_cart[C], colors, C, ri, rj, xl, yl, fname; zoom = true)
        save(joinpath(c_dir, "cart_$(fname)_zoom.pdf"), fig_zoom)
    end

    # ── 3D Cartesian view ────────────────────────────────────────
    fig3 = Figure(size = (900, 800))
    plot_cart_xyz!(fig3[1, 1], all_cart[C], colors, C)
    save(joinpath(c_dir, "cart_xyz.pdf"), fig3)

    fig3_zoom = Figure(size = (900, 800))
    plot_cart_xyz!(fig3_zoom[1, 1], all_cart[C], colors, C; zoom = true)
    save(joinpath(c_dir, "cart_xyz_zoom.pdf"), fig3_zoom)

    fig_cart_all = Figure(size = (1400, 1200))
    plot_cart_xyz!(fig_cart_all[1, 1], all_cart[C], colors, C)
    for (n, (ri, rj, xl, yl, fname)) in enumerate(planes)
        plot_cart_plane!(fig_cart_all[cld(n + 1, 2), iseven(n) ? 1 : 2],
                         all_cart[C], colors, C, ri, rj, xl, yl, fname)
    end
    save(joinpath(c_dir, "cart_all.pdf"), fig_cart_all)

    fig_cart_all_zoom = Figure(size = (1400, 1200))
    plot_cart_xyz!(fig_cart_all_zoom[1, 1], all_cart[C], colors, C; zoom = true)
    for (n, (ri, rj, xl, yl, fname)) in enumerate(planes)
        plot_cart_plane!(fig_cart_all_zoom[cld(n + 1, 2), iseven(n) ? 1 : 2],
                         all_cart[C], colors, C, ri, rj, xl, yl, fname; zoom = true)
    end
    save(joinpath(c_dir, "cart_all_zoom.pdf"), fig_cart_all_zoom)

    # ── KS components vs step index ──────────────────────────────
    ks_names = ["u1", "u2", "u3", "u4", "w1", "w2", "w3", "w4"]
    for (row, name) in enumerate(ks_names)
        fig = Figure(size = (900, 700))
        plot_ks_component!(fig[1, 1], all_ks[C], colors, C, row, name)
        save(joinpath(c_dir, "ks_$(name).pdf"), fig)
    end

    fig_ks_all = Figure(size = (1400, 1800))
    for (row, name) in enumerate(ks_names)
        plot_ks_component!(fig_ks_all[cld(row, 2), isodd(row) ? 1 : 2],
                           all_ks[C], colors, C, row, name)
    end
    save(joinpath(c_dir, "ks_all.pdf"), fig_ks_all)

    println("C = $C done: $(relpath(c_dir, joinpath(@__DIR__, "../../..")))")
end

for C in C_values
    plot_C(C)
end

println("Saved to $(relpath(RESULTS_DIR, joinpath(@__DIR__, "../../..")))")
