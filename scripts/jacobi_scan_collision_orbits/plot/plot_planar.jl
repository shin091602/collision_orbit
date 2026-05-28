using CairoMakie
using JLD2

const DATA_DIR    = joinpath(@__DIR__, "../../../data")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan/planar")
mkpath(RESULTS_DIR)

data     = load(joinpath(DATA_DIR, "jacobi_scan_planar.jld2"))
all_cart = data["cart"]
all_lc   = data["lc"]
C_values = data["C_values"]

fmask(traj, rows) = vec(all(isfinite.(traj[rows, :]), dims = 1))
trajectory_colors(n) = cgrad(:turbo, n, categorical = true)

const planes = [(1, 2, "x", "y", "xy")]

function plot_C(C)
    c_str = "C$(lpad(round(C, digits=2), 4, '0'))"
    ttl   = "C = $(round(C, digits=2))"
    c_dir = joinpath(RESULTS_DIR, c_str)
    colors = trajectory_colors(length(all_lc[C]))
    mkpath(c_dir)

    # ── 2D projection planes ─────────────────────────────────────
    for (ri, rj, xl, yl, fname) in planes
        fig = Figure(size = (900, 700))
        ax = Axis(fig[1, 1], xlabel = xl, ylabel = yl, title = ttl)
        for (i, traj) in enumerate(all_cart[C])
            m = fmask(traj, [ri, rj])
            lines!(ax, traj[ri, m], traj[rj, m]; color = colors[i], linewidth = 0.5)
        end
        save(joinpath(c_dir, "cart_$(fname).pdf"), fig)

        fig_zoom = Figure(size = (900, 700))
        ax_zoom = Axis(fig_zoom[1, 1],
                       xlabel = xl,
                       ylabel = yl,
                       title = "$(ttl), zoom")
        CairoMakie.xlims!(ax_zoom, 0.8, 1.2)
        CairoMakie.ylims!(ax_zoom, -0.2, 0.2)
        for (i, traj) in enumerate(all_cart[C])
            m = fmask(traj, [ri, rj])
            lines!(ax_zoom, traj[ri, m], traj[rj, m]; color = colors[i], linewidth = 0.5)
        end
        save(joinpath(c_dir, "cart_$(fname)_zoom.pdf"), fig_zoom)
    end

    # ── LC components vs step index ──────────────────────────────
    for (row, name) in enumerate(["u1", "u2", "w1", "w2"])
        fig = Figure(size = (900, 700))
        ax = Axis(fig[1, 1], xlabel = "step", ylabel = name, title = ttl)
        for (i, traj) in enumerate(all_lc[C])
            lines!(ax, axes(traj, 2), traj[row, :]; color = colors[i], linewidth = 0.5)
        end
        save(joinpath(c_dir, "lc_$(name).pdf"), fig)
    end

    println("C = $C done: $(relpath(c_dir, joinpath(@__DIR__, "../../..")))")
end

for C in C_values
    plot_C(C)
end

println("Saved to $(relpath(RESULTS_DIR, joinpath(@__DIR__, "../../..")))")
