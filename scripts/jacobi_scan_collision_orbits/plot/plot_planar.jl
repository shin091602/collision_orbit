using JLD2, Plots

const DATA_DIR    = joinpath(@__DIR__, "../../../data")
const RESULTS_DIR = joinpath(@__DIR__, "../../../results/jacobi_scan")
mkpath(RESULTS_DIR)

data     = load(joinpath(DATA_DIR, "jacobi_scan_planar.jld2"))
all_cart = data["cart"]
all_lc   = data["lc"]
C_values = data["C_values"]

clrs = cgrad(:viridis, length(C_values), categorical = true)
fmask(traj, rows) = vec(all(isfinite.(traj[rows, :]), dims = 1))

const planes = [(1, 2, "x",  "y",  "xy"),
                (1, 3, "x",  "ẋ",  "x_xdot"),
                (2, 4, "y",  "ẏ",  "y_ydot"),
                (3, 4, "ẋ", "ẏ",  "xdot_ydot")]

function plot_C(j, C)
    c_str = "C$(lpad(round(C, digits=2), 4, '0'))"
    col   = clrs[j]
    ttl   = "C = $(round(C, digits=2))"

    # ── 3D phase portrait (x, y, ẋ) ─────────────────────────────
    fig = plot3d(; xlabel = "x", ylabel = "y", zlabel = "ẋ", title = ttl, legend = false)
    for traj in all_cart[C]
        m = fmask(traj, [1, 2, 3])
        plot!(fig, traj[1, m], traj[2, m], traj[3, m]; color = col, lw = 0.5)
    end
    savefig(fig, joinpath(RESULTS_DIR, "cart_3d_$(c_str).pdf"))

    # ── 2D projection planes ─────────────────────────────────────
    for (ri, rj, xl, yl, fname) in planes
        fig = plot(; xlabel = xl, ylabel = yl, title = ttl, legend = false)
        for traj in all_cart[C]
            m = fmask(traj, [ri, rj])
            plot!(fig, traj[ri, m], traj[rj, m]; color = col, lw = 0.5)
        end
        savefig(fig, joinpath(RESULTS_DIR, "cart_$(fname)_$(c_str).pdf"))
    end

    # ── LC components vs step index ──────────────────────────────
    for (row, name) in enumerate(["u1", "u2", "w1", "w2"])
        fig = plot(; xlabel = "step", ylabel = name, title = ttl, legend = false)
        for traj in all_lc[C]
            plot!(fig, traj[row, :]; color = col, lw = 0.5)
        end
        savefig(fig, joinpath(RESULTS_DIR, "lc_$(name)_$(c_str).pdf"))
    end

    println("C = $C done")
end

for (j, C) in enumerate(C_values)
    plot_C(j, C)
end

println("Saved → results/jacobi_scan/")
