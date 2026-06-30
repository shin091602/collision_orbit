## 地球近点ポアンカレマップ + 月衝突軌道オーバーレイ
## Input:  data/CL_ver2/planar/poincare/crossings_C=X.XX_random.arrow      (背景: ランダム軌道)
##         data/CL_ver2/planar/poincare/crossings_lc_moon_C=X.XX.arrow     (前景: 月衝突軌道)
## Output: results/CL_ver2/planar/poincare/

using CairoMakie
using Arrow
using DataFrames
using Printf

data_dir    = joinpath(@__DIR__, "../../../data/CL_ver2/planar/poincare")
results_dir = joinpath(@__DIR__, "../../../results/CL_ver2/planar/poincare")
mkpath(results_dir)

C_values = [2.90, 2.95, 3.00, 3.05, 3.10, 3.15]

for C_val in C_values
    rand_path = joinpath(data_dir, @sprintf("crossings_C=%.2f_random.arrow",   C_val))
    moon_path = joinpath(data_dir, @sprintf("crossings_lc_moon_C=%.2f.arrow",  C_val))

    has_rand = isfile(rand_path)
    has_moon = isfile(moon_path)
    (!has_rand && !has_moon) && (println("No data for C=$(@sprintf("%.2f",C_val)), skip"); continue)

    fig = Figure(size = (1024, 1024))
    ax  = Axis(fig[1, 1];
               xlabel = "theta1  [rad]",
               ylabel = "Semi-major axis a",
               title  = @sprintf("Earth-perigee Poincare map  (C = %.2f)", C_val))

    if has_rand
        df = DataFrame(Arrow.Table(rand_path))
        scatter!(ax, mod2pi.(df.theta1), df.a;
                 color      = (:gray50, 0.35),
                 markersize = 1.5,
                 label      = "random orbits")
        println("C=$(@sprintf("%.2f",C_val))  random : $(nrow(df)) crossings")
    end

    if has_moon
        df = DataFrame(Arrow.Table(moon_path))
        scatter!(ax, mod2pi.(df.theta1), df.a;
                 color      = (:red, 0.8),
                 markersize = 2.5,
                 label      = "Moon collision (LC)")
        println("C=$(@sprintf("%.2f",C_val))  LC moon: $(nrow(df)) crossings")
    end

    xlims!(ax, 0.0, 2.0pi)
    ylims!(ax, 0.2, 0.8)

    Legend(fig[1, 2],
           ax;
           framevisible = false)

    outpath = joinpath(results_dir, @sprintf("poincare_overlay_C=%.2f.png", C_val))
    save(outpath, fig)
    println("  Saved: $(basename(outpath))")
end

println("Done.")
