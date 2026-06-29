## Poincare map plot — planar CR3BP (random orbits)
## Input:  data/CL_ver2/planar/poincare/crossings_C=X.XX_random.arrow
## Output: results/CL_ver2/planar/poincare/

using CairoMakie
using Arrow
using DataFrames
using Printf

μ = 0.01215058426994

data_dir    = joinpath(@__DIR__, "../../../data/CL_ver2/planar/poincare")
results_dir = joinpath(@__DIR__, "../../../results/CL_ver2/planar/poincare")
mkpath(results_dir)

# 近地点通過時の軌道長半径を計算（近地点では動径速度=0なので速度は接線方向のみ）
# 慣性系速度: prograde → v_i = v_rot + r1, retrograde → v_i = |v_rot - r1|
function semimajor_axis(r1, X, Y, C, prograde)
    r2    = sqrt((X - 1.0 + μ)^2 + Y^2)
    Ω     = (X^2 + Y^2) / 2.0 + (1.0 - μ) / r1 + μ / r2
    v_rot = sqrt(max(0.0, 2.0Ω - C))
    v_i   = prograde ? v_rot + r1 : abs(v_rot - r1)
    return 1.0 / (2.0 / r1 - v_i^2 / (1.0 - μ))
end

C_values = [2.90, 2.95, 3.00, 3.05, 3.10, 3.15]

for C_val in C_values
    fname = @sprintf("crossings_C=%.2f_random.arrow", C_val)
    fpath = joinpath(data_dir, fname)
    if !isfile(fpath)
        println("Not found: $fname  (skip)")
        continue
    end

    df = DataFrame(Arrow.Table(fpath))
    println("C=$(@sprintf("%.2f", C_val)): $(nrow(df)) crossings")

    # a 列が保存済みならそのまま使用、なければ導出
    a_vals = if "a" in names(df)
        df.a
    else
        semimajor_axis.(df.r1, df.X, df.Y, df.C, df.prograde)
    end

    fig = Figure(; size=(1024, 1024))
    ax  = Axis(fig[1, 1];
               xlabel = "theta1 [rad]",
               ylabel = "Semi-major axis a",
               title  = @sprintf("Poincare map  (C = %.2f, random orbits)", C_val))

    scatter!(ax, mod2pi.(df.theta1), a_vals; color=:black, markersize=1.5)

    xlims!(ax, 0.0, 2.0pi)
    ylims!(ax, 0.2, 0.8)

    outpath = joinpath(results_dir, @sprintf("poincare_C=%.2f_random.png", C_val))
    save(outpath, fig)
    println("  Saved: $(basename(outpath))")
end

println("Done.")
