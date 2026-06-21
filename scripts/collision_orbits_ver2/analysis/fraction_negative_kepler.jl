using DataFrames
using Arrow
using CairoMakie

base_dir    = joinpath(@__DIR__, "../../../data/CL_ver2/planar")
results_dir = joinpath(@__DIR__, "../../../results/CL_ver2/planar")
mkpath(results_dir)

df       = DataFrame(Arrow.Table(joinpath(base_dir, "metadata_analyzed.arrow")))
C_values = sort(unique(df.jacobi_constant))
N_C      = length(C_values)

altitudes_km = [0, 50, 100, 200]

# --- 割合データの計算と保存 ---
frac_data = DataFrame(jacobi_constant = C_values)

for alt in altitudes_km
    col   = Symbol("E_kep_$(alt)km")
    fracs = [count(<(0), skipmissing(df[df.jacobi_constant .== C, col])) /
             sum(df.jacobi_constant .== C)  for C in C_values]
    frac_data[!, Symbol("frac_neg_$(alt)km")] = fracs
end

Arrow.write(joinpath(base_dir, "negative_kepler_fractions.arrow"), frac_data)
println("割合データ保存完了")

# phi グリッドはどの高度でも共通なのでループ外で1回だけ構築
df_s     = sort(df, [:jacobi_constant, :phi_0])
N_phi    = nrow(df_s) ÷ N_C
phi_vals = df_s[1:N_phi, :phi_0]

# --- 高度ごとにプロット ---
for alt in altitudes_km
    frac_col = Symbol("frac_neg_$(alt)km")
    ekep_col = Symbol("E_kep_$(alt)km")

    # Plot 1: ヤコビ定数 vs 割合
    fig1 = Figure(size = (700, 450))
    ax1  = Axis(fig1[1, 1],
        xlabel = "Jacobi constant C",
        ylabel = "Fraction of orbits with Eₖₑₚ < 0",
        title  = "Negative Kepler energy fraction  (altitude = $(alt) km)",
    )
    lines!(ax1, frac_data.jacobi_constant, frac_data[!, frac_col],
           color = :royalblue, linewidth = 2)
    save(joinpath(results_dir, "plot1_frac_neg_$(alt)km.png"), fig1, px_per_unit = 2)

    is_neg_vec = [coalesce(x < 0, false) for x in df_s[!, ekep_col]]
    is_neg_mat = Float64.(permutedims(reshape(is_neg_vec, N_phi, N_C)))

    fig2 = Figure(size = (950, 520))
    ax2  = Axis(fig2[1, 1],
        xlabel = "Jacobi constant C",
        ylabel = "Initial LC angle φ  [rad]",
        title  = "Orbits with Eₖₑₚ < 0  (altitude = $(alt) km)",
        yticks = (collect(0:π/4:π), ["0", "π/4", "π/2", "3π/4", "π"]),
    )
    hm = heatmap!(ax2, C_values, phi_vals, is_neg_mat,
                  colormap = [:lightgray, :steelblue], colorrange = (0.0, 1.0))
    Colorbar(fig2[1, 2], hm,
             ticks = ([0.0, 1.0], ["No", "Yes"]), width = 18)
    save(joinpath(results_dir, "plot2_phi_range_$(alt)km.png"), fig2, px_per_unit = 2)

    # Plot 3: φ vs 割合（全ヤコビ定数にわたる平均）
    # is_neg_mat の形状は (N_C, N_phi): 列方向に平均 → 各 φ での割合
    frac_by_phi = vec(sum(is_neg_mat, dims=1)) ./ N_C

    fig3 = Figure(size = (700, 450))
    ax3  = Axis(fig3[1, 1],
        xlabel = "Initial LC angle φ  [rad]",
        ylabel = "Fraction of C values with Eₖₑₚ < 0",
        title  = "Negative Kepler energy fraction by initial angle  (altitude = $(alt) km)",
        xticks          = (collect(0:π/4:π), ["0", "π/4", "π/2", "3π/4", "π"]),
        xminorticks     = IntervalsBetween(4),
        xminorgridvisible = true,
        yminorticks     = IntervalsBetween(4),
        yminorgridvisible = true,
    )
    lines!(ax3, phi_vals, frac_by_phi, color = :royalblue, linewidth = 2)
    save(joinpath(results_dir, "plot3_frac_by_phi_$(alt)km.png"), fig3, px_per_unit = 2)

    # Plot 4: Plot 3 の深掘り — 分母を「負エネルギー軌道が存在する C」に限定
    C_nonzero_mask  = vec(sum(is_neg_mat, dims=2)) .> 0   # 少なくとも1軌道が負エネルギーの C
    N_C_nonzero     = sum(C_nonzero_mask)
    frac_by_phi_nz  = vec(sum(is_neg_mat[C_nonzero_mask, :], dims=1)) ./ N_C_nonzero

    fig4 = Figure(size = (700, 450))
    ax4  = Axis(fig4[1, 1],
        xlabel = "Initial LC angle φ  [rad]",
        ylabel = "Fraction  (among C with Eₖₑₚ < 0 ∃)",
        title  = "Negative Kepler energy fraction (restricted C)  (altitude = $(alt) km)",
        xticks            = (collect(0:π/4:π), ["0", "π/4", "π/2", "3π/4", "π"]),
        xminorticks       = IntervalsBetween(4),
        xminorgridvisible = true,
        yminorticks       = IntervalsBetween(4),
        yminorgridvisible = true,
    )
    lines!(ax4, phi_vals, frac_by_phi_nz, color = :crimson, linewidth = 2)
    save(joinpath(results_dir, "plot4_frac_by_phi_restricted_$(alt)km.png"), fig4, px_per_unit = 2)

    println("$(alt) km 完了  (N_C_nonzero = $N_C_nonzero / $N_C)")
end

println("すべて完了")
