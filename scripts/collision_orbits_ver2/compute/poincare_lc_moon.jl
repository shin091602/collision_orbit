# LC正則化を用いた月衝突軌道の地球近点ポアンカレ断面計算
#
# 月特異点 u=[0,0]（|w|=√(8μ)、φグリッド）から後方積分し、
# 地球近点通過（(X+μ)·VX + Y·VY = 0）を検出して保存。
# poincare_planar.jl と同一スキーマ → plot_poincare_overlay.jl でオーバーレイ可
#
# 出力: data/CL_ver2/planar/poincare/crossings_lc_moon_C=X.XX.arrow

using DifferentialEquations
using Arrow
using DataFrames
using Printf

include("../../../src/lc_canonical_cr3bp.jl")
include("../../../src/utils_regularization.jl")
include("utils.jl")

# poincare_planar.jl と μ を統一
const mu_lc = 0.01215058426994

# 月特異点での初期条件グリッド
# R = √(8μ): LC正則化ハミルトニアン K|_{u=0} = |w|²/8 - μ = 0 から導出
# φ ∈ [0, π): LC角度 → 月相対物理角度 2φ ∈ [0, 2π) を一様に覆う
const N_grid   = 2000
const tspan_lc = 20.0pi   # LC虚時間の後方積分幅

w_grid, phi_grid = generate_planar_ks_collision_grid(N_grid, mu_lc)

# 地球近点条件: (X+μ)·VX + Y·VY = 0
# lc2cart で [X, Y, px, py] を取得し、VX = px + Y、VY = py - X を計算
function cond_earth_peri(zeta, _, _)
    zeta[1]^2 + zeta[2]^2 < 1e-10 && return 1.0   # 月特異点付近を回避
    cart = lc2cart(zeta, mu_lc)
    X, Y, px, py = cart[1], cart[2], cart[3], cart[4]
    VX = px + Y
    VY = py - X
    return (X + mu_lc) * VX + Y * VY
end

C_values = [2.90, 2.95, 3.00, 3.05, 3.10, 3.15]
out_dir  = joinpath(@__DIR__, "../../../data/CL_ver2/planar/poincare")
mkpath(out_dir)

println("スレッド数: $(Threads.nthreads())")

for C_val in C_values
    println("C = $C_val ...")

    buf_theta1   = Vector{Vector{Float64}}(undef, N_grid)
    buf_a        = Vector{Vector{Float64}}(undef, N_grid)
    buf_prograde = Vector{Vector{Bool}}(undef,   N_grid)
    buf_X        = Vector{Vector{Float64}}(undef, N_grid)
    buf_Y        = Vector{Vector{Float64}}(undef, N_grid)

    Threads.@threads for i in 1:N_grid
        # コールバックはスレッドごとに生成（内部キャッシュを分離）
        cb = ContinuousCallback(cond_earth_peri, (integrator) -> nothing, nothing;
                                save_positions = (true, false))

        w0    = w_grid[i]
        zeta0 = [0.0, 0.0, Float64(w0[1]), Float64(w0[2])]
        prob  = ODEProblem(lc_canonical_cr3bp, zeta0, (0.0, -tspan_lc), [mu_lc, C_val])
        sol   = solve(prob, Vern9();
                      abstol         = 1e-12,
                      reltol         = 1e-12,
                      maxiters       = Int(1e6),
                      callback       = cb,
                      save_everystep = false,
                      saveat         = [-tspan_lc])

        theta1_i   = Float64[]
        a_i        = Float64[]
        prograde_i = Bool[]
        X_i        = Float64[]
        Y_i        = Float64[]

        # sol.u[1] = 月特異点（初期値）、sol.u[end] = 終端、sol.u[2:end-1] = 断面通過
        for u in sol.u[2:end-1]
            u[1]^2 + u[2]^2 < 1e-10 && continue
            cart = lc2cart(u, mu_lc)::Vector{Float64}
            X  = cart[1]; Y  = cart[2]
            px = cart[3]; py = cart[4]
            VX = px + Y
            VY = py - X
            r1 = sqrt((X + mu_lc)^2 + Y^2)
            # 地球相対慣性系速度: vx = VX - Y = px, vy = VY + X + μ = py + μ
            v1_sq = px^2 + (py + mu_lc)^2
            a = 1.0 / (2.0 / r1 - v1_sq / (1.0 - mu_lc))
            push!(theta1_i,   atan(Y, X + mu_lc))
            push!(a_i,        a)
            push!(prograde_i, (X + mu_lc) * VY - Y * VX > 0)
            push!(X_i,        X)
            push!(Y_i,        Y)
        end

        buf_theta1[i]   = theta1_i
        buf_a[i]        = a_i
        buf_prograde[i] = prograde_i
        buf_X[i]        = X_i
        buf_Y[i]        = Y_i
    end

    all_theta1   = reduce(vcat, buf_theta1;   init = Float64[])
    all_a        = reduce(vcat, buf_a;        init = Float64[])
    all_prograde = reduce(vcat, buf_prograde; init = Bool[])
    all_X        = reduce(vcat, buf_X;        init = Float64[])
    all_Y        = reduce(vcat, buf_Y;        init = Float64[])

    df    = DataFrame(theta1   = all_theta1,
                      a        = all_a,
                      C        = fill(C_val, length(all_theta1)),
                      prograde = all_prograde,
                      X        = all_X,
                      Y        = all_Y)
    fname = @sprintf("crossings_lc_moon_C=%.2f.arrow", C_val)
    Arrow.write(joinpath(out_dir, fname), df)
    println("  saved $fname: $(nrow(df)) crossings")
end
println("Done.")
