# 月近点ポアンカレマップ — 平面 CR3BP（ランダム軌道）
# poincare_planar.jl（地球近点版）の月近点バージョン
#
# 断面条件: 月相対動径速度 = 0  →  (X-(1-μ))*VX + Y*VY = 0
# 出力:     data/CL_ver2/planar/poincare_moon/crossings_moon_C=X.XX_random.arrow

using DifferentialEquations
using LinearAlgebra
using StaticArrays
using Arrow
using DataFrames
using Printf

# ───────────────────────────────────────────────
# 1. 定数・問題設定
# ───────────────────────────────────────────────
num_trj = 10000
mu  = 0.01215058426994
t_fin = 2.0pi * 10
x0_temp = zeros(4)
tspan = (0.0, t_fin)

# ───────────────────────────────────────────────
# 2. 運動方程式（地球近点版と同じ）
# ───────────────────────────────────────────────
function eom_pcrtbp!(dx, x, mu, _)
	r1 = SA[x[1] + mu,        x[2]]
	r2 = SA[x[1] - 1.0 + mu,  x[2]]
	nr1 = norm(r1); nr2 = norm(r2)
	dx[1] = x[3]
	dx[2] = x[4]
	dx[3] = -(1.0 - mu) * r1[1] / nr1^3 - mu * r2[1] / nr2^3 + 2.0x[4] + x[1]
	dx[4] = -(1.0 - mu) * r1[2] / nr1^3 - mu * r2[2] / nr2^3 - 2.0x[3] + x[2]
end

# ───────────────────────────────────────────────
# 3. イベント処理
# ───────────────────────────────────────────────
# 断面条件: 月相対動径速度 = 0（月近点を検知）
function condition_moon(x, _, _)
	return (x[1] - 1.0 + mu) * x[3] + x[2] * x[4]
end
affect!(_) = nothing

# 終端条件: 地球または月に近づきすぎた場合に計算を終了する
function condition_too_close(x, _, _)
	r1_sq = (x[1] + mu)^2        + x[2]^2
	r2_sq = (x[1] - 1.0 + mu)^2 + x[2]^2
	return r1_sq > 1.0e-6 && r2_sq > 1.0e-6
end
affect_terminate!(integrator) = terminate!(integrator)

cb  = ContinuousCallback(condition_moon, affect!, nothing, save_positions = (true, false))
cb2 = ContinuousCallback(condition_too_close, affect_terminate!, save_positions = (false, false))
cbs = CallbackSet(cb, cb2)

# ───────────────────────────────────────────────
# 4. ヤコビ定数
# ───────────────────────────────────────────────
function jacobi_C(u)
	r1 = sqrt((u[1] + mu)^2 + u[2]^2)
	r2 = sqrt((u[1] - 1.0 + mu)^2 + u[2]^2)
	Ω  = (u[1]^2 + u[2]^2) / 2.0 + (1.0 - mu) / r1 + mu / r2
	return 2.0Ω - (u[3]^2 + u[4]^2)
end

# ───────────────────────────────────────────────
# 5. アンサンブル設定
# ───────────────────────────────────────────────
C_poincare = 3.15  # ループ内で上書き（solve 中は読み取り専用）

function prob_func(prob, _)
	x0 = 3.0rand()
	φ0 = 2.0pi * rand()
	v0 = sqrt(max(0.0, x0^2 + 2.0 * (1.0 - mu) / abs(x0 + mu) + 2.0 * mu / abs(x0 - 1.0 + mu) - C_poincare))
	u0 = [x0, 0.0, v0 * cos(φ0), v0 * sin(φ0)]
	remake(prob, u0 = u0)
end

function output_func(sol, _)
	theta2_vec   = Float64[]
	r2_vec       = Float64[]
	a_vec        = Float64[]   # 地球相対半長径（常に正: GM_Earth = 1 - mu）
	C_vec        = Float64[]
	prograde_vec = Bool[]
	X_vec        = Float64[]
	Y_vec        = Float64[]
	for u in sol.u[2:end-1]
		r2    = sqrt((u[1] - 1.0 + mu)^2 + u[2]^2)
		r1    = sqrt((u[1] + mu)^2 + u[2]^2)
		# 地球中心慣性系速度（地球近点版と同じ変換）
		vx_i  = u[3] - u[2]
		vy_i  = u[4] + u[1] + mu
		v1_sq = vx_i^2 + vy_i^2
		# 地球周回の軌道長半径（GM_Earth = 1 - mu）
		a     = 1.0 / (2.0 / r1 - v1_sq / (1.0 - mu))
		push!(theta2_vec,   atan(u[2], u[1] - 1.0 + mu))
		push!(r2_vec,       r2)
		push!(a_vec,        a)
		push!(C_vec,        jacobi_C(u))
		push!(prograde_vec, (u[1] - 1.0 + mu) * u[4] - u[2] * u[3] > 0)
		push!(X_vec,        u[1])
		push!(Y_vec,        u[2])
	end
	return (theta2 = theta2_vec, r2 = r2_vec, a = a_vec, C = C_vec,
	        prograde = prograde_vec, X = X_vec, Y = Y_vec), false
end

# ───────────────────────────────────────────────
# 6. C ごとに計算して Arrow ファイルに保存
# ───────────────────────────────────────────────
C_values = [2.90, 2.95, 3.00, 3.05, 3.10, 3.15]
out_dir  = joinpath(@__DIR__, "../../../data/CL_ver2/planar/poincare_moon")
mkpath(out_dir)

prob = ODEProblem(eom_pcrtbp!, x0_temp, tspan, mu)

for C_val in C_values
	global C_poincare = C_val
	println("C = $C_val ...")

	ensemble_prob = EnsembleProblem(prob;
	                                prob_func   = prob_func,
	                                output_func = output_func)
	sim = solve(ensemble_prob, Vern7(), EnsembleThreads();
	            trajectories = num_trj,
	            callback     = cbs,
	            saveat       = t_fin)

	all_theta2   = reduce(vcat, s.theta2   for s in sim.u; init = Float64[])
	all_r2       = reduce(vcat, s.r2       for s in sim.u; init = Float64[])
	all_a        = reduce(vcat, s.a        for s in sim.u; init = Float64[])
	all_C        = reduce(vcat, s.C        for s in sim.u; init = Float64[])
	all_prograde = reduce(vcat, s.prograde for s in sim.u; init = Bool[])
	all_X        = reduce(vcat, s.X        for s in sim.u; init = Float64[])
	all_Y        = reduce(vcat, s.Y        for s in sim.u; init = Float64[])

	df    = DataFrame(theta2 = all_theta2, r2 = all_r2, a = all_a, C = all_C,
	                  prograde = all_prograde, X = all_X, Y = all_Y)
	fname = @sprintf("crossings_moon_C=%.2f_random.arrow", C_val)
	Arrow.write(joinpath(out_dir, fname), df)
	println("  saved $fname: $(nrow(df)) crossings")
end

println("Done.")
