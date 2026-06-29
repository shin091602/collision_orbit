# 0. モジュールのインポート
using DifferentialEquations
using LinearAlgebra
using StaticArrays
using Arrow
using DataFrames
using Printf

# 1. ODEの問題設定
num_trj = 25000
mu = 0.01215058426994
t_fin = 2.0pi * 10
x0_temp = zeros(4)  # プレースホルダー: prob_func の remake で常に上書きされるため値は無意味
tspan = (0.0, t_fin)

# 平面円制限三体問題の運動方程式の定義
function eom_pcrtbp!(dx, x, mu, _)
	r1 = SA[x[1] + mu,       x[2]]
	r2 = SA[x[1] - 1.0 + mu, x[2]]
	nr1 = norm(r1); nr2 = norm(r2)
	dx[1] = x[3]
	dx[2] = x[4]
	dx[3] = -(1.0 - mu) * r1[1] / nr1^3 - mu * r2[1] / nr2^3 + 2.0x[4] + x[1]
	dx[4] = -(1.0 - mu) * r1[2] / nr1^3 - mu * r2[2] / nr2^3 - 2.0x[3] + x[2]
end

# 2. イベント処理の定義
# イベント処理#1: 半径方向速度から近地点を検知
function condition(x, _, _)
	return (x[1] + mu) * x[3] + x[2] * x[4]
end
affect!(_) = nothing

# イベント処理#2: 地球 or 月に近づきすぎると計算を終了する
function condition2(x, _, _)
	r1_sq = (x[1] + mu)^2       + x[2]^2
	r2_sq = (x[1] - 1.0 + mu)^2 + x[2]^2
	return r1_sq > 1.0e-6 && r2_sq > 1.0e-6
end
affect2!(integrator) = terminate!(integrator)

cb  = ContinuousCallback(condition, affect!, nothing, save_positions = (true, false))
cb2 = ContinuousCallback(condition2, affect2!, save_positions = (false, false))
cbs = CallbackSet(cb, cb2)

# 3. ヤコビ定数
function jacobi_C(u)
	r1 = sqrt((u[1] + mu)^2 + u[2]^2)
	r2 = sqrt((u[1] - 1.0 + mu)^2 + u[2]^2)
	Ω  = (u[1]^2 + u[2]^2) / 2.0 + (1.0 - mu) / r1 + mu / r2
	return 2.0Ω - (u[3]^2 + u[4]^2)
end

# 4. アンサンブル・シミュレーションの定義と実行
# C_poincare はループ内で上書き（solve 中は全スレッドで読み取り専用）
C_poincare = 3.15

function prob_func(prob, _)
	x0 = 3.0rand()
	φ0 = 2.0pi * rand()
	v0 = sqrt(max(0.0, x0^2 + 2.0 * (1.0 - mu) / abs(x0 + mu) + 2.0 * mu / abs(x0 - 1.0 + mu) - C_poincare))
	u0 = [x0, 0.0, v0 * cos(φ0), v0 * sin(φ0)]
	remake(prob, u0 = u0)
end

function output_func(sol, _)
	theta1_vec   = Float64[]
	r1_vec       = Float64[]
	a_vec        = Float64[]
	C_vec        = Float64[]
	prograde_vec = Bool[]
	X_vec        = Float64[]
	Y_vec        = Float64[]
	for u in sol.u[2:end-1]
		r1    = sqrt((u[1] + mu)^2 + u[2]^2)
		# 地球中心慣性系での速度（回転系から変換: v_inertial = v_rot + ω×r）
		vx_i  = u[3] - u[2]
		vy_i  = u[4] + u[1] + mu
		v1_sq = vx_i^2 + vy_i^2
		# 地球周回の軌道長半径（GM_Earth = 1 - mu）
		a     = 1.0 / (2.0 / r1 - v1_sq / (1.0 - mu))
		push!(theta1_vec,   atan(u[2], u[1] + mu))
		push!(r1_vec,       r1)
		push!(a_vec,        a)
		push!(C_vec,        jacobi_C(u))
		push!(prograde_vec, (u[1] + mu) * u[4] - u[2] * u[3] > 0)
		push!(X_vec,        u[1])
		push!(Y_vec,        u[2])
	end
	return (theta1 = theta1_vec, r1 = r1_vec, a = a_vec, C = C_vec,
	        prograde = prograde_vec, X = X_vec, Y = Y_vec), false
end

# 5. Cごとに計算してArrowファイルに保存
C_values = [2.90, 2.95, 3.00, 3.05, 3.10, 3.15]
out_dir  = joinpath(@__DIR__, "../../../data/CL_ver2/planar/poincare")
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

	# 全軌道の結果を結合
	all_theta1   = reduce(vcat, s.theta1   for s in sim.u; init = Float64[])
	all_r1       = reduce(vcat, s.r1       for s in sim.u; init = Float64[])
	all_a        = reduce(vcat, s.a        for s in sim.u; init = Float64[])
	all_C        = reduce(vcat, s.C        for s in sim.u; init = Float64[])
	all_prograde = reduce(vcat, s.prograde for s in sim.u; init = Bool[])
	all_X        = reduce(vcat, s.X        for s in sim.u; init = Float64[])
	all_Y        = reduce(vcat, s.Y        for s in sim.u; init = Float64[])

	df    = DataFrame(theta1 = all_theta1, r1 = all_r1, a = all_a, C = all_C,
	                  prograde = all_prograde, X = all_X, Y = all_Y)
	fname = @sprintf("crossings_C=%.2f_random.arrow", C_val)
	Arrow.write(joinpath(out_dir, fname), df)
	println("  saved $fname: $(nrow(df)) crossings")
end

println("Done.")
