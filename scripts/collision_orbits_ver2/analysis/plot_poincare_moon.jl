## 月近点ポアンカレマップ描画 — 平面 CR3BP（ランダム軌道）
## Input:  data/CL_ver2/planar/poincare_moon/crossings_moon_C=X.XX_random.arrow
## Output: results/CL_ver2/planar/poincare_moon/
##
## y 軸: 地球相対軌道長半径 a（ケプラー6要素、GM_Earth = 1 - μ）
##   月フライバイ時の地球周回軌道エネルギーを表す（常に正）

using CairoMakie
using Arrow
using DataFrames
using Printf

mu = 0.01215058426994

data_dir    = joinpath(@__DIR__, "../../../data/CL_ver2/planar/poincare_moon")
results_dir = joinpath(@__DIR__, "../../../results/CL_ver2/planar/poincare_moon")
mkpath(results_dir)

C_values = [2.90, 2.95, 3.00, 3.05, 3.10, 3.15]

for C_val in C_values
	fname = @sprintf("crossings_moon_C=%.2f_random.arrow", C_val)
	fpath = joinpath(data_dir, fname)
	if !isfile(fpath)
		println("Not found: $fname  (skip)")
		continue
	end

	df = DataFrame(Arrow.Table(fpath))
	println("C=$(@sprintf("%.2f", C_val)): $(nrow(df)) crossings")

	fig = Figure(; size=(1024, 1024))
	ax  = Axis(fig[1, 1];
	           xlabel = "theta2 [rad]",
	           ylabel = "Semi-major axis a  [normalized]",
	           title  = @sprintf("Moon-periapsis Poincare map  (C = %.2f, random orbits)", C_val))

	scatter!(ax, mod2pi.(df.theta2), df.a; color=:black, markersize=1.5)

	xlims!(ax, 0.0, 2.0pi)
	ylims!(ax, 0.2, 0.8)

	outpath = joinpath(results_dir, @sprintf("poincare_moon_C=%.2f_random.png", C_val))
	save(outpath, fig)
	println("  Saved: $(basename(outpath))")
end

println("Done.")
