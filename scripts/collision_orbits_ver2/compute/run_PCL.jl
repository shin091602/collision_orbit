using DifferentialEquations
using DataFrames
using Arrow
using Printf

include("../../../src/lc_canonical_cr3bp.jl")
include("../../../src/utils_regularization.jl")
include("utils.jl")

mu = 0.01215058560962404

C_values = [2.94; 3.16]
N_orbits = 2000
N_angles = 18
tspan  = 10pi

base_dir = joinpath(@__DIR__, "../../../data/CL_ver2", "planar")
ts_dir = joinpath(base_dir, "timeseries")
mkpath(ts_dir)

all_meta_data = []
w0_grid = generate_planar_ks_collision_grid(N_orbits, mu)

for C_j in C_values
    println("Jacobi constant C = $(C_j) の計算を開始...")
    # Cartesian(物理)空間の配列
    t_all  = Float64[]
    x_all  = Float64[]
    y_all  = Float64[]
    vx_all = Float64[]
    vy_all = Float64[]
    # LC(正則化)空間の配列を追加
    u1_all = Float64[]
    u2_all = Float64[]
    w1_all = Float64[]
    w2_all = Float64[]

    r_all = Float64[]
    
    id_all = String[]

    for i in 1:N_orbits
        orbit_id = @sprintf("%.3f_%d", C_j, i)
        w0 = w0_grid[i]
        zeta0 = [0.0, 0.0, w0[1], w0[2]]
        prob = ODEProblem(lc_canonical_cr3bp, zeta0, (0.0, -tspan), [mu, C_j])
        sol  = solve(prob, Vern9(); abstol = 1e-13, reltol = 1e-13, maxiters = Int(1e6))
        num_points = length(sol.t)
        traj_cart = zeros(4, num_points)
        traj_lc = zeros(4, num_points)
        traj_r = zeros(num_points)
        for j in 1:num_points
            r = sqrt(sol.u[j][1]^2 + sol.u[j][2]^2)
            traj_r[j] = r
            traj_cart[:, j] = lc2cart(sol.u[j], mu)
            traj_lc[:, j] = sol.u[j]
        end

        push!(all_meta_data,(
            orbit_id = orbit_id,
            jacobi_constant = C_j,
            u1_0 = zeta0[1],
            u2_0 = zeta0[2],
            w1_0 = zeta0[3],
            w2_0 = zeta0[4],
            status = string(sol.retcode)
        ))
        # データの結合
        append!(id_all, fill(orbit_id, num_points))
        append!(t_all, sol.t)

        append!(x_all, traj_cart[1, :])
        append!(y_all, traj_cart[2, :])
        append!(vx_all, traj_cart[3, :])
        append!(vy_all, traj_cart[4, :])
        
        append!(u1_all, traj_lc[1, :])
        append!(u2_all, traj_lc[2, :])
        append!(w1_all, traj_lc[3, :])
        append!(w2_all, traj_lc[4, :])

        append!(r_all, traj_r)
    end
    df_time = DataFrame(
        orbit_id = id_all,
        t = t_all,
        x = x_all,
        y = y_all,
        vx = vx_all,
        vy = vy_all,
        u1 = u1_all,
        u2 = u2_all,
        w1 = w1_all,
        w2 = w2_all,
        r = r_all
    )
    file_path = joinpath(ts_dir, "timeseries_C_$(C_j).arrow")
    Arrow.write(file_path, df_time)
    println("Jacobi constant C = $(C_j) の計算が完了。データを保存しました: $(file_path)")

    df_time = nothing
    t_all = x_all = y_all = vx_all = vy_all = u1_all = u2_all = w1_all = w2_all = r_all = nothing
    id_all = nothing
    GC.gc()
end

meta_path = joinpath(base_dir, "metadata_all.arrow")
df_meta_all = DataFrame(all_meta_data)
Arrow.write(meta_path, df_meta_all)
println("メタデータを保存しました: $(meta_path)")
