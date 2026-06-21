using DifferentialEquations
using DataFrames
using Arrow
using Printf

include("../../../src/lc_canonical_cr3bp.jl")
include("../../../src/utils_regularization.jl")
include("utils.jl")

const mu = 0.01215058560962404
const x_moon = 1.0 - mu

C_values = 2.9:0.001:3.2
const N_orbits = 2000
const tspan  = 10pi

base_dir = joinpath(@__DIR__, "../../../data/CL_ver2", "planar")
ts_dir = joinpath(base_dir, "timeseries")
mkpath(ts_dir)

all_meta_data = []
w0_grid, phi_grid = generate_planar_ks_collision_grid(N_orbits, mu)

println("スレッド数: $(Threads.nthreads())")

for C_j in C_values
    println("Jacobi constant C = $(C_j) の計算を開始...")

    # 軌道ごとの結果を格納するバッファ（型付き・スレッドセーフにインデックスで書き込む）
    buf_orbit_id  = Vector{String}(undef, N_orbits)
    buf_zeta0     = Vector{Vector{Float64}}(undef, N_orbits)
    buf_status    = Vector{String}(undef, N_orbits)
    buf_t         = Vector{Vector{Float64}}(undef, N_orbits)
    buf_traj_cart = Vector{Matrix{Float64}}(undef, N_orbits)
    buf_traj_lc   = Vector{Matrix{Float64}}(undef, N_orbits)
    buf_traj_r    = Vector{Vector{Float64}}(undef, N_orbits)
    buf_traj_theta = Vector{Vector{Float64}}(undef, N_orbits)

    Threads.@threads for i in 1:N_orbits
        orbit_id = @sprintf("%.3f_%d", C_j, i)
        w0 = w0_grid[i]
        zeta0 = [0.0, 0.0, w0[1], w0[2]]
        prob = ODEProblem(lc_canonical_cr3bp, zeta0, (0.0, -tspan), [mu, C_j])
        sol  = solve(prob, Vern9(); abstol = 1e-13, reltol = 1e-13, maxiters = Int(1e6))
        num_points = length(sol.t)
        traj_cart = zeros(4, num_points)
        traj_lc   = zeros(4, num_points)
        traj_r    = zeros(num_points)
        traj_theta = zeros(num_points)
        for j in 1:num_points
            r = sol.u[j][1]^2 + sol.u[j][2]^2   # r_M = u1²+u2²（月からの物理距離）
            traj_r[j] = r
            traj_cart[:, j] = lc2cart(sol.u[j], mu)
            traj_lc[:, j]   = sol.u[j]
            traj_theta[j] = atan(traj_cart[2, j], traj_cart[1, j] - x_moon)
            
        end
        buf_orbit_id[i]  = orbit_id
        buf_zeta0[i]     = zeta0
        buf_status[i]    = string(sol.retcode)
        buf_t[i]         = sol.t
        buf_traj_cart[i] = traj_cart
        buf_traj_lc[i]   = traj_lc
        buf_traj_r[i]    = traj_r
        buf_traj_theta[i] = traj_theta
    end

    # 並列計算後にシーケンシャルに結合
    t_all  = Float64[]
    x_all  = Float64[]
    y_all  = Float64[]
    px_all = Float64[]
    py_all = Float64[]
    u1_all = Float64[]
    u2_all = Float64[]
    w1_all = Float64[]
    w2_all = Float64[]
    r_all  = Float64[]
    theta_all = Float64[]
    id_all = String[]

    for i in 1:N_orbits
        push!(all_meta_data, (
            orbit_id         = buf_orbit_id[i],
            jacobi_constant  = C_j,
            phi_0            = phi_grid[i],
            u1_0             = buf_zeta0[i][1],
            u2_0             = buf_zeta0[i][2],
            w1_0             = buf_zeta0[i][3],
            w2_0             = buf_zeta0[i][4],
            status           = buf_status[i],
        ))
        n = length(buf_t[i])
        append!(id_all, fill(buf_orbit_id[i], n))
        append!(t_all,  buf_t[i])
        append!(x_all,  buf_traj_cart[i][1, :])
        append!(y_all,  buf_traj_cart[i][2, :])
        append!(px_all, buf_traj_cart[i][3, :])
        append!(py_all, buf_traj_cart[i][4, :])
        append!(u1_all, buf_traj_lc[i][1, :])
        append!(u2_all, buf_traj_lc[i][2, :])
        append!(w1_all, buf_traj_lc[i][3, :])
        append!(w2_all, buf_traj_lc[i][4, :])
        append!(r_all,  buf_traj_r[i])
        append!(theta_all, buf_traj_theta[i])
    end

    df_time = DataFrame(
        orbit_id = id_all,
        t  = t_all,
        x  = x_all,
        y  = y_all,
        px = px_all,
        py = py_all,
        u1 = u1_all,
        u2 = u2_all,
        w1 = w1_all,
        w2 = w2_all,
        r  = r_all,
        theta = theta_all,
    )
    file_path = joinpath(ts_dir, "timeseries_C_$(C_j).arrow")
    Arrow.write(file_path, df_time)
    println("Jacobi constant C = $(C_j) の計算が完了。データを保存しました: $(file_path)")

    orbit_results = nothing
    df_time = nothing
    t_all = x_all = y_all = px_all = py_all = u1_all = u2_all = w1_all = w2_all = r_all = theta_all = nothing
    id_all = nothing
    GC.gc()
end

meta_path = joinpath(base_dir, "metadata_all.arrow")
df_meta_all = DataFrame(all_meta_data)
Arrow.write(meta_path, df_meta_all)
println("メタデータを保存しました: $(meta_path)")
