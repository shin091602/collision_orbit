using DataFrames
using Arrow
using Printf

const mu = 0.01215058560962404
const x_moon = 1.0 - mu
const L_km = 384400.0     # 地球-月間の平均距離 (km)
const R_moon_km = 1737.4    # 月の半径 (km)

target_altitudes_km = [0.0, 50.0, 100.0, 200.0]
target_r_norm = (R_moon_km .+ target_altitudes_km) ./ L_km  # 月からの物理距離（正規化）

base_dir = "data/CL_ver2/planar"
ts_dir = joinpath(base_dir, "timeseries")
meta_path = joinpath(base_dir, "metadata_all.arrow")

df_meta = DataFrame(Arrow.Table(meta_path))

for alt in target_altitudes_km
    alt_int = Int(alt)
    df_meta[!, Symbol("E_kep_$(alt_int)km")]   = Vector{Union{Missing, Float64}}(missing, nrow(df_meta))
    # df_meta[!, Symbol("v_in_$(alt_int)km")]   = Vector{Union{Missing, Float64}}(missing, nrow(df_meta))
    df_meta[!, Symbol("v_rot_$(alt_int)km")]   = Vector{Union{Missing, Float64}}(missing, nrow(df_meta))
    df_meta[!, Symbol("theta_$(alt_int)km")]   = Vector{Union{Missing, Float64}}(missing, nrow(df_meta))
    # df_meta[!, Symbol("x_$(alt_int)km")]       = Vector{Union{Missing, Float64}}(missing, nrow(df_meta))
    # df_meta[!, Symbol("y_$(alt_int)km")]       = Vector{Union{Missing, Float64}}(missing, nrow(df_meta))
    # df_meta[!, Symbol("time_s_$(alt_int)km")]  = Vector{Union{Missing, Float64}}(missing, nrow(df_meta)) 
end

orbit_id_to_row = Dict(id => i for (i, id) in enumerate(df_meta.orbit_id))

C_values = unique(df_meta.jacobi_constant)

for C_j in C_values
    println("\nJacobi constant C = $(C_j) の解析を開始します。")
    ts_path = joinpath(ts_dir, "timeseries_C_$(C_j).arrow")
    df_time = DataFrame(Arrow.Table(ts_path), copycols=false)
    gdf = groupby(df_time, :orbit_id)
    orbits = collect(gdf)

    Threads.@threads for df_orb in orbits
        orbit_id = df_orb.orbit_id[1]
        row_idx = get(orbit_id_to_row, orbit_id, nothing)
        if row_idx === nothing; continue; end
        r_array = df_orb.r  # 時系列に保存されている月中心からの物理距離 r_M = u1²+u2²

        for (i, r_target) in enumerate(target_r_norm)
            alt_int = Int(target_altitudes_km[i])

            k = findfirst(>(r_target), r_array)

            if k !== nothing && k > 1
                r_prev = r_array[k-1]
                r_next = r_array[k]
                fraction = (r_target - r_prev) / (r_next - r_prev)
                # 各物理量を補間
                px_val    = df_orb.px[k-1]    + fraction * (df_orb.px[k]    - df_orb.px[k-1])
                py_val    = df_orb.py[k-1]    + fraction * (df_orb.py[k]    - df_orb.py[k-1])
                theta_val = df_orb.theta[k-1] + fraction * (df_orb.theta[k] - df_orb.theta[k-1])
                x_val     = df_orb.x[k-1]     + fraction * (df_orb.x[k]     - df_orb.x[k-1])
                y_val     = df_orb.y[k-1]     + fraction * (df_orb.y[k]     - df_orb.y[k-1])
                s_val     = df_orb.t[k-1]     + fraction * (df_orb.t[k]     - df_orb.t[k-1])

                vx_rot_val = px_val + y_val  # 回転座標系での速度成分
                vy_rot_val = py_val - x_val
                v_rot_sq = vx_rot_val^2 + vy_rot_val^2
                v_rot    = sqrt(v_rot_sq)
                
                v_in_sq = px_val^2 + (py_val + mu - 1)^2
                E_kep    = 0.5 * v_in_sq - (mu / r_target)
                
                df_meta[row_idx, Symbol("E_kep_$(alt_int)km")]  = E_kep
                # df_meta[row_idx, Symbol("v_in_$(alt_int)km")]  = v_in
                df_meta[row_idx, Symbol("v_rot_$(alt_int)km")]  = v_rot
                df_meta[row_idx, Symbol("theta_$(alt_int)km")]  = theta_val
                # df_meta[row_idx, Symbol("x_$(alt_int)km")]      = x_val
                # df_meta[row_idx, Symbol("y_$(alt_int)km")]      = y_val
                # df_meta[row_idx, Symbol("time_s_$(alt_int)km")] = abs(s_val)
            end
        end
    end
end

output_meta_path = joinpath(base_dir, "metadata_analyzed.arrow")
Arrow.write(output_meta_path, df_meta)
println("\nすべての解析が完了しました！")
println("新しいメタデータ保存先: $(output_meta_path)")