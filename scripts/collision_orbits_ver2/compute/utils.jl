# 円制限3体問題におけるKS変換を用いた衝突軌道の初期値グリッド生成
function generate_ks_collision_grid(N::Int, mu::Float64)
    # === 使用例 ===
    # grid = generate_ks_collision_grid(1000, 0.01215)
    
    R = sqrt(8.0 * mu)
    w_grid = Vector{NTuple{4, Float64}}(undef, N)
    # 黄金角 (ラジアン)
    golden_angle = π * (3.0 - sqrt(5.0))

    for i in 0:(N-1)
        # 1. フィボナッチ球面上での物理的な射出方向 (X, Y, Z) を計算
        Z = 1.0 - (2.0 * i + 1.0) / N
        r_xy = sqrt(max(0.0, 1.0 - Z^2)) # maxは数値誤差による負の平方根を防ぐため
        theta_golden = i * golden_angle
        X = r_xy * cos(theta_golden)
        Y = r_xy * sin(theta_golden)
        # 2. (X, Y, Z) を KS空間の角度パラメータ (φ, θ) に逆算
        phi   = acos(X) / 2.0        # 第1成分 X を cos(2φ) に対応させる
        theta = atan(Z, Y)           # 残り (Y,Z) を赤道面の (cosθ, sinθ) に
        # 3. KS運動量ベクトル w を構築 (ゲージを w4=0 に固定)
        w1 = R * cos(phi)
        w2 = R * sin(phi) * cos(theta)
        w3 = R * sin(phi) * sin(theta)
        w4 = 0.0

        w_grid[i+1] = (w1, w2, w3, w4)
    end
    
    return w_grid
end



function generate_planar_ks_collision_grid(N::Int, mu::Float64)
    # === 使用例 ===
    # grid = generate_planar_ks_collision_grid(360, 0.01215)    

    R = sqrt(8.0 * mu)    
    w_grid = Vector{NTuple{4, Float64}}(undef, N)
    d_phi = π / N

    for i in 0:(N-1)
        phi = i * d_phi
        # 平面問題なので w3 = w4 = 0
        w1 = R * cos(phi)
        w2 = R * sin(phi)
        w3 = 0.0
        w4 = 0.0

        w_grid[i+1] = (w1, w2, w3, w4)
    end
    
    return w_grid
end



