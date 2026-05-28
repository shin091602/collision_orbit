using DifferentialEquations
using JLD2

include("../../../src/lc_canonical_cr3bp.jl")
include("../../../src/utils_regularization.jl")

const mu = 0.01215058560962404

# Collision constraint: u = 0, |w|² = 8μ (from regularized energy at r_M = 0)
const W_COLL = sqrt(8 * mu)

C_values = [2.8; 2.9; collect(2.94:0.01:3.20)]
N_angles = 100
tspan  = 20pi   # fictitious time for backward propagation

cb_escape = ContinuousCallback((u, _, _) -> u[1]^2 + u[2]^2 - 5, terminate!)

buf_cart = Vector{Vector{Matrix{Float64}}}(undef, length(C_values))
buf_lc   = Vector{Vector{Matrix{Float64}}}(undef, length(C_values))

Threads.@threads for j in eachindex(C_values)
    C = C_values[j]
    trajs_cart = Vector{Matrix{Float64}}()
    trajs_lc   = Vector{Matrix{Float64}}()

    for k in 0:N_angles - 1
        alpha = 2π * k / N_angles
        zeta0 = [0.0, 0.0, W_COLL * cos(alpha), W_COLL * sin(alpha)]
        prob = ODEProblem(lc_canonical_cr3bp, zeta0, (0.0, -tspan), [mu, C])
        sol  = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10, maxiters = Int(1e6),　callback = cb_escape)
        num_points = length(sol.t)
        traj_cart = zeros(4, num_points)
        traj_lc = zeros(4, num_points)
        for i in 1:num_points
            traj_cart[:, i] = lc2cart(sol.u[i], mu)
            traj_lc[:, i] = sol.u[i]
        end
        push!(trajs_cart, traj_cart)
        push!(trajs_lc,   traj_lc)
    end

    buf_cart[j] = trajs_cart
    buf_lc[j]   = trajs_lc
    println("C = $C: $(length(trajs_cart)) orbits")
end

all_cart = Dict(C_values[j] => buf_cart[j] for j in eachindex(C_values))
all_lc   = Dict(C_values[j] => buf_lc[j]   for j in eachindex(C_values))

const DATA_DIR = joinpath(@__DIR__, "../../../data")
mkpath(DATA_DIR)
save(joinpath(DATA_DIR, "jacobi_scan_planar.jld2"),
     "cart",     all_cart,
     "lc",       all_lc,
     "C_values", C_values,
     "mu",       mu)
println("Saved → data/jacobi_scan/jacobi_scan_planar.jld2")
