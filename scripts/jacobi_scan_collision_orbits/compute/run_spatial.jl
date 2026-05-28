using DifferentialEquations
using JLD2

include("../../../src/ks_canonical_cr3bp.jl")
include("../../../src/utils_regularization.jl")

const mu = 0.01215058560962404

# Collision constraint: u = 0, |w|^2 = 8μ (KS analogue of the LC collision sphere)
const W_COLL = sqrt(8 * mu)

C_values = [2.8; 2.9; collect(2.94:0.01:3.20)]
N_ETA = 2
N_PHI = 3
N_PSI = 3
tspan  = 20pi   # fictitious time for backward propagation

cb_escape = ContinuousCallback((u, _, _) -> u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2 - 5, terminate!)

function collision_directions(n_eta, n_phi, n_psi)
    dirs = Vector{Vector{Float64}}()
    sizehint!(dirs, n_eta * n_phi * n_psi)

    for i in 1:n_eta
        # For S^3 in Hopf coordinates, ρ₁² is uniform on [0, 1].
        rho1 = sqrt((i - 0.5) / n_eta)
        rho2 = sqrt(1 - rho1^2)
        for j in 0:n_phi - 1
            phi = 2π * j / n_phi
            for k in 0:n_psi - 1
                psi = 2π * k / n_psi
                push!(dirs, [
                    rho1 * cos(phi),
                    rho1 * sin(phi),
                    rho2 * cos(psi),
                    rho2 * sin(psi),
                ])
            end
        end
    end

    return dirs
end

directions = collision_directions(N_ETA, N_PHI, N_PSI)

buf_cart = Vector{Vector{Matrix{Float64}}}(undef, length(C_values))
buf_ks   = Vector{Vector{Matrix{Float64}}}(undef, length(C_values))

Threads.@threads for j in eachindex(C_values)
    C = C_values[j]
    trajs_cart = Vector{Matrix{Float64}}()
    trajs_ks   = Vector{Matrix{Float64}}()

    for dir in directions
        zeta0 = vcat(zeros(4), W_COLL .* dir)
        prob = ODEProblem(ks_canonical_cr3bp, zeta0, (0.0, -tspan), [mu, C])
        sol  = solve(prob, Vern9(); abstol = 1e-13, reltol = 1e-13, maxiters = Int(1e6), callback = cb_escape)
        num_points = length(sol.t)
        traj_cart = zeros(6, num_points)
        traj_ks = zeros(8, num_points)
        for i in 1:num_points
            traj_cart[:, i] = ks2cart(sol.u[i], mu)
            traj_ks[:, i] = sol.u[i]
        end
        push!(trajs_cart, traj_cart)
        push!(trajs_ks,   traj_ks)
    end

    buf_cart[j] = trajs_cart
    buf_ks[j]   = trajs_ks
    println("C = $C: $(length(trajs_cart)) orbits")
end

all_cart = Dict(C_values[j] => buf_cart[j] for j in eachindex(C_values))
all_ks   = Dict(C_values[j] => buf_ks[j]   for j in eachindex(C_values))

const DATA_DIR = joinpath(@__DIR__, "../../../data")
const DATA_PATH = joinpath(DATA_DIR, "jacobi_scan_spatial.jld2")
mkpath(DATA_DIR)
save(DATA_PATH,
     "cart",      all_cart,
     "ks",        all_ks,
     "C_values",  C_values,
     "mu",        mu,
     "directions", directions,
     "N_ETA",     N_ETA,
     "N_PHI",     N_PHI,
     "N_PSI",     N_PSI)
println("Saved to $(relpath(DATA_PATH, joinpath(@__DIR__, "../../..")))")
