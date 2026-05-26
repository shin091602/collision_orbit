using DifferentialEquations
using Plots

# Lorenz system
function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - β * u[3]
end

u0 = [1.0, 0.0, 0.0]
p  = (10.0, 28.0, 8/3)
tspan = (0.0, 50.0)

prob = ODEProblem(lorenz!, u0, tspan, p)
sol  = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

println("Solved: $(length(sol.t)) time steps, t ∈ [$(sol.t[1]), $(sol.t[end])]")
println("Final state: $(sol.u[end])")

# 3D phase portrait
plt = plot(sol, idxs=(1,2,3),
           xlabel="x", ylabel="y", zlabel="z",
           title="Lorenz Attractor",
           legend=false, linewidth=0.5, color=:plasma)

savefig(plt, "lorenz_attractor.png")
println("Saved: lorenz_attractor.png")
