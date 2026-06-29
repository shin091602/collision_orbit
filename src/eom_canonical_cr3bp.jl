using LinearAlgebra
using ForwardDiff

function H_hamiltonian(q, p, mu)
    r1 = sqrt((mu + q[1])^2 + q[2]^2 + q[3]^2)
    r2 = sqrt((q[1] - 1 + mu)^2 + q[2]^2 + q[3]^2)

    H = (1/2)*(p[1]^2 + p[2]^2 + p[3]^2) + p[1]*q[2] - p[2]*q[1] - (1 - mu)/r1 - mu/r2
    return H
end


function eom_hamiltonian_normal(dz, z, mu, t)
    q = z[1:3]
    p = z[4:6]
    dH_dp = ForwardDiff.gradient(p -> H_hamiltonian(q, p, mu), p)
    dH_dq = ForwardDiff.gradient(q -> H_hamiltonian(q, p, mu), q)

    dz[1:3] = dH_dp
    dz[4:6] = -dH_dq

    return nothing
end