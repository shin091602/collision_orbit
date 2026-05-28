using LinearAlgebra

function cart2lc(z, mu)
    qM = [z[1] - 1 + mu, z[2]]
    p = [z[3], z[4]]
    r_M = sqrt(qM[1]^2 + qM[2]^2)
    u = [sqrt((qM[1] + r_M) / 2), qM[2] / (2 * sqrt((qM[1] + r_M) / 2))]

    L = hcat([u[1], u[2]], [-u[2], u[1]])
    w = 2 * L' * p

    return vcat(u, w)
end

function lc2cart(zeta, mu)
    u = zeta[1:2]
    w = zeta[3:4]
    L = hcat([u[1], u[2]], [-u[2], u[1]])

    q = [u[1]^2 - u[2]^2 + 1 - mu, 2 * u[1] * u[2]]
    p = (1 / (2*(u[1]^2 + u[2]^2))) * L * w

    z = vcat(q, p)

    return z
end

function cart2ks(z, mu)
    q = z[1:3]
    p = z[4:6]
    r2 = sqrt((q[1] - 1 + mu)^2 + q[2]^2 + q[3]^2)
    u_1 = sqrt((q[1] + r2 + mu - 1) / 2)
    u = [u_1, q[2] / (2*u_1), q[3] / (2*u_1), 0]
    A = hcat(
        [u[1],  u[2],  u[3]],
        [-u[2], u[1],  u[4]],
        [-u[3], -u[4], u[1]],
        [u[4],  -u[3], u[2]],
    )
    w = 2 * A' * p
    zeta = vcat(u, w)
    return zeta
end

function ks2cart(zeta, mu)
    u = zeta[1:4]
    w = zeta[5:8]
    A = hcat(
        [u[1],  u[2],  u[3],  u[4]],
        [-u[2], u[1],  u[4], -u[3]],
        [-u[3], -u[4], u[1],  u[2]],
        [u[4],  -u[3], u[2], -u[1]],
    )
    rM = u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2
    p = (1/(2*rM)) * A * w
    q = A * u + [1-mu, 0.0, 0.0, 0.0]
    z = vcat(q[1:3], p[1:3])
    return z
end