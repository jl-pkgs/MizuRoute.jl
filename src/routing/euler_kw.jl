function _kw_reference_celerity(p::ReachParameters, i::Int, profile::Vector{Float64}, qin::Float64)
    qref = max(qin, maximum(profile), 1.0e-10)
    c = _reach_celerity(p, i, qref)
    if !(c > 0.0) || !isfinite(c)
        c = p.irf_velocity[i]
    end
    return c
end

"""Eulerian kinematic-wave routing following the current mizuRoute ADE path.

mizuRoute's Euler kinematic-wave solver linearizes
`∂Q/∂t + cₖ ∂Q/∂x = 0` and calls the shared implicit advection-diffusion
tridiagonal solver with `D=0`. This implementation uses the same centered,
fully implicit default and preserves mizuRoute's previous-gradient Neumann
condition at the outlet. A uniform lateral source is retained as the physically
intended RHS source term.
"""
function _route_reach!(method::EulerKinematicWave, state::EulerKWState,
    p::ReachParameters, i::Int, qin::Float64, qlat::Float64, dt::Float64)
    q = state.profile[i]
    n = length(q)
    n >= 4 || throw(ArgumentError("Euler KW profile needs at least 4 nodes"))

    # Although the implicit scheme is not CFL-limited, mizuRoute recomputes
    # hydraulic properties in substeps. Retaining CFL-based substepping makes
    # the linearized celerity update more robust for flashy hydrographs.
    dx = p.length[i] / max(n - 2, 1)
    c0 = _kw_reference_celerity(p, i, q, qin + qlat)
    nsub = max(1, ceil(Int, c0 * dt / max(method.cfl * dx, 1.0e-12)))
    dts = dt / nsub

    qin0 = q[1]
    for k in 1:nsub
        f = k / nsub
        qbc = qin0 + f * (max(qin, 0.0) - qin0)
        c = _kw_reference_celerity(p, i, q, qbc + qlat)
        q = _solve_ade(q, p.length[i], dts, qbc, c, 0.0, qlat;
            advection=:central, downstream=:neumann, wc=1.0, wd=1.0,
            include_lateral=true)
    end

    state.profile[i] = q
    qcandidate = q[end]
    qout, vnew = _balance_update(state.volume[i], qin, qlat, qcandidate, dt)
    state.qout[i] = qout
    state.volume[i] = vnew
    return qout
end
