function _dw_reference_hydraulics(p::ReachParameters, i::Int,
    profile::Vector{Float64}, qin::Float64)
    qref = max(qin, maximum(profile), 1.0e-10)
    c = _reach_celerity(p, i, qref)
    d = _reach_diffusivity(p, i, qref)
    if !(c > 0.0) || !isfinite(c)
        c = p.irf_velocity[i]
    end
    if !(d >= 0.0) || !isfinite(d)
        d = p.irf_diffusivity[i]
    end
    return c, d
end

"""Diffusive-wave reach routing using the mizuRoute ADE discretization.

The linearized equation is
`∂Q/∂t + cₖ ∂Q/∂x = D ∂²Q/∂x² + q_lat`, solved with the same tridiagonal
ADE kernel used by Euler-KW. `alpha` and `beta` correspond to the new-time
advection and diffusion weights (`wc`, `wd` in mizuRoute). Upstream discharge
is Dirichlet and the default outlet condition preserves the previous gradient.
"""
function _route_reach!(method::DiffusiveWave, state::DWState,
    p::ReachParameters, i::Int, qin::Float64, qlat::Float64, dt::Float64)
    qold = state.profile[i]
    length(qold) >= 4 || throw(ArgumentError("diffusive-wave profile needs at least 4 nodes"))
    cvel, dcoef = _dw_reference_hydraulics(p, i, qold, qin + qlat)

    qnew = _solve_ade(qold, p.length[i], dt, max(qin, 0.0), cvel, dcoef, qlat;
        advection=:central, downstream=:neumann,
        wc=method.alpha, wd=method.beta, include_lateral=true)
    state.profile[i] = qnew

    qcandidate = qnew[end]
    qout, vnew = _balance_update(state.volume[i], qin, qlat, qcandidate, dt)
    state.qout[i] = qout
    state.volume[i] = vnew
    return qout
end
