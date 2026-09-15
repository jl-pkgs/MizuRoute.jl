@inline function _kw_hydraulics(p::ReachParameters, i::Int, qbar::Float64)
    q = abs(qbar)
    depth = flow_depth(q, p.bottom_width[i], p.side_slope[i], p.slope[i], p.mann_n[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
    c = celerity(q, p.bottom_width[i], p.side_slope[i], p.slope[i], p.mann_n[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
    return depth, c
end

"""Euler kinematic-wave update matching mizuRoute v3.1 `kwe_route.f90`."""
function _route_reach!(::EulerKinematicWave, state::EulerKWState,
    p::ReachParameters, i::Int, qin::Float64, qlat::Float64, dt::Float64)
    qprev = state.profile[i]
    n = length(qprev)
    n >= 4 || throw(ArgumentError("Euler KW profile needs at least 4 nodes"))

    qbar = (qin + qprev[1] + qprev[n - 1]) / 3.0
    _, ck = _kw_hydraulics(p, i, qbar)
    qnode = _solve_ade(qprev, p.length[i], dt, qin, ck, 0.0, qlat;
        advection=:central, downstream=:neumann, wc=1.0, wd=1.0,
        include_lateral=false)

    channel_out = qnode[n - 1]
    if abs(channel_out) > 0.0
        vol = max(0.0, state.volume[i])
        reduction = min((vol + dt * qin) * 0.999 / (channel_out * dt), 1.0)
        @inbounds qnode[2:end] .*= reduction
        channel_out = qnode[n - 1]
    end
    state.volume[i] += (qin - channel_out) * dt
    state.profile[i] = qnode
    # Fortran stores Qnode(nMolecule-1)+Qlat directly; do not add a Julia-only
    # non-negative clamp after the ADE/low-flow limiter.
    state.qout[i] = channel_out + qlat
    return state.qout[i]
end
