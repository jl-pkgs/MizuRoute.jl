function _mc_x_celerity(p::ReachParameters, i::Int, qbar::Float64)
    q = max(qbar, 1.0e-12)
    depth = flow_depth(q, p.bottom_width[i], p.side_slope[i], p.slope[i], p.mann_n[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
    btop = top_width(depth, p.bottom_width[i], p.side_slope[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
    ck = _reach_celerity(p, i, q)
    if !(ck > 0.0) || !isfinite(ck)
        ck = p.irf_velocity[i]
    end
    s = max(p.slope[i], 1.0e-10)
    x = 0.5 * (1.0 - q / max(btop * s * ck * p.length[i], 1.0e-12))
    # Classical Cunge weighting should remain in [0, 0.5] for a monotone
    # three-point Muskingum solution.
    return clamp(x, 0.0, 0.5), ck
end

"""Dynamic Muskingum-Cunge step for one reach.

Uses the mizuRoute technical-note coefficients
O(t+1)=C0 I(t+1)+C1 I(t)+C2 O(t), with Y=0.5 and Cunge's dynamic X.
The time step is subdivided until the Courant number is below `method.cfl`.
"""
function _route_reach!(method::MuskingumCunge, state::MCState,
    p::ReachParameters, i::Int, qin::Float64, qlat::Float64, dt::Float64)
    itarget = max(0.0, qin + qlat)
    iprev = max(0.0, state.input_prev[i])
    oprev = max(0.0, state.output_prev[i])
    qbar0 = max((itarget + iprev + oprev) / 3.0, 1.0e-10)
    _, ck0 = _mc_x_celerity(p, i, qbar0)
    nsub = max(1, ceil(Int, ck0 * dt / max(method.cfl * p.length[i], 1.0e-12)))
    dts = dt / nsub

    o = oprev
    iold = iprev
    for k in 1:nsub
        f = k / nsub
        inew = iprev + f * (itarget - iprev)
        qbar = max((inew + iold + o) / 3.0, 1.0e-10)
        x, ck = _mc_x_celerity(p, i, qbar)
        cn = ck * dts / p.length[i]
        den = 1.0 - x + 0.5 * cn
        c0 = (-x + 0.5 * cn) / den
        c1 = ( x + 0.5 * cn) / den
        c2 = (1.0 - x - 0.5 * cn) / den
        onew = c0 * inew + c1 * iold + c2 * o
        o = max(0.0, isfinite(onew) ? onew : 0.0)
        iold = inew
    end

    qout, vnew = _balance_update(state.volume[i], qin, qlat, o, dt)
    state.qout[i] = qout
    state.volume[i] = vnew
    state.input_prev[i] = itarget
    state.output_prev[i] = qout
    return qout
end
