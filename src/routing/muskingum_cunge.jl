@inline function _mc_x_celerity(p::ReachParameters, i::Int, qbar::Float64)
    q = abs(qbar)
    depth = flow_depth(q, p.bottom_width[i], p.side_slope[i], p.slope[i], p.mann_n[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
    btop = top_width(depth, p.bottom_width[i], p.side_slope[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
    ck = celerity(q, p.bottom_width[i], p.side_slope[i], p.slope[i], p.mann_n[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
    s = p.slope[i]
    x = (q > 1.0e-50 && ck > 0.0 && s > 0.0) ?
        0.5 * (1.0 - q / (btop * s * ck * p.length[i])) : 0.0
    return x, ck
end

"""Muskingum-Cunge update matching mizuRoute v3.1 `mc_route.f90`."""
function _route_reach!(::MuskingumCunge, state::MCState,
    p::ReachParameters, i::Int, qin::Float64, qlat::Float64, dt::Float64)
    iprev = state.input_prev[i]
    oprev = state.output_prev[i]
    qbar0 = (iprev + qin + oprev) / 3.0

    channel_out = 0.0
    if qbar0 > 1.0e-50
        _, ck0 = _mc_x_celerity(p, i, qbar0)
        cn0 = ck0 * dt / p.length[i]
        nsub = cn0 > 1.0 ? ceil(Int, cn0) : 1
        dts = dt / nsub
        qouts = zeros(Float64, nsub)
        oldout = oprev
        oldin = iprev
        for k in 1:nsub
            newin = qin
            qbar = (newin + oldin + oldout) / 3.0
            if qbar > 1.0e-50
                x, ck = _mc_x_celerity(p, i, qbar)
                cn = ck * dts / p.length[i]
                den = 1.0 - x + 0.5 * cn
                c0 = (-x + 0.5 * cn) / den
                c1 = ( x + 0.5 * cn) / den
                c2 = (1.0 - x - 0.5 * cn) / den
                qouts[k] = max(0.0, c0 * newin + c1 * oldin + c2 * oldout)
            else
                qouts[k] = 0.0
            end
            oldout = qouts[k]
            oldin = newin
        end
        channel_out = sum(qouts) / nsub
        if abs(channel_out) > 0.0
            reduction = min((state.volume[i] / dt + qin) * 0.999 / channel_out, 1.0)
            channel_out *= reduction
        end
    end

    state.volume[i] += (qin - channel_out) * dt
    state.input_prev[i] = qin
    state.output_prev[i] = channel_out
    state.qout[i] = max(0.0, channel_out + qlat)
    return state.qout[i]
end