function _packet_celerity(p::ReachParameters, i::Int, qchar::Float64)
    c = _reach_celerity(p, i, max(qchar, 1.0e-12))
    if !(c > 0.0) || !isfinite(c)
        # A dry/flat reach should not create infinite travel time. Use the
        # user-specified IRF velocity as a conservative fallback.
        return p.irf_velocity[i]
    end
    return c
end

function _merge_packets!(packets::Vector{WavePacket{Float64}}, method::LagrangianKWT)
    length(packets) <= 1 && return packets
    merged = WavePacket{Float64}[]
    push!(merged, packets[1])
    for k in 2:length(packets)
        a = merged[end]
        b = packets[k]
        rel = abs(a.qchar - b.qchar) / max(abs(a.qchar), abs(b.qchar), 1.0e-12)
        if rel <= method.merge_rtol
            total = a.volume + b.volume
            if total > 0.0
                a.distance = (a.distance * a.volume + b.distance * b.volume) / total
                a.qchar = (a.qchar * a.volume + b.qchar * b.volume) / total
                a.volume = total
            end
        else
            push!(merged, b)
        end
    end
    if length(merged) > method.max_packets
        # Coarsen deterministically by merging the two oldest/closest-to-outlet
        # packets until the cap is met.
        while length(merged) > method.max_packets
            a = merged[1]
            b = merged[2]
            total = a.volume + b.volume
            q = total > 0 ? (a.qchar * a.volume + b.qchar * b.volume) / total : 0.0
            d = total > 0 ? (a.distance * a.volume + b.distance * b.volume) / total : 0.0
            merged[1] = WavePacket(total, d, q)
            deleteat!(merged, 2)
        end
    end
    empty!(packets)
    append!(packets, merged)
    return packets
end

"""Lagrangian characteristic routing.

Each time-step inflow is represented by a conservative water packet. Packets
move at kinematic celerity dQ/dA. This preserves the characteristic/Lagrangian
interpretation of mizuRoute KWT while replacing the original Fortran wave-array
bookkeeping with a compact Julia representation.
"""
function _route_reach!(method::LagrangianKWT, state::KWTState, p::ReachParameters,
    i::Int, qin::Float64, qlat::Float64, dt::Float64)
    packets = state.packets[i]
    outvol = 0.0

    # Advect existing characteristics.
    keep = WavePacket{Float64}[]
    for packet in packets
        c = _packet_celerity(p, i, packet.qchar)
        packet.distance -= c * dt
        if packet.distance <= 0.0
            outvol += packet.volume
        else
            push!(keep, packet)
        end
    end
    empty!(packets)
    append!(packets, keep)

    # Add the current characteristic. Lateral inflow is treated as reach input
    # in the compact core; the Euler/DW solvers distribute it spatially.
    qnew = max(0.0, qin + qlat)
    vin = qnew * dt
    if vin > 0.0
        cnew = _packet_celerity(p, i, qnew)
        dist = p.length[i] - cnew * dt
        if dist <= 0.0
            outvol += vin
        else
            push!(packets, WavePacket(vin, dist, qnew))
        end
    end
    _merge_packets!(packets, method)

    qcandidate = outvol / dt
    qout, vnew = _balance_update(state.volume[i], qin, qlat, qcandidate, dt)
    # If balance capping occurred, retain the difference in reach storage.
    state.qout[i] = qout
    state.volume[i] = vnew
    return qout
end
