"""One-step water-balance residual [m³].

Positive `abstraction` denotes water removed from the reach; negative values are
injections. A numerically conservative routing step should keep this close to
floating-point roundoff.
"""
function water_balance_error(vold::Real, vnew::Real, qin::Real, qlat::Real,
    qout::Real, dt::Real; abstraction::Real=0.0)
    Float64(vnew) - Float64(vold) -
        (Float64(qin) + Float64(qlat) - Float64(qout) - Float64(abstraction)) * Float64(dt)
end

function _apply_water_management(qin::Float64, qlat::Float64, volume::Float64,
    request::Float64, dt::Float64)
    # mizuRoute convention: positive = abstraction, negative = injection.
    request == 0.0 && return qin, qlat, volume, 0.0
    if request < 0.0
        return qin, qlat - request, volume, request
    end

    remaining = request
    take_storage = min(remaining * dt, volume)
    volume2 = volume - take_storage
    remaining -= take_storage / dt

    take_qin = min(remaining, qin)
    qin2 = qin - take_qin
    remaining -= take_qin

    take_qlat = min(remaining, qlat)
    qlat2 = qlat - take_qlat
    remaining -= take_qlat

    actual = request - remaining
    return qin2, qlat2, volume2, actual
end

@inline function _balance_update(vold::Float64, qin::Float64, qlat::Float64,
    qout_candidate::Float64, dt::Float64)
    available = max(0.0, vold + (qin + qlat) * dt)
    qmax = available / dt
    qout = clamp(qout_candidate, 0.0, qmax)
    vnew = max(0.0, available - qout * dt)
    return qout, vnew
end
