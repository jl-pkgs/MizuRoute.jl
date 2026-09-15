@inline function _irf_green(x::Float64, t::Float64, c::Float64, d::Float64)
    (t <= 0.0 || d <= 0.0 || c <= 0.0) && return 0.0
    pot = ((c * t - x)^2) / (4.0 * d * t)
    pot > 69.0 && return 0.0
    return x / (2.0 * sqrt(pi * d * t)) * exp(-pot)
end

"""Build the reach IRF exactly as mizuRoute `make_uh` does.

The Green function is sampled at one-hour resolution for at most 240 hours,
normalized, convolved with a rectangular runoff pulse whose duration is the
routing time step, truncated at cumulative probability 0.9999, and finally
aggregated to the model time step.
"""
function _irf_weights(length::Real, velocity::Real, diff::Real, dt::Real;
    horizon_factor::Real=6.0, min_steps::Integer=8)
    x = Float64(length)
    c = Float64(velocity)
    d = Float64(diff)
    Δt = Float64(dt)
    dtu = 3600.0
    ntmax = 240
    ntsub = ceil(Int, Δt / dtu)

    uhm = zeros(Float64, ntmax)
    if c > 0.0 && d > 0.0
        for ihr in 1:ntmax
            uhm[ihr] = _irf_green(x, ihr * dtu, c, d)
        end
        s = sum(uhm)
        s > 0.0 && (uhm ./= s)
    end

    cum = 0.0
    ihr_last = ntmax
    for ihr in 1:ntmax
        cum += uhm[ihr]
        ihr_last = ihr
        cum > 0.99999 && break
    end
    cum = 0.0
    ihr_start = ntmax
    for ihr in ntmax:-1:1
        cum += uhm[ihr]
        ihr_start = ihr
        cum > 0.99999 && break
    end

    fr = zeros(Float64, ntmax)
    fr[1:min(ntsub, ntmax)] .= 1.0 / ntsub
    uhq = zeros(Float64, ntmax)
    for jhr in 1:ntmax
        q0 = 0.0
        for ihr in ihr_start:ihr_last
            lag = jhr - ihr
            if lag > 0
                lag <= ntsub && (q0 += fr[lag] * uhm[ihr])
            else
                break
            end
        end
        uhq[jhr] = q0
    end
    s = sum(uhq)
    if !(s > 0.0)
        n = max(Int(min_steps), 1)
        w = zeros(Float64, n)
        w[clamp(round(Int, x / max(c * Δt, 1.0e-12)) + 1, 1, n)] = 1.0
        return w
    end
    uhq ./= s

    cum = 0.0
    ihr_last = ntmax
    for ihr in 1:ntmax
        cum += uhq[ihr]
        ihr_last = ihr
        cum > 0.9999 && break
    end
    uhq ./= cum

    ntdh = div(ihr_last + ntsub - 1, ntsub)
    seguh = zeros(Float64, ntdh)
    for jhr in 1:ihr_last
        itagg = div(jhr + ntsub - 1, ntsub)
        seguh[itagg] += uhq[jhr]
    end
    return seguh
end

"""mizuRoute v3.1 IRF reach update.

Only upstream discharge is delayed by the reach IRF. Local lateral flow is
inserted at the bottom of the reach (the default `hw_drain_point=2`) and is not
included in channel storage.
"""
function _route_reach!(::IRF, state::IRFState, p::ReachParameters,
    i::Int, qin::Float64, qlat::Float64, dt::Float64)
    future = state.history[i]
    kernel = state.kernel[i]
    @inbounds for k in eachindex(kernel)
        future[k] += kernel[k] * qin
    end

    channel_out = min((max(0.0, state.volume[i]) / dt + qin) * 0.999, future[1])
    state.volume[i] -= (channel_out - qin) * dt
    qout = channel_out + qlat

    @inbounds for k in 1:(length(future) - 1)
        future[k] = future[k + 1]
    end
    future[end] = 0.0

    state.qout[i] = max(0.0, qout)
    return state.qout[i]
end