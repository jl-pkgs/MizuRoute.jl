@inline function _irf_green(x::Float64, t::Float64, c::Float64, d::Float64)
    (t <= 0.0 || d <= 0.0 || c <= 0.0) && return 0.0
    pot = ((c * t - x)^2) / (4.0 * d * t)
    pot > 69.0 && return 0.0
    # mizuRoute public_var uses this literal value rather than the language
    # intrinsic, so retain it for bit-level parity in UH construction.
    π_mizu = 3.14159265359
    return x / (2.0 * sqrt(π_mizu * d * t)) * exp(-pot)
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
        # Keep explicit scalar accumulation order to mirror the Fortran loop.
        inte = 0.0
        for ihr in 1:ntmax
            uhm[ihr] = _irf_green(x, ihr * dtu, c, d)
            inte += uhm[ihr]
        end
        inte > 0.0 && (uhm ./= inte)
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
    inte = 0.0
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
        inte += q0
    end
    if !(inte > 0.0)
        n = max(Int(min_steps), 1)
        w = zeros(Float64, n)
        w[clamp(round(Int, x / max(c * Δt, 1.0e-12)) + 1, 1, n)] = 1.0
        return w
    end
    uhq ./= inte

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

"""mizuRoute IRF reach update.

Only upstream discharge is delayed by the reach IRF. Local lateral flow is
inserted at the bottom of the reach and is not included in channel storage.

`method.volume_limiter` controls version-specific behavior. mizuRoute v1.2
serial IRF used pure UH convolution with no volume limiter; newer versions add
reach-volume protection.
"""
function _route_reach!(method::IRF, state::IRFState, p::ReachParameters,
    i::Int, qin::Float64, qlat::Float64, dt::Float64)
    future = state.history[i]
    kernel = state.kernel[i]
    @inbounds for k in eachindex(kernel)
        future[k] += kernel[k] * qin
    end

    vold = state.volume[i]
    channel_out = if method.volume_limiter === :none
        future[1]
    elseif method.volume_limiter === :legacy
        min(vold / dt + 0.999 * qin, future[1])
    else
        min((max(0.0, vold) / dt + qin) * 0.999, future[1])
    end
    state.volume[i] = vold + (qin - channel_out) * dt
    qout = channel_out + qlat

    @inbounds for k in 1:(length(future) - 1)
        future[k] = future[k + 1]
    end
    future[end] = 0.0

    # Fortran v1.2 and current mizuRoute do not clamp final IRF discharge here.
    state.qout[i] = qout
    return qout
end
