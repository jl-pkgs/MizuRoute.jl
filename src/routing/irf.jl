@inline function _irf_green(x::Float64, t::Float64, c::Float64, d::Float64)
    t <= 0.0 && return 0.0
    d <= 0.0 && return 0.0
    den = 2.0 * t * sqrt(pi * t * d)
    return x / den * exp(-((c * t - x)^2) / (4.0 * d * t))
end

"""Build a discrete reach IRF from the mizuRoute diffusive-wave Green function."""
function _irf_weights(length::Real, velocity::Real, diff::Real, dt::Real;
    horizon_factor::Real=6.0, min_steps::Integer=8)
    x = Float64(length)
    c = Float64(velocity)
    d = Float64(diff)
    Δt = Float64(dt)
    mean_t = x / max(c, 1.0e-12)
    n = max(Int(min_steps), ceil(Int, Float64(horizon_factor) * max(mean_t, Δt) / Δt))
    w = zeros(Float64, n)
    if d <= 1.0e-14
        lag = mean_t / Δt
        k0 = clamp(floor(Int, lag) + 1, 1, n)
        frac = clamp(lag - floor(lag), 0.0, 1.0)
        w[k0] += 1.0 - frac
        if k0 < n
            w[k0 + 1] += frac
        else
            w[k0] += frac
        end
    else
        for k in 1:n
            t = (k - 0.5) * Δt
            w[k] = _irf_green(x, t, c, d) * Δt
        end
    end
    s = sum(w)
    if !(s > 0.0) || !isfinite(s)
        fill!(w, 0.0)
        k = clamp(round(Int, mean_t / Δt) + 1, 1, n)
        w[k] = 1.0
    else
        w ./= s
    end
    return w
end

function _route_reach!(::IRF, state::IRFState, p::ReachParameters,
    i::Int, qin::Float64, qlat::Float64, dt::Float64)
    input = max(0.0, qin + qlat)
    h = state.history[i]
    pop!(h)
    pushfirst!(h, input)
    qcandidate = sum(state.kernel[i] .* h)
    qout, vnew = _balance_update(state.volume[i], qin, qlat, qcandidate, dt)
    state.qout[i] = qout
    state.volume[i] = vnew
    return qout
end
