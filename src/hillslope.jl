# Dependency-free Lanczos approximation for log Γ(x), x > 0.
const _LANCZOS_COEF = (
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
)

function _loggamma_pos(x::Float64)
    x > 0 || throw(ArgumentError("gamma shape must be positive"))
    if x < 0.5
        return log(pi) - log(sinpi(x)) - _loggamma_pos(1.0 - x)
    end
    z = x - 1.0
    a = 0.99999999999980993
    for (k, c) in pairs(_LANCZOS_COEF)
        a += c / (z + k)
    end
    t = z + 7.5
    return 0.5 * log(2pi) + (z + 0.5) * log(t) - t + log(a)
end

@inline function _gamma_pdf(t::Float64, shape::Float64, scale::Float64)
    t <= 0 && return 0.0
    exp((shape - 1.0) * log(t) - t / scale - _loggamma_pos(shape) - shape * log(scale))
end

"""Discrete gamma unit-hydrograph weights.

The continuous gamma PDF is sampled at interval midpoints and then normalized
so the discrete ordinates sum to one. `shape` is dimensionless and `scale` is
in seconds.
"""
function gamma_uh_weights(shape::Real, scale::Real, dt::Real;
    horizon_factor::Real=8.0, min_steps::Integer=8)
    a = Float64(shape)
    θ = Float64(scale)
    Δt = Float64(dt)
    a > 0 || throw(ArgumentError("shape must be positive"))
    θ > 0 || throw(ArgumentError("scale must be positive"))
    Δt > 0 || throw(ArgumentError("dt must be positive"))
    mean_t = a * θ
    n = max(Int(min_steps), ceil(Int, Float64(horizon_factor) * mean_t / Δt))
    w = Vector{Float64}(undef, n)
    for k in 1:n
        t = (k - 0.5) * Δt
        w[k] = _gamma_pdf(t, a, θ) * Δt
    end
    s = sum(w)
    s > 0 || throw(ArgumentError("gamma UH collapsed numerically; adjust dt/shape/scale"))
    w ./= s
    return w
end

mutable struct GammaUHRouter{T<:AbstractFloat}
    kernel::Vector{T}
    history::Vector{Vector{T}}
end

function GammaUHRouter(nreach::Integer, dt::Real; shape::Real=2.5, scale::Real=3600.0,
    horizon_factor::Real=8.0, min_steps::Integer=8)
    k = gamma_uh_weights(shape, scale, dt; horizon_factor=horizon_factor, min_steps=min_steps)
    h = [zeros(Float64, length(k)) for _ in 1:nreach]
    GammaUHRouter(k, h)
end

"""Apply gamma-UH hillslope delay to a vector of instantaneous lateral inflows."""
function route_hillslope!(out::AbstractVector, router::GammaUHRouter, qinst::AbstractVector)
    length(out) == length(router.history) == length(qinst) || throw(ArgumentError("reach count mismatch"))
    for i in eachindex(qinst)
        h = router.history[i]
        pop!(h)
        pushfirst!(h, Float64(qinst[i]))
        out[i] = sum(router.kernel .* h)
    end
    return out
end

function route_hillslope!(router::GammaUHRouter, qinst::AbstractVector)
    out = similar(Float64.(qinst))
    route_hillslope!(out, router, qinst)
end
