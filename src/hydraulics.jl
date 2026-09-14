const _HYD_EPS = 1.0e-12

"""Top width of a compound trapezoidal channel [m]."""
function top_width(y::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=Inf)
    y1 = max(Float64(y), 0.0)
    b1 = Float64(b)
    zc1 = Float64(zc)
    bf = Float64(bankfull_depth)
    zf1 = Float64(zf)
    if y1 <= bf
        return b1 + 2.0 * zc1 * y1
    end
    bbank = b1 + 2.0 * zc1 * bf
    return bbank + 2.0 * zf1 * (y1 - bf)
end

"""Wetted perimeter [m] of the compound cross section."""
function wetted_perimeter(y::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=Inf)
    y1 = max(Float64(y), 0.0)
    b1 = Float64(b)
    zc1 = Float64(zc)
    bf = Float64(bankfull_depth)
    zf1 = Float64(zf)
    if y1 <= bf
        return b1 + 2.0 * y1 * hypot(1.0, zc1)
    end
    pbank = b1 + 2.0 * bf * hypot(1.0, zc1)
    return pbank + 2.0 * (y1 - bf) * hypot(1.0, zf1)
end

"""Flow cross-sectional area [m²]."""
function flow_area(y::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=Inf)
    y1 = max(Float64(y), 0.0)
    b1 = Float64(b)
    zc1 = Float64(zc)
    bf = Float64(bankfull_depth)
    zf1 = Float64(zf)
    if y1 <= bf
        return y1 * (b1 + zc1 * y1)
    end
    abank = bf * (b1 + zc1 * bf)
    bbank = b1 + 2.0 * zc1 * bf
    dy = y1 - bf
    return abank + dy * (bbank + zf1 * dy)
end

"""Invert cross-sectional area to water depth [m]."""
function water_height(area::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=Inf)
    a = max(Float64(area), 0.0)
    b1 = Float64(b)
    zc1 = Float64(zc)
    bf = Float64(bankfull_depth)
    zf1 = Float64(zf)
    if isinf(bf)
        if zc1 == 0.0
            return a / b1
        end
        return (-b1 + sqrt(max(0.0, b1^2 + 4.0 * zc1 * a))) / (2.0 * zc1)
    end
    abank = flow_area(bf, b1, zc1; zf=zf1, bankfull_depth=bf)
    if a <= abank
        if zc1 == 0.0
            return a / b1
        end
        return (-b1 + sqrt(max(0.0, b1^2 + 4.0 * zc1 * a))) / (2.0 * zc1)
    end
    bbank = top_width(bf, b1, zc1; zf=zf1, bankfull_depth=bf)
    if zf1 == 0.0
        return bf + (a - abank) / bbank
    end
    disc = bbank^2 + 4.0 * zf1 * (a - abank)
    return bf + (-bbank + sqrt(max(0.0, disc))) / (2.0 * zf1)
end

"""Hydraulic radius A/P [m]."""
function hydraulic_radius(y::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=Inf)
    a = flow_area(y, b, zc; zf=zf, bankfull_depth=bankfull_depth)
    p = wetted_perimeter(y, b, zc; zf=zf, bankfull_depth=bankfull_depth)
    return a / max(p, _HYD_EPS)
end

"""Manning uniform-flow discharge [m³/s]."""
function manning_discharge(y::Real, b::Real, zc::Real, slope::Real, mann_n::Real;
    zf::Real=1000.0, bankfull_depth::Real=Inf)
    y <= 0 && return 0.0
    a = flow_area(y, b, zc; zf=zf, bankfull_depth=bankfull_depth)
    r = hydraulic_radius(y, b, zc; zf=zf, bankfull_depth=bankfull_depth)
    s = max(Float64(slope), 0.0)
    n = Float64(mann_n)
    return (a * r^(2.0 / 3.0) * sqrt(s)) / n
end

"""Normal flow depth [m] obtained by robust bracketed bisection of Manning's equation."""
function flow_depth(q::Real, b::Real, zc::Real, slope::Real, mann_n::Real;
    zf::Real=1000.0, bankfull_depth::Real=Inf, rtol::Real=1.0e-8, maxiter::Integer=100)
    target = max(Float64(q), 0.0)
    target == 0.0 && return 0.0
    slope <= 0 && return 0.0
    lo = 0.0
    hi = isfinite(bankfull_depth) ? max(1.0, Float64(bankfull_depth)) : 1.0
    while manning_discharge(hi, b, zc, slope, mann_n; zf=zf, bankfull_depth=bankfull_depth) < target
        hi *= 2.0
        hi > 1.0e6 && throw(ArgumentError("could not bracket flow depth for q=$target"))
    end
    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        qm = manning_discharge(mid, b, zc, slope, mann_n; zf=zf, bankfull_depth=bankfull_depth)
        abs(qm - target) <= max(1.0, target) * rtol && return mid
        if qm < target
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

"""Channel storage [m³] for uniform depth over a reach."""
storage(y::Real, length::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=Inf) =
    flow_area(y, b, zc; zf=zf, bankfull_depth=bankfull_depth) * Float64(length)

"""Kinematic-wave celerity dQ/dA [m/s], evaluated numerically from Manning flow."""
function celerity(q::Real, b::Real, zc::Real, slope::Real, mann_n::Real;
    zf::Real=1000.0, bankfull_depth::Real=Inf)
    q1 = max(Float64(q), 0.0)
    q1 <= _HYD_EPS && return 0.0
    y = flow_depth(q1, b, zc, slope, mann_n; zf=zf, bankfull_depth=bankfull_depth)
    dy = max(1.0e-6, 1.0e-4 * max(y, 1.0))
    y0 = max(0.0, y - dy)
    y2 = y + dy
    q0 = manning_discharge(y0, b, zc, slope, mann_n; zf=zf, bankfull_depth=bankfull_depth)
    q2 = manning_discharge(y2, b, zc, slope, mann_n; zf=zf, bankfull_depth=bankfull_depth)
    a0 = flow_area(y0, b, zc; zf=zf, bankfull_depth=bankfull_depth)
    a2 = flow_area(y2, b, zc; zf=zf, bankfull_depth=bankfull_depth)
    return max(0.0, (q2 - q0) / max(a2 - a0, _HYD_EPS))
end

"""Diffusive-wave diffusivity [m²/s].

Using D = K²/(2QB) and Q = K√S gives D = Q/(2BS) under the
uniform-flow approximation used for routing parameters.
"""
function diffusivity(q::Real, b::Real, zc::Real, slope::Real, mann_n::Real;
    zf::Real=1000.0, bankfull_depth::Real=Inf)
    q1 = max(Float64(q), 0.0)
    q1 <= _HYD_EPS && return 0.0
    s = max(Float64(slope), 1.0e-10)
    y = flow_depth(q1, b, zc, slope, mann_n; zf=zf, bankfull_depth=bankfull_depth)
    bt = top_width(y, b, zc; zf=zf, bankfull_depth=bankfull_depth)
    return q1 / (2.0 * max(bt, _HYD_EPS) * s)
end

@inline function _reach_celerity(p::ReachParameters, i::Int, q::Real)
    celerity(q, p.bottom_width[i], p.side_slope[i], p.slope[i], p.mann_n[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
end

@inline function _reach_diffusivity(p::ReachParameters, i::Int, q::Real)
    diffusivity(q, p.bottom_width[i], p.side_slope[i], p.slope[i], p.mann_n[i];
        zf=p.floodplain_slope[i], bankfull_depth=p.bankfull_depth[i])
end
