const _HYD_EPS = 1.0e-50
const _CONST23 = 2.0 / 3.0
const _CONST53 = 5.0 / 3.0
const _CONST103 = 10.0 / 3.0
const _FLOW_DEPTH_RTOL = 0.005

"""Top width of the mizuRoute compound trapezoidal section [m]."""
function top_width(y::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=1.0e6)
    yy = max(Float64(y), 0.0); bb=Float64(b); zz=Float64(zc)
    zff=Float64(zf); bd=Float64(bankfull_depth)
    yy <= bd && return bb + 2.0 * yy * zz
    return bb + 2.0 * bd * zz + 2.0 * zff * (yy - bd)
end

"""Wetted perimeter [m] matching mizuRoute `Pwet`."""
function wetted_perimeter(y::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=1.0e6)
    yy=max(Float64(y),0.0); bb=Float64(b); zz=Float64(zc)
    zff=Float64(zf); bd=Float64(bankfull_depth)
    yy <= bd && return bb + 2.0 * yy * sqrt(1.0 + zz^2)
    return bb + 2.0 * bd * sqrt(1.0 + zz^2) + 2.0 * (yy-bd) * sqrt(1.0 + zff^2)
end

"""Flow area [m²] matching mizuRoute `flow_area`."""
function flow_area(y::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=1.0e6)
    yy=max(Float64(y),0.0); bb=Float64(b); zz=Float64(zc)
    zff=Float64(zf); bd=Float64(bankfull_depth)
    yy <= bd && return yy * (bb + zz * yy)
    ab = bd * (bb + zz * bd)
    bt = top_width(yy,bb,zz;zf=zff,bankfull_depth=bd)
    btb = top_width(bd,bb,zz;zf=zff,bankfull_depth=bd)
    return ab + (yy-bd) * (bt+btb) / 2.0
end

"""Invert cross-sectional area to water height [m]."""
function water_height(area::Real, b::Real, zc::Real; zf::Real=1000.0, bankfull_depth::Real=1.0e6)
    a=max(Float64(area),0.0); bb=Float64(b); zz=Float64(zc)
    zff=Float64(zf); bd=Float64(bankfull_depth)
    ab=flow_area(bd,bb,zz;zf=zff,bankfull_depth=bd)
    if a > ab
        btb=top_width(bd,bb,zz;zf=zff,bankfull_depth=bd)
        disc=btb^2 - 4.0*zff*(ab-a)
        return bd + (-btb + sqrt(max(0.0,disc))) / (2.0*zff)
    elseif zz == 0.0
        return a / bb
    else
        return (-bb + sqrt(bb^2 + 4.0*a*zz)) / (2.0*zz)
    end
end

hydraulic_radius(y::Real,b::Real,zc::Real;zf::Real=1000.0,bankfull_depth::Real=1.0e6) =
    flow_area(y,b,zc;zf=zf,bankfull_depth=bankfull_depth) /
    wetted_perimeter(y,b,zc;zf=zf,bankfull_depth=bankfull_depth)

"""Manning uniform flow exactly following mizuRoute `uniformFlow`."""
function manning_discharge(y::Real,b::Real,zc::Real,slope::Real,mann_n::Real;
    zf::Real=1000.0,bankfull_depth::Real=1.0e6)
    yy=Float64(y); yy<=0.0 && return 0.0
    bb=Float64(b); zz=Float64(zc); s=Float64(slope); n=Float64(mann_n)
    zff=Float64(zf); bd=Float64(bankfull_depth)
    if yy <= bd
        a=flow_area(yy,bb,zz;zf=zff,bankfull_depth=bd)
        p=wetted_perimeter(yy,bb,zz;zf=zff,bankfull_depth=bd)
        return a*(a/p)^_CONST23*sqrt(s)/n
    end
    ab=flow_area(bd,bb,zz;zf=zff,bankfull_depth=bd)
    pb=wetted_perimeter(bd,bb,zz;zf=zff,bankfull_depth=bd)
    btb=top_width(bd,bb,zz;zf=zff,bankfull_depth=bd)
    ye=yy-bd
    ach=ab+btb*ye
    qch=ach*(ach/pb)^_CONST23*sqrt(s)/n
    afp=ye*zff*ye/2.0
    pfp=ye*sqrt(1.0+zff^2)
    qfp=afp>0 ? 2.0*(afp*(afp/pfp)^_CONST23*sqrt(s)/n) : 0.0
    return qch+qfp
end

"""Normal depth with the Newton iteration and tolerance used by mizuRoute v3.1."""
function flow_depth(q::Real,b::Real,zc::Real,slope::Real,mann_n::Real;
    zf::Real=1000.0,bankfull_depth::Real=1.0e6,rtol::Real=_FLOW_DEPTH_RTOL,maxiter::Integer=100)
    Q=Float64(q); Q<=_HYD_EPS && return 0.0
    bb=Float64(b); zz=Float64(zc); S=Float64(slope); n=Float64(mann_n)
    S<=0.0 && return 0.0
    zff=Float64(zf); bd=Float64(bankfull_depth)
    ab=flow_area(bd,bb,zz;zf=zff,bankfull_depth=bd)
    pb=wetted_perimeter(bd,bb,zz;zf=zff,bankfull_depth=bd)
    btb=top_width(bd,bb,zz;zf=zff,bankfull_depth=bd)
    qbf=ab*(ab/pb)^_CONST23*sqrt(S)/n
    err=100.0
    y=0.0
    if Q < qbf
        coef1=(sqrt(S)/n/Q)^3
        coef2=2.0*sqrt(zz^2+1.0)
        y0=(1.0/(coef1*bb^3))^(1.0/5.0)
        for _ in 1:maxiter
            a=flow_area(y0,bb,zz;zf=zff,bankfull_depth=bd)
            bt=top_width(y0,bb,zz;zf=zff,bankfull_depth=bd)
            p=wetted_perimeter(y0,bb,zz;zf=zff,bankfull_depth=bd)
            h=coef1*a^5/p^2-1.0
            dh=coef1*(5.0*a^4*bt*p-2.0*coef2*a^5)/p^3
            y=y0-h/dh
            if !(isfinite(y) && y>0.0); return max(y0,0.0); end
            err=abs((y-y0)/y)
            y0=y
            err <= rtol && break
        end
    else
        y0=bd+2.0
        coef1=sqrt(S)/n/pb^_CONST23
        coef2=2.0*(zff/2.0)^_CONST53*sqrt(S)/n/(zff^2+1.0)^(1.0/3.0)
        for _ in 1:maxiter
            ye=y0-bd
            h=coef1*(ab+btb*ye)^_CONST53 + coef2*ye^(_CONST103-_CONST23) - Q
            dh=coef1*_CONST53*btb*(ab+btb*ye)^_CONST23 +
               coef2*(_CONST103-_CONST23)*ye^_CONST53
            y=y0-h/dh
            if !(isfinite(y) && y>bd); return max(y0,bd); end
            err=abs((y-y0)/y)
            y0=y
            err <= rtol && break
        end
    end
    return y
end

storage(y::Real,length::Real,b::Real,zc::Real;zf::Real=1000.0,bankfull_depth::Real=1.0e6) =
    flow_area(y,b,zc;zf=zf,bankfull_depth=bankfull_depth)*Float64(length)

"""Wave celerity [m/s] using the exact mizuRoute hydraulic formula."""
function celerity(q::Real,b::Real,zc::Real,slope::Real,mann_n::Real;
    zf::Real=1000.0,bankfull_depth::Real=1.0e6)
    Q=abs(Float64(q)); Q<=_HYD_EPS && return 0.0
    y=flow_depth(Q,b,zc,slope,mann_n;zf=zf,bankfull_depth=bankfull_depth)
    y<=0.0 && return 0.0
    a=flow_area(y,b,zc;zf=zf,bankfull_depth=bankfull_depth)
    p=wetted_perimeter(y,b,zc;zf=zf,bankfull_depth=bankfull_depth)
    bt=top_width(y,b,zc;zf=zf,bankfull_depth=bankfull_depth)
    sf=(Q*Float64(mann_n)/a/(a/p)^_CONST23)^2
    return _CONST53*sf^0.3*Q^0.4/bt^0.4/Float64(mann_n)^0.6
end

"""Diffusive-wave diffusivity [m²/s] matching mizuRoute v3.1."""
function diffusivity(q::Real,b::Real,zc::Real,slope::Real,mann_n::Real;
    zf::Real=1000.0,bankfull_depth::Real=1.0e6)
    Q=abs(Float64(q)); Q<=_HYD_EPS && return 0.0
    y=flow_depth(Q,b,zc,slope,mann_n;zf=zf,bankfull_depth=bankfull_depth)
    y<=0.0 && return 0.0
    a=flow_area(y,b,zc;zf=zf,bankfull_depth=bankfull_depth)
    p=wetted_perimeter(y,b,zc;zf=zf,bankfull_depth=bankfull_depth)
    bt=top_width(y,b,zc;zf=zf,bankfull_depth=bankfull_depth)
    sf=(Q*Float64(mann_n)/a/(a/p)^_CONST23)^2
    return Q/sf/bt/2.0
end

@inline _reach_celerity(p::ReachParameters,i::Int,q::Real) =
    celerity(q,p.bottom_width[i],p.side_slope[i],p.slope[i],p.mann_n[i];
        zf=p.floodplain_slope[i],bankfull_depth=p.bankfull_depth[i])
@inline _reach_diffusivity(p::ReachParameters,i::Int,q::Real) =
    diffusivity(q,p.bottom_width[i],p.side_slope[i],p.slope[i],p.mann_n[i];
        zf=p.floodplain_slope[i],bankfull_depth=p.bankfull_depth[i])