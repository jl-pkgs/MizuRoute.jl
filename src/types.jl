abstract type AbstractRoutingMethod end

"""No in-channel delay; simply accumulates upstream and local lateral inflow."""
struct Accumulation <: AbstractRoutingMethod end

"""Reach-specific impulse response function routing.

`legacy_volume_limiter=true` reproduces the older serial mizuRoute IRF volume
limiter used by the bundled Cameo `ForComparison` output. The default `false`
matches the current mizuRoute main branch. The other keyword fields are retained
for API compatibility; the kernel itself uses mizuRoute's fixed hourly `make_uh`.
"""
struct IRF <: AbstractRoutingMethod
    horizon_factor::Float64
    min_steps::Int
    legacy_volume_limiter::Bool
end
IRF(; horizon_factor::Real=6.0, min_steps::Integer=8,
    legacy_volume_limiter::Bool=false) =
    IRF(Float64(horizon_factor), Int(min_steps), legacy_volume_limiter)

"""mizuRoute/TopNet Lagrangian kinematic-wave tracking.

`max_packets` corresponds to mizuRoute's `MAXQPAR` cap on the wave array.
`merge_rtol` is retained for backward API compatibility but is not used by the
exact Goring/TopNet shock-merging algorithm.
"""
struct LagrangianKWT <: AbstractRoutingMethod
    max_packets::Int
    merge_rtol::Float64
end
LagrangianKWT(; max_packets::Integer=20, merge_rtol::Real=0.0) =
    LagrangianKWT(Int(max_packets), Float64(merge_rtol))

"""Eulerian kinematic-wave routing using mizuRoute's shared implicit ADE solver."""
struct EulerKinematicWave <: AbstractRoutingMethod
    cells_per_reach::Int
    cfl::Float64
end
EulerKinematicWave(; cells_per_reach::Integer=20, cfl::Real=1.0) =
    EulerKinematicWave(Int(cells_per_reach), Float64(cfl))

"""Muskingum-Cunge routing with dynamic Cunge X and Courant-based substepping."""
struct MuskingumCunge <: AbstractRoutingMethod
    cfl::Float64
end
MuskingumCunge(; cfl::Real=0.9) = MuskingumCunge(Float64(cfl))

"""Diffusive-wave routing using the weighted centered finite-difference form
reported in the mizuRoute technical note.

`alpha` weights the new-time advection term and `beta` the new-time diffusion
term. Values near one are robust for large routing time steps.
"""
struct DiffusiveWave <: AbstractRoutingMethod
    cells_per_reach::Int
    alpha::Float64
    beta::Float64
end
DiffusiveWave(; cells_per_reach::Integer=20, alpha::Real=1.0, beta::Real=1.0) =
    DiffusiveWave(Int(cells_per_reach), Float64(alpha), Float64(beta))

"""Hydraulic and IRF parameters for all river reaches.

All vectors have one element per river reach. Scalar keyword values are expanded
across the network by the convenience constructor `ReachParameters(n; ...)`.

Units:
- `length`: m
- `slope`: m/m
- `bottom_width`: m
- `side_slope`: horizontal/vertical
- `floodplain_slope`: horizontal/vertical
- `bankfull_depth`: m
- `mann_n`: s m^(-1/3)
- `irf_velocity`: m/s
- `irf_diffusivity`: m^2/s
"""
struct ReachParameters{T<:AbstractFloat}
    length::Vector{T}
    slope::Vector{T}
    bottom_width::Vector{T}
    side_slope::Vector{T}
    floodplain_slope::Vector{T}
    bankfull_depth::Vector{T}
    mann_n::Vector{T}
    irf_velocity::Vector{T}
    irf_diffusivity::Vector{T}
end

_expand_parameter(x::Number, n::Integer) = fill(Float64(x), n)
function _expand_parameter(x::AbstractVector, n::Integer)
    length(x) == n || throw(ArgumentError("parameter length $(length(x)) != number of reaches $n"))
    Float64.(x)
end

function ReachParameters(n::Integer;
    length,
    slope,
    bottom_width,
    side_slope=0.0,
    floodplain_slope=1000.0,
    bankfull_depth=1.0e6,
    mann_n=0.035,
    irf_velocity=1.0,
    irf_diffusivity=1000.0,
)
    n > 0 || throw(ArgumentError("number of reaches must be positive"))
    p = ReachParameters(
        _expand_parameter(length, n),
        _expand_parameter(slope, n),
        _expand_parameter(bottom_width, n),
        _expand_parameter(side_slope, n),
        _expand_parameter(floodplain_slope, n),
        _expand_parameter(bankfull_depth, n),
        _expand_parameter(mann_n, n),
        _expand_parameter(irf_velocity, n),
        _expand_parameter(irf_diffusivity, n),
    )
    _validate_parameters(p)
    return p
end

function _validate_parameters(p::ReachParameters)
    n = length(p.length)
    fields = (
        p.slope, p.bottom_width, p.side_slope, p.floodplain_slope,
        p.bankfull_depth, p.mann_n, p.irf_velocity, p.irf_diffusivity,
    )
    all(length(x) == n for x in fields) || throw(ArgumentError("all reach parameter vectors must have equal length"))
    all(>(0), p.length) || throw(ArgumentError("reach length must be positive"))
    all(>=(0), p.slope) || throw(ArgumentError("reach slope must be non-negative"))
    all(>(0), p.bottom_width) || throw(ArgumentError("bottom width must be positive"))
    all(>=(0), p.side_slope) || throw(ArgumentError("side slope must be non-negative"))
    all(>(0), p.floodplain_slope) || throw(ArgumentError("floodplain slope must be positive"))
    all(>(0), p.bankfull_depth) || throw(ArgumentError("bankfull depth must be positive"))
    all(>(0), p.mann_n) || throw(ArgumentError("Manning n must be positive"))
    all(>(0), p.irf_velocity) || throw(ArgumentError("IRF velocity must be positive"))
    all(>=(0), p.irf_diffusivity) || throw(ArgumentError("IRF diffusivity must be non-negative"))
    return p
end

mutable struct BasicState{T<:AbstractFloat}
    qout::Vector{T}
    qin::Vector{T}
    volume::Vector{T}
    wm_actual::Vector{T}
end

mutable struct IRFState{T<:AbstractFloat}
    qout::Vector{T}
    qin::Vector{T}
    volume::Vector{T}
    wm_actual::Vector{T}
    history::Vector{Vector{T}}
    kernel::Vector{Vector{T}}
end

"""One mizuRoute KWT `FPOINT` wave entry.

`q` is unit-width discharge [m²/s], and `tentry`/`texit` are seconds on the
model time axis. `routed` corresponds to Fortran `RF`; `qmod` corresponds to
`QM` and is reserved for water-management compatibility.
"""
mutable struct KWTPoint{T<:AbstractFloat}
    q::T
    tentry::T
    texit::T
    routed::Bool
    qmod::T
end

mutable struct KWTState{T<:AbstractFloat}
    qout::Vector{T}
    qin::Vector{T}
    volume::Vector{T}
    wm_actual::Vector{T}
    waves::Vector{Vector{KWTPoint{T}}}
    qlat_prev::Vector{T}
end

mutable struct EulerKWState{T<:AbstractFloat}
    qout::Vector{T}
    qin::Vector{T}
    volume::Vector{T}
    wm_actual::Vector{T}
    profile::Vector{Vector{T}}
end

mutable struct MCState{T<:AbstractFloat}
    qout::Vector{T}
    qin::Vector{T}
    volume::Vector{T}
    wm_actual::Vector{T}
    input_prev::Vector{T}
    output_prev::Vector{T}
end

mutable struct DWState{T<:AbstractFloat}
    qout::Vector{T}
    qin::Vector{T}
    volume::Vector{T}
    wm_actual::Vector{T}
    profile::Vector{Vector{T}}
end

mutable struct RoutingModel{M<:AbstractRoutingMethod,S,T<:AbstractFloat}
    method::M
    network
    params::ReachParameters{T}
    dt::T
    time::T
    state::S
    active::BitVector
    headwater_drain_point::Int
end
