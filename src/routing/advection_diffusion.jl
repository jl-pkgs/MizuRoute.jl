"""Thomas algorithm for a tridiagonal system.

`lower[i]` is A[i,i-1], `diag[i]` is A[i,i], and `upper[i]` is A[i,i+1].
The unused entries `lower[1]` and `upper[end]` are ignored.
"""
function _tdma(lower::Vector{Float64}, diag::Vector{Float64},
    upper::Vector{Float64}, rhs::Vector{Float64})
    n = length(diag)
    (length(lower) == n && length(upper) == n && length(rhs) == n) ||
        throw(ArgumentError("tridiagonal length mismatch"))
    n >= 2 || throw(ArgumentError("tridiagonal system needs at least two rows"))

    d = copy(diag)
    r = copy(rhs)
    @inbounds for i in 2:n
        abs(d[i - 1]) > 1.0e-14 || throw(ArgumentError("singular tridiagonal system"))
        m = lower[i] / d[i - 1]
        d[i] -= m * upper[i - 1]
        r[i] -= m * r[i - 1]
    end

    abs(d[n]) > 1.0e-14 || throw(ArgumentError("singular tridiagonal system"))
    x = similar(rhs)
    x[n] = r[n] / d[n]
    @inbounds for i in (n - 1):-1:1
        abs(d[i]) > 1.0e-14 || throw(ArgumentError("singular tridiagonal system"))
        x[i] = (r[i] - upper[i] * x[i + 1]) / d[i]
    end
    return x
end

"""Solve the linearized 1-D advection-diffusion equation used by mizuRoute.

The discretization mirrors `route/build/src/advection_diffusion.f90`:

    dQ/dt + c dQ/dx = D d²Q/dx² + s_lat

with an upstream Dirichlet condition. The default downstream boundary is the
mizuRoute Neumann condition that preserves the previous-step outlet gradient.
`wc=wd=1` gives the default fully implicit scheme. `advection=:central` is the
mizuRoute default; `:upwind` is also provided.

Important for numerical regression: the current mizuRoute Fortran routine
accepts `FluxLat` but does not add it to the right-hand side, and it does not
clip negative values after the tridiagonal solve. Set `include_lateral=false`
to reproduce that behavior (the river-routing kernels do so).
"""
function _solve_ade(qprev::Vector{Float64}, reach_length::Float64, dt::Float64,
    qupstream::Float64, ck::Float64, dk::Float64, qlat::Float64=0.0;
    advection::Symbol=:central, downstream::Symbol=:neumann,
    wc::Float64=1.0, wd::Float64=1.0, include_lateral::Bool=true)

    n = length(qprev)
    n >= 4 || throw(ArgumentError("ADE solver requires at least 4 molecular nodes"))
    reach_length > 0 || throw(ArgumentError("reach_length must be positive"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
    0.0 <= wc <= 1.0 || throw(ArgumentError("wc must be in [0,1]"))
    0.0 <= wd <= 1.0 || throw(ArgumentError("wd must be in [0,1]"))
    advection in (:central, :upwind) || throw(ArgumentError("advection must be :central or :upwind"))
    downstream in (:neumann, :absorbing) || throw(ArgumentError("downstream must be :neumann or :absorbing"))

    # Match mizuRoute exactly: Nx=nMolecule-1 and dx=L/(Nx-1), i.e. one
    # computational sub-segment extends beyond the nominal outlet.
    nx = n - 1
    dx = reach_length / (nx - 1)
    cd = dk * dt / dx^2
    ca = ck * dt / dx

    lower = zeros(Float64, n)
    diag  = zeros(Float64, n)
    upper = zeros(Float64, n)
    rhs   = zeros(Float64, n)

    # Upstream Dirichlet boundary. Do not clamp: Fortran assigns FluxUpstream
    # directly, and exact regression must preserve the same algebra.
    diag[1] = 1.0
    rhs[1] = qupstream

    # `FluxLat` is presently unused by mizuRoute's Fortran ADE RHS. Keep an
    # opt-in source term for standalone use, but all parity kernels disable it.
    source = include_lateral ? 2.0 * dt * ck * qlat / reach_length : 0.0

    if advection === :upwind
        @inbounds for j in 2:(n - 1)
            lower[j] = -wc * ca - wd * cd
            diag[j]  = 1.0 + wc * ca + 2.0 * wd * cd
            upper[j] = -wd * cd
            rhs[j] = ((1.0 - wc) * ca + (1.0 - wd) * cd) * qprev[j - 1] +
                     (1.0 - (1.0 - wc) * ca - 2.0 * (1.0 - wd) * cd) * qprev[j] +
                     (1.0 - wd) * cd * qprev[j + 1] + 0.5 * source
        end
    else
        @inbounds for j in 2:(n - 1)
            lower[j] = -wc * ca - 2.0 * wd * cd
            diag[j]  = 2.0 + 4.0 * wd * cd
            upper[j] =  wc * ca - 2.0 * wd * cd
            rhs[j] = ((1.0 - wc) * ca + 2.0 * (1.0 - wd) * cd) * qprev[j - 1] +
                     (2.0 - 4.0 * (1.0 - wd) * cd) * qprev[j] -
                     ((1.0 - wc) * ca - 2.0 * (1.0 - wd) * cd) * qprev[j + 1] + source
        end
    end

    if downstream === :absorbing
        lower[n] = -wc * ca
        diag[n] = 1.0 + wc * ca
        rhs[n] = (1.0 - (1.0 - wc) * ca) * qprev[n] +
                 (1.0 - wc) * ca * qprev[n - 1]
    else
        # mizuRoute Neumann BC retains the previous gradient, rather than
        # forcing a zero gradient at every step.
        lower[n] = -1.0
        diag[n] = 1.0
        rhs[n] = qprev[n] - qprev[n - 1]
    end

    # Do not post-process the solution. The Fortran TDMA returns the raw
    # tridiagonal solution; clipping here changes KW/DW propagation.
    return _tdma(lower, diag, upper, rhs)
end
