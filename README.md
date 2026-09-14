# MizuRoute.jl

[![CI](https://github.com/jl-pkgs/MizuRoute.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/jl-pkgs/MizuRoute.jl/actions/workflows/ci.yml)

A lightweight Julia river-routing kernel that reimplements the principal
algorithms documented in ESCOMP/mizuRoute for direct coupling to hydrological
and land-surface models.

## Included core algorithms

- runoff accumulation (diagnostic baseline)
- gamma-distribution hillslope unit hydrograph
- reach IRF routing using the diffusive-wave Green function
- Lagrangian kinematic-wave characteristic routing (KWT-style)
- Eulerian kinematic wave through the shared implicit ADE solver
- Muskingum-Cunge with dynamic Cunge weighting and CFL substepping
- weighted implicit finite-difference diffusive wave through the shared ADE solver
- river-network construction and topological traversal
- compound trapezoidal/floodplain hydraulics
- conservative grid/HRU-to-reach runoff mapping
- water management sign convention compatible with mizuRoute
- stepwise water-balance accounting

This package deliberately does **not** port mizuRoute's MPI, PIO, CESM/CTSM
coupling, NetCDF control-file parser, restart machinery, logging framework, or
other infrastructure. The objective is a small routing kernel that can be
embedded in another Julia hydrological model.

## Quick start

```julia
using MizuRoute

net = RiverNetwork(
    [101, 102, 103],       # reach IDs
    [103, 103, 0],         # downstream reach IDs; 0 = outlet
)

p = ReachParameters(3;
    length=[4000.0, 3500.0, 8000.0],
    slope=[0.002, 0.003, 0.001],
    bottom_width=[8.0, 6.0, 15.0],
    side_slope=1.0,
    bankfull_depth=2.5,
    mann_n=0.035,
    irf_velocity=1.2,
    irf_diffusivity=300.0,
)

river = RoutingModel(MuskingumCunge(), net, p; dt=3600.0)

qlat = [2.0, 1.0, 0.5]  # m3/s local lateral inflow
step!(river, qlat)
Q = discharge(river)
```

For a gridded hydrological model whose runoff is in m/s:

```julia
qlat = map_runoff(
    length(net),
    runoff_depth,
    cell_area,
    cell_to_reach,
)
step!(river, qlat)
```

## Choosing a routing method

```julia
IRF()
LagrangianKWT()
EulerKinematicWave(cells_per_reach=12)
MuskingumCunge()
DiffusiveWave(cells_per_reach=12, alpha=1.0, beta=1.0)
```

## Tests

From the package directory:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The supplied test suite covers topology, hydraulic geometry, Manning inversion,
gamma-UH normalization, remapping, the ADE solver, non-negativity, finite values,
and stepwise water conservation for all routing methods. See `TEST_STATUS.md`
for the execution status of this artifact build.

## Documentation

The technical documentation is written in Typst:

```bash
cd docs
typst compile main.typ MizuRoute.pdf
```

Typst >= 0.15.1 is recommended. The document intentionally uses no external
Typst Universe packages, so it can compile offline.

## Fidelity note

The equations and algorithmic structure follow the current mizuRoute technical
note and source organization, but this is an independent Julia architecture,
not a line-by-line Fortran translation. In particular, the compact Lagrangian
KWT packet bookkeeping differs from the historical TopNet/mizuRoute internal
wave data structure. Bit-for-bit regression against the official mizuRoute
Cameo test case is therefore a separate validation milestone. See `docs/validation.typ`, `PORTING_MAP.md`, and `TEST_STATUS.md`.

## References

- Mizukami et al. (2016), GMD, doi:10.5194/gmd-9-2223-2016
- Mizukami et al. (2021), JAMES, doi:10.1029/2020MS002434
- Cortés-Salazar et al. (2023), HESS, doi:10.5194/hess-27-3505-2023
- ESCOMP/mizuRoute documentation: https://mizuroute.readthedocs.io/
