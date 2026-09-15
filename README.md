# MizuRoute.jl

[![CI](https://github.com/jl-pkgs/MizuRoute.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/jl-pkgs/MizuRoute.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/jl-pkgs/MizuRoute.jl/graph/badge.svg)](https://codecov.io/gh/jl-pkgs/MizuRoute.jl)

A lightweight Julia river-routing kernel that reimplements the principal
algorithms documented in ESCOMP/mizuRoute for direct coupling to hydrological
and land-surface models.

This project is an independent Julia reimplementation of core river-routing
ideas and equations documented by ESCOMP/mizuRoute. mizuRoute is developed by
NCAR/ESCOMP and contributors and is distributed under the Apache License 2.0.

The implementation in this package does not copy the mizuRoute MPI/PIO/CESM
infrastructure and is not an official ESCOMP product. See docs/validation.typ
for the current fidelity and validation status.


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

## Inputs required by each method

All methods require a `RiverNetwork`, a routing interval `dt` [s], and lateral
inflow `qlat` [m³/s] for every reach and time step. The routing methods use the
following fields from `ReachParameters`:

| Method                 | Reach data used                                                                                 |
| ---------------------- | ----------------------------------------------------------------------------------------------- |
| `Accumulation()`       | None; it only sums upstream and lateral inflow                                                  |
| `IRF()`                | `length`, `irf_velocity`, `irf_diffusivity`                                                     |
| `LagrangianKWT()`      | `length`, `slope`, `mann_n`                                                                     |
| `EulerKinematicWave()` | `length`, `slope`, `bottom_width`, `side_slope`, `floodplain_slope`, `bankfull_depth`, `mann_n` |
| `MuskingumCunge()`     | `length`, `slope`, `bottom_width`, `side_slope`, `floodplain_slope`, `bankfull_depth`, `mann_n` |
| `DiffusiveWave()`      | `length`, `slope`, `bottom_width`, `side_slope`, `floodplain_slope`, `bankfull_depth`, `mann_n` |

`ReachParameters` currently requires `length`, `slope`, and `bottom_width` even
when the selected method does not use all three. Every field may be either one
scalar shared by all reaches or a vector of length `nreach`; vector order must
match `reach_id`. Gridded input additionally requires runoff depth [m/s], cell
area [m²], and an internal reach index for each cell, as shown above.

### Where flow direction enters

MizuRoute.jl does not read a D8/D∞ flow-direction raster or route water across
land grid cells. Flow direction is used during preprocessing to derive:

1. `cell_to_reach`: the receiving reach for each runoff cell or HRU;
2. `downstream_id`: the downstream reach for each river reach.

At runtime, `map_runoff` immediately aggregates `runoff_depth * cell_area` into
local reach inflow `qlat`; optional `GammaUHRouter` adds hillslope travel-time
delay, and the selected routing method then routes flow through `RiverNetwork`.
If the host hydrological model already routes runoff to the channel, pass its
channel-entry discharge directly as `qlat` and skip both `map_runoff` and
`GammaUHRouter`. `qlat` must exclude upstream-reach discharge, which MizuRoute
adds internally.

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
