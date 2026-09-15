# mizuRoute → MizuRoute.jl porting map

This project intentionally ports the **routing physics and numerical kernels**,
not the CESM/MPI/PIO application framework.

| mizuRoute source | MizuRoute.jl | Porting status |
|---|---|---|
| `network_topo.f90`, `process_ntopo.f90` | `src/network.jl` | core upstream/downstream graph and topological processing retained; metadata machinery omitted |
| `hydraulic.f90` | `src/hydraulics.jl` | compound-section geometry, Manning flow, depth inversion, celerity and diffusivity |
| `basinUH.f90`, `process_param.f90` | `src/hillslope.jl` | gamma unit-hydrograph core |
| `irf_route.f90` | `src/routing/irf.jl` | diffusive-wave Green-function IRF convolution |
| `kwt_route.f90` | `src/routing/kwt.jl` | Lagrangian particle/characteristic algorithm; Julia state bookkeeping is compact rather than Fortran-identical |
| `kwe_route.f90` | `src/routing/euler_kw.jl` | linearized Euler kinematic wave |
| `advection_diffusion.f90` | `src/routing/advection_diffusion.jl` | centered/upwind weighted ADE + Thomas solver, shared by KW/DW |
| `mc_route.f90` | `src/routing/muskingum_cunge.jl` | dynamic Cunge X, hydraulic celerity, substeps |
| `dfw_route.f90` | `src/routing/diffusive_wave.jl` | linearized diffusive wave through shared ADE solver |
| `process_remap.f90` | `src/remap.jl` | model-grid/HRU runoff to reach mapping in a lightweight array interface |
| `water_balance.f90` | `src/balance.jl` | reach-scale conservation check and non-negative storage update |
| water-management logic in route modules | `src/balance.jl`, `src/model.jl` | positive abstraction / negative injection convention and withdrawal priority |
| MPI / PIO / restart / control parser / CESM glue | intentionally omitted | application infrastructure, not routing physics |

## Deliberate numerical differences

1. The KWT implementation preserves the documented Lagrangian particle concept,
   including travel time based on kinematic celerity and within-step movement,
   but does not clone the historical `FPOINT/KWAVE` memory layout.
2. Current mizuRoute `solve_ade` accepts `FluxLat` and documents it in the PDE,
   while the current Fortran RHS does not explicitly add it. The Julia kernel
   includes the documented uniform lateral source by default; it can be disabled
   internally for bug-for-bug ADE experiments.
3. Package data structures and coupling API are idiomatic Julia rather than a
   transliteration of Fortran derived types.
