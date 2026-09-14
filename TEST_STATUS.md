# Test and validation status

## Included Julia tests

Run from the package root:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite contains tests for:

- river-network topological ordering and cycle rejection;
- compound cross-section geometry and area/depth inversion;
- Manning discharge/depth round trip;
- gamma-UH normalization;
- runoff grid-to-reach remapping;
- shared advection-diffusion tridiagonal solver invariants;
- all routing methods on a multi-reach network;
- non-negativity and finite-value checks;
- reach-by-reach, time-step water-balance closure;
- `route_series` and `reset!` API behavior.

## GitHub Actions execution status

The repository CI has been executed successfully on GitHub Actions. The Julia
package loads and `Pkg.test()` passes on Julia **1.10, 1.11, and 1.12**. The
Typst technical documentation also compiles successfully with **Typst 0.15.1**
and is uploaded as the `MizuRoute-technical-documentation` workflow artifact.

Static source checks and an independent Python numerical mirror/sanity test are
also included in `validation/`.

## Remaining reference validation

A stronger validation milestone is a numerical regression against the official
mizuRoute Cameo/testCase using identical network, forcing, routing interval and
parameters. That is intentionally distinguished from internal unit testing,
especially for the historical KWT particle bookkeeping.
