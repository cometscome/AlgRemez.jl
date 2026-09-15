# Contributing

Bug reports, numerical regression cases, documentation improvements, and pull
requests are welcome.

## Development setup

Clone the repository, start Julia in its root directory, and instantiate the
project:

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Run the complete test suite before submitting a change:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests use `AlgRemez_jll` only as a test-time compatibility oracle.

## Numerical changes

Changes to the exchange, root-finding, or partial-fraction calculations should
include tests covering:

- Float64 ULP distance from `AlgRemez_jll` for `algorithm=:legacy`;
- extrema and dense-grid errors for `algorithm=:refined`;
- positive shifts and forward/inverse round trips;
- the degree-10 and degree-15 LatticeQCD RHMC profiles;
- generation time and allocations using `benchmark/rhmc.jl`.

Do not loosen numerical tolerances solely to make a regression pass. Explain
the numerical reason and application-level effect of any intentional change.

## Scope and compatibility

The public API follows semantic versioning. Avoid changing the fields or
meaning of `AlgRemez_coeffs`, because LatticeQCD code relies on that interface.
The legacy algorithm should remain behaviorally compatible with the original
C++ implementation; algorithmic experiments belong behind a separate option.
