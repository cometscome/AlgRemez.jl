# Test strategy

`AlgRemez_jll` is a test-only compatibility oracle for the legacy C++
implementation. `reference.jl` runs it in a separate temporary working
directory because the executable always writes `approx.dat` and `error.dat`.

The compatibility checks compare every Float64 normalization, residue, and
shift directly with the JLL, allowing at most four ULPs across platforms. On
the current host the faithful Julia path differs by at most one ULP, and the
standard case has exactly the same Remez iteration count. The legacy
`error.dat` is still only a coarse diagnostic grid and must not be used as a
refined minimax certificate.

The approximately `2.550823e-3` relative error in `runtests.jl` is the golden
result for the deliberately low `(5, 5)` README example over the very wide
interval `[0.0004, 64]`. It is not an accuracy requirement. Native Julia
accuracy tests must specify a required error first and use a degree high enough
to meet it.

The Julia tests cover both the faithful legacy path and the independently
refined path:

1. the analytic `(0, 0)` constant approximation;
2. the exact linear case and interpolation residuals at the exchanged zeros;
3. a refined certificate with `n + d + 2` alternating, equal-magnitude error extrema;
4. an independent dense grid plus positive partial-fraction shifts, placing the
   standard case's poles outside its positive approximation interval;
5. evaluation-based differential checks against `STANDARD_REFERENCE`;
6. agreement between coefficient and partial-fraction forms, including the
   reciprocal approximation;
7. power-times-exponential targets, including negative integer powers;
8. direct and inverse target/approximation evaluation;
9. the `AlgRemez_coeffs`/`calc_coefficients` compatibility API;
10. the legacy initial-point, search, Gaussian-elimination, Newton-deflation,
    and PFE path, including the standard-case iteration count;
11. the legacy `approx.dat`, `error.dat`, and command-line entry point.

`rhmc.jl` additionally exercises all six production profiles used by
LatticeQCD for `Nf = 1, 2, 3`: degree-15 `x^(Nf/16)` action approximations and
degree-10 `x^(Nf/8)` molecular-dynamics approximations on `[0.0004, 64]`. It
applies the generated PFE to a positive-definite matrix through shifted linear
solves, checks positive and inverse powers, verifies positive shifts, tests the
negative-order coefficient swap, and checks their round trip. Every Float64
normalization, residue, and shift is also compared with `AlgRemez_jll` and must
be within four ULPs of the legacy output.

The current standard-case grid is independent of the extrema search and the
partial-fraction reconstruction is independent of coefficient-form evaluation.
Additional high-degree accuracy tiers should state an application error target
explicitly; the `(5, 5)` regression is not such a target.

Do not call the JLL for `(n, d) = (0, 0)`: its partial-fraction routine writes
through zero-length arrays and the executable crashes.
