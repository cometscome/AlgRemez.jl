# AlgRemez.jl

`AlgRemez.jl` computes arbitrary-precision rational minimax approximations for
positive real functions of the form

```text
x^(y/z) * exp(sum(a[j] * x^k[j]))
```

It is a native Julia port of the AlgRemez implementation by M. A. Clark and
A. D. Kennedy. The default compatibility path reproduces the original Remez
exchange, root-finding, and partial-fraction algorithms. A separately refined
algorithm is also available. The public package contains no C or C++ source and
does not compile native code; `AlgRemez_jll` is used only by the test suite as a
reference executable.

The name follows the historical AlgRemez program and `AlgRemez_jll`. This
package is distinct from the registered `Remez.jl`: its compatibility API,
rational partial fractions, and tests are designed for the RHMC workflow.

The package is intended in particular for Rational Hybrid Monte Carlo (RHMC),
where an approximation is applied as

```text
α0 * b + sum(α[i] * (A + β[i]I)^(-1) * b)
```

All generated RHMC shifts are tested for positivity on the supported production
profiles.

## Installation

After registration in the Julia General registry:

```julia
import Pkg
Pkg.add("AlgRemez")
```

Before registration, install a checked-out copy or its Git URL:

```julia
import Pkg
Pkg.develop(path="/path/to/AlgRemez.jl")
```

AlgRemez requires Julia 1.10 or later.

## RHMC-compatible interface

`calc_coefficients` is a drop-in replacement for the historical
`AlgRemez_jll` wrapper:

```julia
using AlgRemez

positive, negative = calc_coefficients(
    1, 16, 15,
    4.0e-4, 64.0;
    precision=42,
)

positive.α0
positive.α
positive.β
positive.n

f = fittedfunction(positive)
f(0.1)
```

The returned approximation has the form

```text
α0 + sum(α[i] / (x + β[i]), i=1:n)
```

`positive` approximates `x^(y/z)` and `negative` approximates its reciprocal.
Their fields are `Float64`, matching the coefficient types expected by
LatticeQCD and LatticeDiracOperators.

The compatibility API interprets `precision` as decimal digits, as the C++ API
did. GMP rounds the requested storage to a machine-limb boundary, so
`precision=42` corresponds to 192 effective bits on a 64-bit host.

## Arbitrary-precision interface

Use `remez` when the polynomial coefficients and error information should
remain in `BigFloat`:

```julia
approximation = remez(
    1, 2, 5, 5,
    4.0e-4, 64.0;
    precision=192,       # bits
    tolerance=1.0e-15,
)

r = approximation(big"0.1")
f = target_value(approximation, big"0.1")
δ = relative_error(approximation, big"0.1")

positive_pfe = partial_fraction(approximation)
negative_pfe = partial_fraction(approximation; inverse=true)
```

The main arguments are:

- `y`, `z`: the target exponent `y/z`;
- `n`, `d`: numerator and denominator degrees;
- `lower`, `upper`: a strictly positive approximation interval;
- `precision`: `BigFloat` working precision in bits;
- `tolerance`: tolerated spread among error extrema.

Equal and unequal numerator/denominator degrees are supported for rational-form
evaluation. Partial-fraction conversion requires `n == d` and distinct real
poles. The safe constant case `(n, d) == (0, 0)` is supported.

### Algorithm selection

The default is the original calculation:

```julia
approximation = remez(args...; algorithm=:legacy)
pfe = partial_fraction(approximation; algorithm=:legacy)
```

This path uses the original initial points, scaled-pivot Gaussian elimination,
bounded ten-step extrema search, relaxation update, Newton deflation, and
residue construction. It rounds interval bounds through `Float64`, matching the
C++ interface.

The refined implementation is selected explicitly:

```julia
approximation = remez(args...; algorithm=:refined, precision=256)
pfe = partial_fraction(approximation; algorithm=:refined)
```

It uses logarithmic initial points, golden-section extrema searches, and
derivative-root isolation. String and `BigFloat` interval endpoints retain
their full precision in this mode.

### Power times exponential targets

The extended target from the original class is available with parallel arrays:

```julia
approximation = remez(
    1, 2, 4, 3,
    0.01, 10.0;
    exponential_coefficients=[0.1, -0.02],
    exponential_powers=[1, -1],
)
```

This represents
`sqrt(x) * exp(0.1x - 0.02/x)`.

## Accuracy and compatibility

The test suite compares every Float64 normalization, residue, and shift from
the RHMC compatibility API with `AlgRemez_jll`. The cross-platform allowance is
four ULPs; the current development host differs by at most one ULP. The original
and Julia implementations also take the same number of Remez iterations in the
reference cases.

The often quoted maximum relative error `2.550823e-3` belongs to the deliberately
low degree `(5, 5)` square-root example on the wide interval `[0.0004, 64]`. It
is a regression value, not a recommended RHMC accuracy. Choose the degree from
the required approximation error and spectral interval.

The LatticeQCD-oriented tests cover `Nf = 1, 2, 3` for:

- degree-15 action approximations `x^(Nf/16)`;
- degree-10 molecular-dynamics approximations `x^(Nf/8)`;
- positive shifted systems;
- forward/inverse application and round trips;
- negative-power coefficient dispatch.

## Command-line compatibility

The original seven-argument interface and output files are available through:

```sh
julia --project=. bin/algremez.jl 1 2 5 5 0.0004 64 42
```

This writes `approx.dat` and `error.dat` in the historical format.

## Testing and benchmarking

Run the complete test suite with:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run JIT-warmed coefficient-generation benchmarks with:

```sh
julia --project=. benchmark/rhmc.jl
julia --project=. benchmark/rhmc.jl --all
```

On the development host, representative Julia/C++ generation times were
0.47/0.65 seconds for degree 5, 2.70/2.97 seconds for degree 10, and 8.55/8.89
seconds for degree 15. Timings vary by host and coefficient generation is a
one-time setup cost rather than a per-trajectory RHMC cost.

More implementation details are in [README_JULIA.md](README_JULIA.md), and the
original C++ code analysis is available as
[AlgRemez_analysis.pdf](docs/AlgRemez_analysis.pdf).
Maintainer release steps are documented in
[docs/RELEASING.md](docs/RELEASING.md).

## Citation and license

If this software contributes to published work, cite the original AlgRemez
implementation:

> M. A. Clark and A. D. Kennedy, AlgRemez (2005),
> <https://github.com/mikeaclark/AlgRemez>.

Please also cite the version of `AlgRemez.jl` used in the calculation. The
software is distributed under the MIT license; see [LICENSE](LICENSE).
