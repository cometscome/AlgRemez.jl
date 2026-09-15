# AlgRemez.jl

The native Julia implementation covers the calculations exposed by the original
C++ `AlgRemez` class: unequal or equal numerator/denominator degrees, the
power-times-exponential target, direct and inverse evaluation, and positive and
inverse partial-fraction expansions. All core arithmetic, coefficients, error
zeros, and extrema use `BigFloat`.

```julia
using AlgRemez

approximation = remez(
    1, 2, 5, 5,
    "0.0004", "64";
    precision=192,
    tolerance=1.0e-15,
)

value = approximation(big"0.1")
error = relative_error(approximation, big"0.1")

positive = partial_fraction(approximation)
negative = partial_fraction(approximation; inverse=true)

inverse_approximation = inverse_value(approximation, big"0.1")
exact_target = target_value(approximation, big"0.1")
inverse_target = inverse_target_value(approximation, big"0.1")
```

The positional arguments are exponent numerator `y`, exponent denominator `z`,
numerator degree `n`, denominator degree `d`, and the lower and upper interval
bounds. `precision` is measured in bits. The default `algorithm=:legacy`
faithfully follows the original initial points, scaled-pivot Gaussian
elimination, ten-step extrema search, zero update, Newton deflation, and PFE
construction. Like the C++ API, this mode rounds interval bounds through
`Float64`. `y` may be zero (for a constant power factor), while `z` must be
positive.

An independently refined variant remains available with `algorithm=:refined`.
It uses logarithmic initial points, golden-section extrema searches, and
derivative-root isolation. In that mode, string or `BigFloat` endpoints avoid
rounding through `Float64`.

The extended target supported by the C++ overload

```text
x^(y/z) * exp(sum(a[j] * x^k[j]))
```

is selected with two parallel keyword arrays. Integer powers may be positive,
zero, or negative because the approximation interval is strictly positive.

```julia
approximation = remez(
    1, 2, 4, 3, "0.01", "10";
    exponential_coefficients=["0.1", "-0.02"],
    exponential_powers=[1, -1],
)
```

`RationalApproximation` contains ascending-order numerator and denominator
coefficients together with error zeros, extrema, signed extrema errors,
iteration count, and maximum relative error. Its monic denominator includes the
leading coefficient as its last entry. The legacy extrema use the original
bounded marching search; independently refined extrema are provided only by
`algorithm=:refined`.

`partial_fraction` supports equal numerator and denominator degrees with
distinct real poles. It returns

```text
constant + sum(residues[i] / (x + shifts[i]))
```

Unlike the legacy C++ executable, `(n, d) = (0, 0)` is supported safely. The
legacy `AlgRemez_jll` executable is used only by the test suite as a
compatibility oracle. Pass `algorithm=:refined` to `partial_fraction` to select
derivative-root isolation instead of the default legacy Newton-deflation path.

## Compatibility with the existing Julia wrapper

Code in LatticeDiracOperators.jl can continue to use the previous interface:

```julia
coeff_plus, coeff_minus = calc_coefficients(
    1, 2, 10, 4.0e-4, 64.0;
    precision=42,
)
f = fittedfunction(coeff_plus)
```

`AlgRemez_coeffs` retains the fields `α0`, `α`, `β`, and `n`, and deliberately
contains `Float64` values to remain source-compatible. Use `remez` followed by
`partial_fraction` when the coefficients must remain at arbitrary precision.
For compatibility, `calc_coefficients` measures `precision` in decimal digits;
`remez` measures it in bits. `calc_coefficients` reproduces GMP's precision
rounding: the original request for 42 decimal digits becomes 139 requested bits
and then 192 effective bits on a 64-bit host. It also uses the same legacy
exchange and PFE algorithms as the C++ implementation.

## Legacy files and command line

`write_legacy_files(approximation)` writes the original `approx.dat` and
`error.dat` formats. The command-line entry point accepts the same seven
arguments as the C++ executable:

```sh
julia --project=. bin/algremez.jl 1 2 5 5 0.0004 64 42
```

The file-producing interface requires equal degrees because both the original
format and the original `getPFE`/`getIPFE` routines do. Unequal-degree rational
approximations are fully available through `remez` and direct evaluation.

The legacy class operations map to the Julia API as follows:

| Original C++ operation | Native Julia operation |
| --- | --- |
| `generateApprox(n,d,pnum,pden)` | `remez(pnum,pden,n,d,lower,upper)` |
| exponential `generateApprox` overload | `remez(...; exponential_coefficients, exponential_powers)` |
| `evaluateApprox(x)` | `approximation(x)` |
| `evaluateInverseApprox(x)` | `inverse_value(approximation,x)` |
| `evaluateFunc(x)` | `target_value(approximation,x)` |
| `evaluateInverseFunc(x)` | `inverse_target_value(approximation,x)` |
| `getPFE(...)` | `partial_fraction(approximation)` |
| `getIPFE(...)` | `partial_fraction(approximation; inverse=true)` |
| `setBounds(lower,upper)` | pass the new bounds to the next `remez` call |

## LatticeQCD RHMC verification

The test suite reproduces the coefficient and shifted-solve contract used by
`LatticeDiracOperators.Rhmc.RHMC`. For `Nf = 1, 2, 3`, it checks both the
degree-15 action approximation `x^(Nf/16)` and the degree-10 molecular-dynamics
approximation `x^(Nf/8)` on `[0.0004, 64]`. These tests require:

- all PFE shifts to be positive, so every shifted `D†D + β[i]` system remains
  positive definite;
- action-profile relative errors below `5e-9`;
- molecular-dynamics-profile relative errors below `5e-6`;
- forward/inverse PFE round-trip error below `5e-12`;
- correct coefficient swapping for negative rational powers.

The test applies each PFE to a positive-definite diagonal matrix using the same
`α0*b + sum(α[i]*(A + β[i]*I)\b)` algebra implemented with `shiftedcg` in
LatticeQCD. Exact shifted solves are used so failures isolate coefficient
generation rather than iterative-solver tolerance.

All Float64 values written by the compatibility API (`α0`, residues, and
shifts) are compared directly with `AlgRemez_jll` using ULP distance. The
allowed cross-platform tolerance is four ULPs; the current environment produces
at most one ULP difference for the README case and all six RHMC profiles. The
standard case also requires exactly the same Remez iteration count as the JLL.

## Performance

Run the reproducible, JIT-warmed benchmark with:

```sh
julia --project=. benchmark/rhmc.jl
```

Pass `--all` to calculate all six RHMC profiles. A JIT-warmed run on the current
test machine showed comparable or slightly faster coefficient generation than
the legacy C++ executable:

| Profile | Julia | Legacy C++ | Julia/C++ | Julia allocation |
| --- | ---: | ---: | ---: | ---: |
| degree-5 `x^(1/2)` | 0.47 s | 0.65 s | 0.72x | 0.19 GiB |
| RHMC MD, degree-10 `x^(1/4)` | 2.70 s | 2.97 s | 0.91x | 1.63 GiB |
| RHMC action, degree-15 `x^(1/16)` | 8.55 s | 8.89 s | 0.96x | 6.34 GiB |

These timings exclude Julia compilation and vary by host. RHMC evolution uses
the resulting Float64 PFE coefficients; these generation times are startup
costs rather than per-trajectory costs.

Run the tests with:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```
