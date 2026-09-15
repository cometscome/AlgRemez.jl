using AlgRemez
using Printf

const LOWER = 4.0e-4
const UPPER = 64.0

function benchmark_profile(label, power::Rational, degree::Int)
    y = numerator(power)
    z = denominator(power)
    GC.gc()
    result = @timed calc_coefficients(
        y,
        z,
        degree,
        LOWER,
        UPPER;
        precision=42,
    )
    positive, negative = result.value
    grid = exp.(range(log(LOWER), log(UPPER); length=4_001))
    positive_error = maximum(abs(positive(x) / x^(y / z) - 1) for x in grid)
    negative_error = maximum(abs(negative(x) * x^(y / z) - 1) for x in grid)
    @printf(
        "%-12s x^(%d/%d), n=%2d: %8.3f s, %6.3f GiB, errors=(%.6e, %.6e)\n",
        label,
        y,
        z,
        degree,
        result.time,
        result.bytes / 2.0^30,
        positive_error,
        negative_error,
    )
end

# Remove compilation from the reported profiles.
calc_coefficients(1, 2, 1, 1.0, 4.0; precision=42, tolerance="1e-6")

if "--all" in ARGS
    benchmark_profile("Nf=1 action", 1 // 16, 15)
    benchmark_profile("Nf=2 action", 1 // 8, 15)
    benchmark_profile("Nf=3 action", 3 // 16, 15)
    benchmark_profile("Nf=1 MD", 1 // 8, 10)
    benchmark_profile("Nf=2 MD", 1 // 4, 10)
    benchmark_profile("Nf=3 MD", 3 // 8, 10)
else
    benchmark_profile("README", 1 // 2, 5)
    benchmark_profile("RHMC MD", 1 // 4, 10)
    benchmark_profile("RHMC action", 1 // 16, 15)
end
