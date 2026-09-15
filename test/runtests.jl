using AlgRemez
using LinearAlgebra
using Test

include("reference.jl")
using .AlgRemezReference

const STANDARD_INPUT = (
    y=1,
    z=2,
    n=5,
    d=5,
    lower=4.0e-4,
    upper=64.0,
    precision=42,
)

const STANDARD_REFERENCE = jll_reference(
    STANDARD_INPUT.y,
    STANDARD_INPUT.z,
    STANDARD_INPUT.n,
    STANDARD_INPUT.d,
    STANDARD_INPUT.lower,
    STANDARD_INPUT.upper;
    precision=STANDARD_INPUT.precision,
)

@testset "native analytic cases" begin
    original_precision = precision(BigFloat)
    constant = remez(
        1,
        2,
        0,
        0,
        "1",
        "4";
        algorithm=:refined,
        precision=192,
        tolerance="1e-20",
    )

    @test precision(BigFloat) == original_precision
    @test length(constant.numerator) == 1
    @test constant.denominator == BigFloat[1]
    @test isapprox(constant.numerator[1], big"4" / 3; atol=big"1e-20", rtol=0)
    @test isapprox(
        constant.maximum_relative_error,
        big"1" / 3;
        atol=big"1e-20",
        rtol=0,
    )
    @test constant.extrema == BigFloat[1, 4]
    @test signbit(constant.extrema_errors[1]) != signbit(constant.extrema_errors[2])
    @test isempty(constant.exponential_coefficients)
    @test isempty(constant.exponential_powers)
    @test isapprox(
        inverse_value(constant, big"2") * constant(big"2"),
        big"1";
        atol=big"1e-50",
    )
    @test isapprox(
        inverse_target_value(constant, big"2") * target_value(constant, big"2"),
        big"1";
        atol=big"1e-50",
    )

    constant_pfe = partial_fraction(constant)
    inverse_constant_pfe = partial_fraction(constant; inverse=true)
    @test isempty(constant_pfe.residues)
    @test isapprox(constant_pfe.constant, big"4" / 3; atol=big"1e-20", rtol=0)
    @test isapprox(inverse_constant_pfe.constant, big"3" / 4; atol=big"1e-20", rtol=0)

    legacy_constant = AlgRemez_coeffs(constant_pfe)
    @test legacy_constant isa AlgRemez_coeffs
    @test legacy_constant.n == 0
    @test fittedfunction(legacy_constant)(2.0) == legacy_constant(2.0)
    @test fittedfunction(constant_pfe)(big"2") == constant_pfe(big"2")

    coeff_plus, coeff_minus = calc_coefficients(
        1,
        2,
        0,
        "1",
        "4";
        precision=20,
        tolerance="1e-10",
    )
    @test isapprox(coeff_plus.α0, 4 / 3; rtol=1e-9)
    @test isapprox(coeff_minus.α0, 3 / 4; rtol=1e-9)

    exact_unit = remez(
        0,
        1,
        0,
        0,
        "0.1",
        "10";
        algorithm=:refined,
        precision=96,
        tolerance="1e-10",
    )
    @test exact_unit.numerator == BigFloat[1]
    @test iszero(exact_unit.maximum_relative_error)

    exact_linear = remez(
        1,
        1,
        1,
        0,
        "0.1",
        "10";
        algorithm=:refined,
        precision=128,
        tolerance="1e-20",
    )
    @test exact_linear.numerator == BigFloat[0, 1]
    @test exact_linear.denominator == BigFloat[1]
    @test iszero(exact_linear.maximum_relative_error)
    @test iszero(relative_error(exact_linear, big"3.25"))

    unequal_degrees = remez(
        1,
        2,
        2,
        1,
        "1",
        "4";
        algorithm=:refined,
        precision=128,
        tolerance="1e-10",
    )
    @test length(unequal_degrees.extrema) == 5
    @test all(
        signbit(unequal_degrees.extrema_errors[i]) !=
        signbit(unequal_degrees.extrema_errors[i + 1])
        for i in 1:(length(unequal_degrees.extrema_errors) - 1)
    )
    @test_throws ArgumentError partial_fraction(unequal_degrees)
end

@testset "power times exponential target" begin
    approximation = remez(
        1,
        2,
        1,
        1,
        "1",
        "2";
        algorithm=:refined,
        precision=128,
        tolerance="1e-10",
        exponential_coefficients=["0.1", "-0.02"],
        exponential_powers=[1, -1],
    )

    @test all(
        isapprox(value, expected; rtol=big"1e-35") for (value, expected) in zip(
            approximation.exponential_coefficients,
            BigFloat[big"0.1", big"-0.02"],
        )
    )
    @test approximation.exponential_powers == [1, -1]
    for x in BigFloat[big"1", big"1.5", big"2"]
        expected = sqrt(x) * exp(big"0.1" * x - big"0.02" / x)
        @test isapprox(target_value(approximation, x), expected; rtol=big"1e-35")
        @test isapprox(
            inverse_target_value(approximation, x),
            inv(expected);
            rtol=big"1e-35",
        )
        @test isapprox(
            inverse_value(approximation, x) * approximation(x),
            big"1";
            rtol=big"1e-35",
        )
    end

    magnitudes = abs.(approximation.extrema_errors)
    spread = (maximum(magnitudes) - minimum(magnitudes)) / maximum(magnitudes)
    @test spread <= big"1e-10"
    @test all(
        signbit(approximation.extrema_errors[i]) !=
        signbit(approximation.extrema_errors[i + 1])
        for i in 1:(length(approximation.extrema_errors) - 1)
    )
end

@testset "native (5,5) minimax certificate and JLL differential" begin
    reference = STANDARD_REFERENCE
    native = remez(
        1,
        2,
        5,
        5,
        "0.0004",
        "64";
        algorithm=:refined,
        precision=256,
        tolerance="1e-15",
        max_iterations=2_000,
    )

    @test native.numerator[end] != 0
    @test native.denominator[end] == 1
    @test length(native.error_zeros) == 11
    @test length(native.extrema) == 12
    @test issorted(native.error_zeros)
    @test issorted(native.extrema)
    @test maximum(abs(relative_error(native, x)) for x in native.error_zeros) < big"1e-45"
    @test all(
        signbit(native.extrema_errors[i]) != signbit(native.extrema_errors[i + 1])
        for i in 1:(length(native.extrema_errors) - 1)
    )

    magnitudes = abs.(native.extrema_errors)
    amplitude_spread = (maximum(magnitudes) - minimum(magnitudes)) / maximum(magnitudes)
    @test amplitude_spread <= big"1e-12"
    @test isapprox(
        native.maximum_relative_error,
        big"0.00255082266915";
        rtol=big"1e-10",
    )

    xs = log_grid(STANDARD_INPUT.lower, STANDARD_INPUT.upper, 2_001)
    sampled_maximum = maximum(abs(relative_error(native, x)) for x in xs)
    @test sampled_maximum <= native.maximum_relative_error * (1 + big"1e-12")
    @test isapprox(sampled_maximum, native.maximum_relative_error; rtol=big"1e-9")

    differential_error = maximum(
        abs(native(x) - BigFloat(reference.positive(x))) / sqrt(BigFloat(x)) for x in xs
    )
    @test differential_error < big"5e-10"

    positive_pfe = partial_fraction(native; algorithm=:refined)
    negative_pfe = partial_fraction(native; inverse=true, algorithm=:refined)
    @test length(positive_pfe) == 5
    @test length(negative_pfe) == 5
    @test issorted(positive_pfe.shifts)
    @test issorted(negative_pfe.shifts)
    @test all(>(0), positive_pfe.shifts)
    @test all(>(0), negative_pfe.shifts)

    coefficient_pfe_defect = maximum(abs(positive_pfe(x) / native(x) - 1) for x in xs)
    inverse_pfe_defect = maximum(abs(negative_pfe(x) * native(x) - 1) for x in xs)
    @test coefficient_pfe_defect < big"1e-40"
    @test inverse_pfe_defect < big"1e-40"

    @test isapprox(
        positive_pfe.constant,
        BigFloat(reference.positive.constant);
        rtol=big"1e-9",
    )
    @test isapprox(
        negative_pfe.constant,
        BigFloat(reference.negative.constant);
        rtol=big"1e-9",
    )

    compatible = AlgRemez_coeffs(positive_pfe)
    @test compatible.n == 5
    @test compatible.α0 isa Float64
    @test compatible.α isa Vector{Float64}
    @test compatible.β isa Vector{Float64}
    @test isapprox(compatible(0.1), Float64(positive_pfe(big"0.1")); rtol=1e-14)

    compatible_inverse = AlgRemez_coeffs(negative_pfe)
    native_positive_values = vcat(compatible.α0, compatible.α, compatible.β)
    legacy_positive_values = vcat(
        reference.positive.constant,
        reference.positive.residues,
        reference.positive.shifts,
    )
    native_negative_values = vcat(
        compatible_inverse.α0,
        compatible_inverse.α,
        compatible_inverse.β,
    )
    legacy_negative_values = vcat(
        reference.negative.constant,
        reference.negative.residues,
        reference.negative.shifts,
    )
    @test maximum(ulp_distance.(native_positive_values, legacy_positive_values)) <= 4
    @test maximum(ulp_distance.(native_negative_values, legacy_negative_values)) <= 4
end

@testset "faithful legacy algorithm and JLL coefficient identity" begin
    reference = STANDARD_REFERENCE
    @test AlgRemez._legacy_precision_bits(42) == 192
    @test AlgRemez._legacy_precision_bits(20) == 128
    legacy_native = remez(
        STANDARD_INPUT.y,
        STANDARD_INPUT.z,
        STANDARD_INPUT.n,
        STANDARD_INPUT.d,
        STANDARD_INPUT.lower,
        STANDARD_INPUT.upper;
        algorithm=:legacy,
        precision=192,
        tolerance=1.0e-15,
    )
    @test legacy_native.iterations == reference.iterations
    @test legacy_native.precision == 192
    @test isapprox(
        Float64(legacy_native.maximum_relative_error),
        reference.reported_error;
        rtol=1.0e-6,
    )

    positive, negative = calc_coefficients(
        STANDARD_INPUT.y,
        STANDARD_INPUT.z,
        STANDARD_INPUT.n,
        STANDARD_INPUT.lower,
        STANDARD_INPUT.upper;
        precision=STANDARD_INPUT.precision,
    )
    native_positive_values = vcat(positive.α0, positive.α, positive.β)
    reference_positive_values = vcat(
        reference.positive.constant,
        reference.positive.residues,
        reference.positive.shifts,
    )
    native_negative_values = vcat(negative.α0, negative.α, negative.β)
    reference_negative_values = vcat(
        reference.negative.constant,
        reference.negative.residues,
        reference.negative.shifts,
    )
    @test maximum(ulp_distance.(native_positive_values, reference_positive_values)) <= 1
    @test maximum(ulp_distance.(native_negative_values, reference_negative_values)) <= 1
end

@testset "legacy file and command-line interfaces" begin
    approximation = remez(1, 2, 0, 0, "1", "4"; precision=96, tolerance="1e-10")
    mktempdir() do directory
        paths = write_legacy_files(approximation; directory=directory, ratio=2)
        @test paths.approximation == joinpath(directory, "approx.dat")
        @test paths.error == joinpath(directory, "error.dat")
        approximation_text = read(paths.approximation, String)
        error_lines = readlines(paths.error)
        @test count(contains("Approximation to f(x)"), eachline(IOBuffer(approximation_text))) == 2
        @test count(contains("alpha[0]"), eachline(IOBuffer(approximation_text))) == 2
        @test length(error_lines) == 2
        @test all(length(split(line)) == 2 for line in error_lines)
    end

    mktempdir() do directory
        command_result = redirect_stdout(devnull) do
            legacy_main(["1", "2", "0", "0", "1", "4", "20"]; directory=directory)
        end
        @test command_result isa RationalApproximation
        @test isfile(joinpath(directory, "approx.dat"))
        @test isfile(joinpath(directory, "error.dat"))
    end
end

include("rhmc.jl")

@testset "AlgRemez_jll reference" begin
    reference = STANDARD_REFERENCE

    @test length(reference.positive) == STANDARD_INPUT.n
    @test length(reference.negative) == STANDARD_INPUT.n
    @test length(reference.error_samples) > 1_000
    @test reference.iterations > 0
    @test isapprox(reference.reported_error, 2.550823e-3; rtol=1.0e-6)

    @test all(isfinite, reference.positive.residues)
    @test all(isfinite, reference.positive.shifts)
    @test all(isfinite, reference.negative.residues)
    @test all(isfinite, reference.negative.shifts)
    @test issorted(reference.positive.shifts)
    @test issorted(reference.negative.shifts)

    # Coarse golden values detect an unexpected change of the packaged legacy
    # reference while allowing harmless last-bit differences across platforms.
    @test isapprox(reference.positive.constant, 18.735935560108796; rtol=1.0e-12)
    @test isapprox(reference.negative.constant, 0.053373368882049743; rtol=1.0e-12)
    @test isapprox(first(reference.positive.shifts), 1.2694829377677037e-3; rtol=1.0e-11)
    @test isapprox(first(reference.negative.shifts), 2.0859203185118118e-4; rtol=1.0e-11)
end

@testset "legacy (5,5) regression — not an accuracy target" begin
    reference = STANDARD_REFERENCE
    xs = log_grid(STANDARD_INPUT.lower, STANDARD_INPUT.upper, 2_001)

    relative_errors = [reference.positive(x) / sqrt(x) - 1 for x in xs]
    degree_five_reference_error = maximum(abs, relative_errors)

    # This deliberately records the legacy low-degree result.  It does not
    # declare 2.55e-3 to be an acceptable production accuracy.  Accuracy
    # acceptance tests must state a target and choose a sufficient degree.
    @test 2.54e-3 < degree_five_reference_error < 2.56e-3

    # Both PFE sections are generated from the same rational approximation and
    # must reconstruct reciprocal functions.  This checks the parser and the
    # root/residue post-processing without requiring coefficient equality.
    reciprocal_defect = maximum(
        abs(reference.positive(x) * reference.negative(x) - 1) for x in xs
    )
    @test reciprocal_defect < 5.0e-12

    # The legacy error grid crosses all K = n+d+1 interpolation zeros.  It is
    # not accurate enough to certify extrema, but it is useful as a smoke test.
    @test sign_change_count(reference.error_samples) ==
          STANDARD_INPUT.n + STANDARD_INPUT.d + 1

    file_error_defect = maximum(
        abs(reference.positive(sample.x) / sqrt(sample.x) - 1 - sample.error)
        for sample in reference.error_samples
    )
    @test file_error_defect < 5.0e-9
end

@testset "reference runner rejects unsafe legacy inputs" begin
    @test_throws ArgumentError jll_reference(1, 2, 0, 0, 1.0, 4.0)
    @test_throws ArgumentError jll_reference(1, 2, 2, 1, 1.0e-3, 1.0)
    @test_throws ArgumentError jll_reference(1, 0, 2, 2, 1.0e-3, 1.0)
    @test_throws ArgumentError jll_reference(1, 2, 2, 2, 0.0, 1.0)
end

@testset "native input validation" begin
    @test_throws ArgumentError remez(-1, 2, 2, 2, 1, 4)
    @test_throws ArgumentError remez(1, 0, 2, 2, 1, 4)
    @test_throws ArgumentError remez(1, 2, -1, 2, 1, 4)
    @test_throws ArgumentError remez(1, 2, 2, -1, 1, 4)
    @test_throws ArgumentError remez(1, 2, 2, 2, 0, 4)
    @test_throws ArgumentError remez(1, 2, 2, 2, 4, 1)
    @test_throws ArgumentError remez(1, 2, 2, 2, 1, 4; precision=32)
    @test_throws ArgumentError remez(1, 2, 2, 2, 1, 4; tolerance=2)
    @test_throws ArgumentError remez(1, 2, 2, 2, 1, 4; algorithm=:unknown)
    @test_throws DimensionMismatch remez(
        1,
        2,
        2,
        2,
        1,
        4;
        exponential_coefficients=[1],
        exponential_powers=Int[],
    )
    @test_throws DimensionMismatch AlgRemez_coeffs(1, [2], [3, 4], 1)
    @test_throws ArgumentError legacy_main(String[])
end
