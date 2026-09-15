# This fixture reproduces the coefficient ownership and sign dispatch in
# LatticeDiracOperators.Rhmc.RHMC without adding that large package as a test
# dependency. The matrix application below mirrors its shiftedcg + weighted-sum
# contract, with exact shifted solves so this test isolates the Remez/PFE layer.
struct _RHMCFixture
    y::Int64
    z::Int64
    coeffs::AlgRemez_coeffs
    coeffs_inverse::AlgRemez_coeffs
end

function _rhmc_fixture(
    order::Rational,
    positive::AlgRemez_coeffs,
    negative::AlgRemez_coeffs,
)
    y = numerator(order)
    z = denominator(order)
    if y > 0
        return _RHMCFixture(y, z, positive, negative)
    end
    return _RHMCFixture(y, z, negative, positive)
end

function _apply_rhmc_pfe(
    coefficients::AlgRemez_coeffs,
    operator::AbstractMatrix,
    source::AbstractVector,
)
    result = coefficients.α0 .* source
    for i in 1:coefficients.n
        shifted_solution = (operator + coefficients.β[i] * I) \ source
        result .+= coefficients.α[i] .* shifted_solution
    end
    return result
end

@testset "LatticeQCD RHMC production profiles" begin
    lower = 4.0e-4
    upper = 64.0
    spectrum = exp.(range(log(lower), log(upper); length=513))
    operator = Diagonal(spectrum)
    source = [1 + 0.1 * sin(i) for i in eachindex(spectrum)]

    # LatticeQCD uses x^(Nf/16), degree 15 for the accept/reject action and
    # x^(Nf/8), degree 10 for molecular-dynamics evolution when Nf = 1,2,3.
    profiles = [
        (nf=nf, role=:action, power=nf // 16, degree=15, error_limit=5.0e-9)
        for nf in 1:3
    ]
    append!(
        profiles,
        [
            (nf=nf, role=:md, power=nf // 8, degree=10, error_limit=5.0e-6)
            for nf in 1:3
        ],
    )

    for profile in profiles
        y = numerator(profile.power)
        z = denominator(profile.power)
        # Exercise the unchanged LatticeDiracOperators call signature and its
        # default strict exchange tolerance.
        positive, negative = calc_coefficients(
            y,
            z,
            profile.degree,
            lower,
            upper;
            precision=42,
        )
        legacy = jll_reference(
            y,
            z,
            profile.degree,
            profile.degree,
            lower,
            upper;
            precision=42,
        )

        @test positive isa AlgRemez_coeffs
        @test negative isa AlgRemez_coeffs
        @test positive.n == profile.degree
        @test negative.n == profile.degree
        @test all(isfinite, positive.α)
        @test all(isfinite, positive.β)
        @test all(isfinite, negative.α)
        @test all(isfinite, negative.β)
        @test all(>(0), positive.β)
        @test all(>(0), negative.β)

        native_positive_values = vcat(positive.α0, positive.α, positive.β)
        legacy_positive_values = vcat(
            legacy.positive.constant,
            legacy.positive.residues,
            legacy.positive.shifts,
        )
        native_negative_values = vcat(negative.α0, negative.α, negative.β)
        legacy_negative_values = vcat(
            legacy.negative.constant,
            legacy.negative.residues,
            legacy.negative.shifts,
        )
        @test maximum(
            ulp_distance.(native_positive_values, legacy_positive_values),
        ) <= 4
        @test maximum(
            ulp_distance.(native_negative_values, legacy_negative_values),
        ) <= 4
        @test maximum(abs(positive(x) / legacy.positive(x) - 1) for x in spectrum) <
              2.0e-13
        @test maximum(abs(negative(x) / legacy.negative(x) - 1) for x in spectrum) <
              2.0e-13

        rhmc = _rhmc_fixture(profile.power, positive, negative)
        exact_positive = @. spectrum^(y / z) * source
        exact_negative = @. spectrum^(-y / z) * source
        applied_positive = _apply_rhmc_pfe(rhmc.coeffs, operator, source)
        applied_negative = _apply_rhmc_pfe(rhmc.coeffs_inverse, operator, source)
        positive_error = maximum(abs.(applied_positive ./ exact_positive .- 1))
        negative_error = maximum(abs.(applied_negative ./ exact_negative .- 1))
        @test positive_error < profile.error_limit
        @test negative_error < profile.error_limit

        # The inverse PFE must undo the forward PFE mode by mode. This is the
        # coefficient pairing used between pseudofermion sampling and action
        # evaluation in RHMC.
        round_trip = _apply_rhmc_pfe(
            rhmc.coeffs_inverse,
            operator,
            applied_positive,
        )
        @test maximum(abs.(round_trip ./ source .- 1)) < 5.0e-12

        # A negative Rational order is represented by swapping the exact same
        # two coefficient objects, as in LatticeDiracOperators.Rhmc.RHMC.
        inverse_rhmc = _rhmc_fixture(-profile.power, positive, negative)
        @test inverse_rhmc.coeffs === negative
        @test inverse_rhmc.coeffs_inverse === positive
        applied_swapped = _apply_rhmc_pfe(inverse_rhmc.coeffs, operator, source)
        @test maximum(abs.(applied_swapped ./ exact_negative .- 1)) <
              profile.error_limit
    end
end
