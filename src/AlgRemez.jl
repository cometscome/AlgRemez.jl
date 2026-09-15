"""
    AlgRemez

Arbitrary-precision rational minimax approximation using a Remez exchange
algorithm.
"""
module AlgRemez

using LinearAlgebra
using Printf

export AlgRemez_coeffs,
       PartialFractionApproximation,
       RationalApproximation,
       RemezConvergenceError,
       calc_coefficients,
       fittedfunction,
       inverse_target_value,
       inverse_value,
       legacy_main,
       partial_fraction,
       relative_error,
       remez,
       target_value,
       write_legacy_files

"""Raised when the exchange iteration cannot meet its requested tolerance."""
struct RemezConvergenceError{T} <: Exception
    iterations::Int
    spread::T
end

function Base.showerror(io::IO, error::RemezConvergenceError)
    print(
        io,
        "Remez exchange did not converge after ",
        error.iterations,
        " iterations (error-amplitude spread = ",
        error.spread,
        ')',
    )
end

"""
    RationalApproximation

Result of a power or power-times-exponential Remez approximation. Polynomial
coefficients are in ascending order, so `numerator[j + 1]` and
`denominator[j + 1]` multiply `x^j`. The denominator is monic.
"""
struct RationalApproximation{T<:AbstractFloat}
    numerator::Vector{T}
    denominator::Vector{T}
    power_num::Int
    power_den::Int
    exponential_coefficients::Vector{T}
    exponential_powers::Vector{Int}
    lower::T
    upper::T
    maximum_relative_error::T
    iterations::Int
    error_zeros::Vector{T}
    extrema::Vector{T}
    extrema_errors::Vector{T}
    precision::Int

    function RationalApproximation(
        numerator::Vector{T},
        denominator::Vector{T},
        power_num::Integer,
        power_den::Integer,
        exponential_coefficients::Vector{T},
        exponential_powers::Vector{Int},
        lower::T,
        upper::T,
        maximum_relative_error::T,
        iterations::Integer,
        error_zeros::Vector{T},
        extrema::Vector{T},
        extrema_errors::Vector{T},
        precision::Integer,
    ) where {T<:AbstractFloat}
        isempty(numerator) && throw(ArgumentError("the numerator cannot be empty"))
        isempty(denominator) && throw(ArgumentError("the denominator cannot be empty"))
        length(extrema) == length(extrema_errors) ||
            throw(DimensionMismatch("extrema and extrema error counts differ"))
        length(exponential_coefficients) == length(exponential_powers) ||
            throw(DimensionMismatch("exponential coefficient and power counts differ"))
        return new{T}(
            numerator,
            denominator,
            Int(power_num),
            Int(power_den),
            exponential_coefficients,
            exponential_powers,
            lower,
            upper,
            maximum_relative_error,
            Int(iterations),
            error_zeros,
            extrema,
            extrema_errors,
            Int(precision),
        )
    end
end

Base.broadcastable(approximation::RationalApproximation) = Ref(approximation)

"""Evaluate a polynomial whose coefficients are stored in ascending order."""
function _evaluate_polynomial(x, coefficients::AbstractVector)
    value = coefficients[end]
    for i in (lastindex(coefficients) - 1):-1:firstindex(coefficients)
        value = muladd(value, x, coefficients[i])
    end
    return value
end

function (approximation::RationalApproximation{T})(x::Real) where {T}
    xx = convert(promote_type(T, typeof(x)), x)
    numerator = _evaluate_polynomial(xx, approximation.numerator)
    denominator = _evaluate_polynomial(xx, approximation.denominator)
    return numerator / denominator
end

function (approximation::RationalApproximation{BigFloat})(x::Real)
    return setprecision(BigFloat, approximation.precision) do
        xx = BigFloat(x)
        numerator = _evaluate_polynomial(xx, approximation.numerator)
        denominator = _evaluate_polynomial(xx, approximation.denominator)
        return numerator / denominator
    end
end

"""Evaluate the target function represented by an approximation result."""
function target_value(approximation::RationalApproximation{T}, x::Real) where {T}
    S = promote_type(T, typeof(x))
    xx = convert(S, x)
    exponent = convert(S, approximation.power_num) / convert(S, approximation.power_den)
    value = xx == one(S) ? one(S) : xx^exponent
    if !isempty(approximation.exponential_coefficients)
        exponent_sum = zero(S)
        for i in eachindex(
            approximation.exponential_coefficients,
            approximation.exponential_powers,
        )
            exponent_sum += convert(S, approximation.exponential_coefficients[i]) *
                            xx^approximation.exponential_powers[i]
        end
        value *= exp(exponent_sum)
    end
    return value
end

function target_value(approximation::RationalApproximation{BigFloat}, x::Real)
    return setprecision(BigFloat, approximation.precision) do
        _target(
            BigFloat(x),
            approximation.power_num,
            approximation.power_den,
            approximation.exponential_coefficients,
            approximation.exponential_powers,
        )
    end
end

"""Evaluate the reciprocal of the rational approximation."""
inverse_value(approximation::RationalApproximation, x::Real) = inv(approximation(x))

"""Evaluate the reciprocal of the target function."""
inverse_target_value(approximation::RationalApproximation, x::Real) =
    inv(target_value(approximation, x))

"""Return the signed relative error `R(x) / f(x) - 1`."""
function relative_error(approximation::RationalApproximation, x::Real)
    target = target_value(approximation, x)
    return approximation(x) / target - one(target)
end

function relative_error(approximation::RationalApproximation{BigFloat}, x::Real)
    return setprecision(BigFloat, approximation.precision) do
        _relative_error(
            BigFloat(x),
            approximation.numerator,
            approximation.denominator,
            approximation.power_num,
            approximation.power_den,
            approximation.exponential_coefficients,
            approximation.exponential_powers,
        )
    end
end

"""
    AlgRemez_coeffs(α0, α, β, n)

Float64 partial-fraction coefficients compatible with the former
`AlgRemez_jll` wrapper used by LatticeDiracOperators.jl.
"""
struct AlgRemez_coeffs
    α0::Float64
    α::Vector{Float64}
    β::Vector{Float64}
    n::Int64

    function AlgRemez_coeffs(α0::Real, α::AbstractVector, β::AbstractVector, n::Integer)
        n >= 0 || throw(ArgumentError("n must be nonnegative"))
        length(α) == length(β) == n ||
            throw(DimensionMismatch("α and β must both contain n entries"))
        return new(Float64(α0), Float64.(α), Float64.(β), Int64(n))
    end
end

Base.length(coefficients::AlgRemez_coeffs) = coefficients.n
Base.broadcastable(coefficients::AlgRemez_coeffs) = Ref(coefficients)

function Base.show(io::IO, ::MIME"text/plain", coefficients::AlgRemez_coeffs)
    println(io, "f(x) = α0 + sum(α[i] / (x + β[i]), i=1:n)")
    println(io, "n: ", coefficients.n)
    println(io, "α0: ", coefficients.α0)
    println(io, "α: ", coefficients.α)
    print(io, "β: ", coefficients.β)
end

function (coefficients::AlgRemez_coeffs)(x::Real)
    value = coefficients.α0
    for i in 1:coefficients.n
        value += coefficients.α[i] / (x + coefficients.β[i])
    end
    return value
end

"""Return a callable closure for legacy `AlgRemez_coeffs`."""
fittedfunction(coefficients::AlgRemez_coeffs) = x -> coefficients(x)

"""
    PartialFractionApproximation

Representation `constant + sum(residues[i] / (x + shifts[i]))`.
"""
struct PartialFractionApproximation{T<:AbstractFloat}
    constant::T
    residues::Vector{T}
    shifts::Vector{T}

    function PartialFractionApproximation(
        constant::T,
        residues::Vector{T},
        shifts::Vector{T},
    ) where {T<:AbstractFloat}
        length(residues) == length(shifts) ||
            throw(DimensionMismatch("residue and shift counts differ"))
        return new{T}(constant, residues, shifts)
    end
end

Base.length(expansion::PartialFractionApproximation) = length(expansion.residues)
Base.broadcastable(expansion::PartialFractionApproximation) = Ref(expansion)

function (expansion::PartialFractionApproximation{T})(x::Real) where {T}
    S = promote_type(T, typeof(x))
    xx = convert(S, x)
    value = convert(S, expansion.constant)
    for i in eachindex(expansion.residues, expansion.shifts)
        value += convert(S, expansion.residues[i]) /
                 (xx + convert(S, expansion.shifts[i]))
    end
    return value
end

function (expansion::PartialFractionApproximation{BigFloat})(x::Real)
    return setprecision(BigFloat, precision(expansion.constant)) do
        xx = BigFloat(x)
        value = expansion.constant
        for i in eachindex(expansion.residues, expansion.shifts)
            value += expansion.residues[i] / (xx + expansion.shifts[i])
        end
        return value
    end
end

AlgRemez_coeffs(expansion::PartialFractionApproximation) = AlgRemez_coeffs(
    expansion.constant,
    expansion.residues,
    expansion.shifts,
    length(expansion),
)

"""Return a callable closure for a native partial-fraction expansion."""
fittedfunction(expansion::PartialFractionApproximation) = x -> expansion(x)

_as_bigfloat(x::AbstractString) = parse(BigFloat, x)
_as_bigfloat(x::Real) = BigFloat(x)

function _target(
    x::BigFloat,
    power_num::Int,
    power_den::Int,
    exponential_coefficients::Vector{BigFloat},
    exponential_powers::Vector{Int},
)
    value = x == one(x) ? one(x) : x^(BigFloat(power_num) / BigFloat(power_den))
    if !isempty(exponential_coefficients)
        exponent_sum = zero(BigFloat)
        for i in eachindex(exponential_coefficients, exponential_powers)
            exponent_sum += exponential_coefficients[i] * x^exponential_powers[i]
        end
        value *= exp(exponent_sum)
    end
    return value
end

function _relative_error(
    x::BigFloat,
    numerator::Vector{BigFloat},
    denominator::Vector{BigFloat},
    power_num::Int,
    power_den::Int,
    exponential_coefficients::Vector{BigFloat},
    exponential_powers::Vector{Int},
)
    rational = _evaluate_polynomial(x, numerator) / _evaluate_polynomial(x, denominator)
    target = _target(
        x,
        power_num,
        power_den,
        exponential_coefficients,
        exponential_powers,
    )
    return rational / target - one(BigFloat)
end

function _initial_points(lower::BigFloat, upper::BigFloat, count::Int)
    log_lower = log(lower)
    log_width = log(upper) - log_lower
    pi_value = BigFloat(pi)

    extrema = Vector{BigFloat}(undef, count + 1)
    extrema[1] = lower
    for i in 1:(count - 1)
        t = (one(BigFloat) - cos(pi_value * BigFloat(i) / BigFloat(count))) / 2
        extrema[i + 1] = exp(log_lower + log_width * t)
    end
    extrema[end] = upper

    zeros = Vector{BigFloat}(undef, count)
    for i in 0:(count - 1)
        angle = pi_value * BigFloat(2i + 1) / BigFloat(2count)
        t = (one(BigFloat) - cos(angle)) / 2
        zeros[i + 1] = exp(log_lower + log_width * t)
    end
    return zeros, extrema
end

function _solve_interpolation(
    zeros::Vector{BigFloat},
    numerator_degree::Int,
    denominator_degree::Int,
    power_num::Int,
    power_den::Int,
    exponential_coefficients::Vector{BigFloat},
    exponential_powers::Vector{Int},
)
    equation_count = numerator_degree + denominator_degree + 1
    matrix = Matrix{BigFloat}(undef, equation_count, equation_count)
    right_hand_side = Vector{BigFloat}(undef, equation_count)

    for row in 1:equation_count
        x = zeros[row]
        target = _target(
            x,
            power_num,
            power_den,
            exponential_coefficients,
            exponential_powers,
        )

        power = one(BigFloat)
        for degree in 0:numerator_degree
            matrix[row, degree + 1] = power
            power *= x
        end

        power = one(BigFloat)
        for degree in 0:(denominator_degree - 1)
            matrix[row, numerator_degree + 2 + degree] = -target * power
            power *= x
        end
        right_hand_side[row] = target * power
    end

    parameters = matrix \ right_hand_side
    numerator = parameters[1:(numerator_degree + 1)]
    denominator = vcat(
        parameters[(numerator_degree + 2):end],
        one(BigFloat),
    )
    return numerator, denominator
end

function _golden_maximum(
    objective,
    left::BigFloat,
    right::BigFloat,
    position_tolerance::BigFloat,
    maximum_iterations::Int,
)
    left == right && return left, objective(left)

    inverse_phi = (sqrt(BigFloat(5)) - one(BigFloat)) / 2
    a = left
    b = right
    c = b - inverse_phi * (b - a)
    d = a + inverse_phi * (b - a)
    fc = objective(c)
    fd = objective(d)

    for _ in 1:maximum_iterations
        scale = max(abs(a), abs(b), sqrt(eps(BigFloat)))
        b - a <= position_tolerance * scale && break
        if fc < fd
            a = c
            c = d
            fc = fd
            d = a + inverse_phi * (b - a)
            fd = objective(d)
        else
            b = d
            d = c
            fd = fc
            c = b - inverse_phi * (b - a)
            fc = objective(c)
        end
    end

    candidates = (left, c, d, right)
    values = map(objective, candidates)
    best = argmax(values)
    return candidates[best], values[best]
end

function _find_extrema(
    zeros::Vector{BigFloat},
    lower::BigFloat,
    upper::BigFloat,
    numerator::Vector{BigFloat},
    denominator::Vector{BigFloat},
    power_num::Int,
    power_den::Int,
    exponential_coefficients::Vector{BigFloat},
    exponential_powers::Vector{Int},
    position_tolerance::BigFloat,
    maximum_iterations::Int,
)
    boundaries = vcat(lower, zeros, upper)
    extrema = Vector{BigFloat}(undef, length(boundaries) - 1)
    errors = similar(extrema)

    signed_error(x) = _relative_error(
        x,
        numerator,
        denominator,
        power_num,
        power_den,
        exponential_coefficients,
        exponential_powers,
    )

    for i in eachindex(extrema)
        location, _ = _golden_maximum(
            x -> abs(signed_error(x)),
            boundaries[i],
            boundaries[i + 1],
            position_tolerance,
            maximum_iterations,
        )
        extrema[i] = location
        errors[i] = signed_error(location)
    end
    return extrema, errors
end

function _alternating_signs(errors::Vector{BigFloat})
    all(!iszero, errors) || return false
    return all(signbit(errors[i]) != signbit(errors[i + 1]) for i in 1:(length(errors) - 1))
end

function _default_tolerance(precision::Int)
    decimal_digits = floor(Int, precision * log10(2))
    requested_digits = clamp(decimal_digits ÷ 3, 12, 24)
    return BigFloat(10)^(-requested_digits)
end

"""
    _refined_remez(y, z, n, d, lower, upper; kwargs...) -> RationalApproximation

Construct a relative-error minimax rational approximation of type `(n, d)` to
`x^(y/z) * exp(sum(a[j] * x^k[j]))` on the positive interval
`[lower, upper]`.

Keyword arguments:

  * `precision=256`: working precision in bits;
  * `tolerance=nothing`: relative spread allowed among error extrema;
  * `max_iterations=10_000`: maximum exchange iterations;
  * `extrema_iterations=160`: maximum golden-section iterations per interval.
  * `exponential_coefficients=[]`: the coefficients `a` of an optional
    exponential factor;
  * `exponential_powers=[]`: the corresponding integer powers `k`.

String or `BigFloat` endpoints avoid an initial `Float64` rounding. The result
retains `BigFloat` coefficients and a minimax certificate consisting of the
error zeros, extrema, and signed errors at those extrema.
"""
function _refined_remez(
    y::Integer,
    z::Integer,
    n::Integer,
    d::Integer,
    lower::Union{Real,AbstractString},
    upper::Union{Real,AbstractString};
    precision::Integer=256,
    tolerance::Union{Nothing,Real,AbstractString}=nothing,
    max_iterations::Integer=10_000,
    extrema_iterations::Integer=160,
    exponential_coefficients=nothing,
    exponential_powers=nothing,
)
    y >= 0 || throw(ArgumentError("the exponent numerator must be nonnegative"))
    z > 0 || throw(ArgumentError("the exponent denominator must be positive"))
    n >= 0 || throw(ArgumentError("the numerator degree must be nonnegative"))
    d >= 0 || throw(ArgumentError("the denominator degree must be nonnegative"))
    precision >= 64 || throw(ArgumentError("precision must be at least 64 bits"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    extrema_iterations > 0 || throw(ArgumentError("extrema_iterations must be positive"))

    return setprecision(BigFloat, Int(precision)) do
        lower_big = _as_bigfloat(lower)
        upper_big = _as_bigfloat(upper)
        0 < lower_big < upper_big || throw(ArgumentError("expected 0 < lower < upper"))

        coefficients_big = exponential_coefficients === nothing ? BigFloat[] :
                           _as_bigfloat.(collect(exponential_coefficients))
        powers_int = exponential_powers === nothing ? Int[] :
                     Int.(collect(exponential_powers))
        length(coefficients_big) == length(powers_int) || throw(
            DimensionMismatch(
                "exponential_coefficients and exponential_powers must have equal length",
            ),
        )
        all(isfinite, coefficients_big) ||
            throw(ArgumentError("exponential coefficients must be finite"))

        target_tolerance = tolerance === nothing ?
                           _default_tolerance(Int(precision)) : _as_bigfloat(tolerance)
        0 < target_tolerance < 1 ||
            throw(ArgumentError("tolerance must be between zero and one"))
        position_tolerance = sqrt(target_tolerance) / 10

        power_num = Int(y)
        power_den = Int(z)
        equation_count = Int(n + d + 1)
        zeros, _ = _initial_points(lower_big, upper_big, equation_count)

        relaxation = BigFloat(0.25)
        previous_spread = BigFloat(Inf)
        last_spread = previous_spread

        for iteration in 1:Int(max_iterations)
            numerator, denominator = _solve_interpolation(
                zeros,
                Int(n),
                Int(d),
                power_num,
                power_den,
                coefficients_big,
                powers_int,
            )
            extrema, extrema_errors = _find_extrema(
                zeros,
                lower_big,
                upper_big,
                numerator,
                denominator,
                power_num,
                power_den,
                coefficients_big,
                powers_int,
                position_tolerance,
                Int(extrema_iterations),
            )
            magnitudes = abs.(extrema_errors)
            maximum_error = maximum(magnitudes)
            minimum_error = minimum(magnitudes)

            if iszero(maximum_error)
                return RationalApproximation(
                    numerator,
                    denominator,
                    power_num,
                    power_den,
                    copy(coefficients_big),
                    copy(powers_int),
                    lower_big,
                    upper_big,
                    maximum_error,
                    iteration,
                    copy(zeros),
                    extrema,
                    extrema_errors,
                    Int(precision),
                )
            end

            last_spread = (maximum_error - minimum_error) / maximum_error
            if last_spread <= target_tolerance && _alternating_signs(extrema_errors)
                return RationalApproximation(
                    numerator,
                    denominator,
                    power_num,
                    power_den,
                    copy(coefficients_big),
                    copy(powers_int),
                    lower_big,
                    upper_big,
                    maximum_error,
                    iteration,
                    copy(zeros),
                    extrema,
                    extrema_errors,
                    Int(precision),
                )
            end

            if last_spread >= previous_spread
                relaxation /= 2
                relaxation > sqrt(eps(BigFloat)) ||
                    throw(RemezConvergenceError(iteration, last_spread))
            end
            previous_spread = last_spread

            new_zeros = similar(zeros)
            for i in eachindex(zeros)
                ratio = iszero(magnitudes[i + 1]) ? BigFloat(0.25) :
                        magnitudes[i] / magnitudes[i + 1] - one(BigFloat)
                ratio = min(ratio, BigFloat(0.25))
                step = relaxation * ratio * (extrema[i + 1] - extrema[i])
                candidate = zeros[i] - step
                if candidate <= extrema[i]
                    candidate = (extrema[i] + zeros[i]) / 2
                elseif candidate >= extrema[i + 1]
                    candidate = (extrema[i + 1] + zeros[i]) / 2
                end
                new_zeros[i] = candidate
            end
            zeros = new_zeros
        end
        throw(RemezConvergenceError(Int(max_iterations), last_spread))
    end
end

function _legacy_target(
    x::BigFloat,
    power_num::Int,
    power_den::Int,
    exponential_coefficients::Vector{BigFloat},
    exponential_powers::Vector{Int},
)
    exponent = BigFloat(power_num) / BigFloat(power_den)
    value = x == one(BigFloat) ? one(BigFloat) : x^exponent
    if !isempty(exponential_coefficients)
        exponent_sum = zero(BigFloat)
        for i in eachindex(exponential_coefficients, exponential_powers)
            exponent_sum += exponential_coefficients[i] * x^exponential_powers[i]
        end
        value *= exp(exponent_sum)
    end
    return value
end

function _legacy_approximation_value(
    x::BigFloat,
    parameters::Vector{BigFloat},
    numerator_degree::Int,
    denominator_degree::Int,
)
    numerator = parameters[numerator_degree + 1]
    for degree in (numerator_degree - 1):-1:0
        numerator = x * numerator + parameters[degree + 1]
    end

    if iszero(denominator_degree)
        denominator = one(BigFloat)
    else
        denominator = x + parameters[numerator_degree + denominator_degree + 1]
        for index in (numerator_degree + denominator_degree - 1):-1:(numerator_degree + 1)
            denominator = x * denominator + parameters[index + 1]
        end
    end
    return numerator / denominator
end

function _legacy_signed_error(
    x::BigFloat,
    parameters::Vector{BigFloat},
    numerator_degree::Int,
    denominator_degree::Int,
    power_num::Int,
    power_den::Int,
    exponential_coefficients::Vector{BigFloat},
    exponential_powers::Vector{Int},
)
    target = _legacy_target(
        x,
        power_num,
        power_den,
        exponential_coefficients,
        exponential_powers,
    )
    error = _legacy_approximation_value(
        x,
        parameters,
        numerator_degree,
        denominator_degree,
    ) - target
    iszero(target) || (error /= target)
    return error
end

"""Scaled partial-pivot Gaussian elimination used by the original C++ code."""
function _legacy_simq!(
    matrix::Matrix{BigFloat},
    right_hand_side::Vector{BigFloat},
)
    count = length(right_hand_side)
    size(matrix) == (count, count) || throw(DimensionMismatch("matrix must be square"))
    permutation = collect(1:count)
    solution = Vector{BigFloat}(undef, count)

    for row in 1:count
        row_norm = zero(BigFloat)
        for column in 1:count
            row_norm = max(row_norm, abs(matrix[row, column]))
        end
        iszero(row_norm) && throw(ArgumentError("legacy Gaussian elimination found a zero row"))
        solution[row] = one(BigFloat) / row_norm
    end

    if count > 1
        for column in 1:(count - 1)
            largest = zero(BigFloat)
            pivot_index = column
            for index in column:count
                row = permutation[index]
                candidate = abs(matrix[row, column]) * solution[row]
                if candidate > largest
                    largest = candidate
                    pivot_index = index
                end
            end
            iszero(largest) &&
                throw(ArgumentError("legacy Gaussian elimination found a zero pivot"))
            if pivot_index != column
                permutation[column], permutation[pivot_index] =
                    permutation[pivot_index], permutation[column]
            end

            pivot_row = permutation[column]
            pivot = matrix[pivot_row, column]
            for index in (column + 1):count
                row = permutation[index]
                elimination_multiplier = -matrix[row, column] / pivot
                matrix[row, column] = -elimination_multiplier
                for trailing_column in (column + 1):count
                    matrix[row, trailing_column] = matrix[row, trailing_column] +
                                                   elimination_multiplier *
                                                   matrix[pivot_row, trailing_column]
                end
            end
        end
    end

    iszero(matrix[permutation[end], end]) &&
        throw(ArgumentError("legacy Gaussian elimination found a zero final pivot"))

    solution[1] = right_hand_side[permutation[1]]
    for index in 2:count
        row = permutation[index]
        accumulated = zero(BigFloat)
        for column in 1:(index - 1)
            accumulated += matrix[row, column] * solution[column]
        end
        solution[index] = right_hand_side[row] - accumulated
    end

    solution[end] /= matrix[permutation[end], end]
    for index in (count - 1):-1:1
        row = permutation[index]
        accumulated = zero(BigFloat)
        for column in (index + 1):count
            accumulated += matrix[row, column] * solution[column]
        end
        solution[index] =
            (solution[index] - accumulated) / matrix[row, index]
    end
    return solution
end

function _legacy_equations(
    zeros::Vector{BigFloat},
    numerator_degree::Int,
    denominator_degree::Int,
    power_num::Int,
    power_den::Int,
    exponential_coefficients::Vector{BigFloat},
    exponential_powers::Vector{Int},
)
    equation_count = numerator_degree + denominator_degree + 1
    matrix = Matrix{BigFloat}(undef, equation_count, equation_count)
    right_hand_side = Vector{BigFloat}(undef, equation_count)

    for row in 1:equation_count
        x = zeros[row]
        target = _legacy_target(
            x,
            power_num,
            power_den,
            exponential_coefficients,
            exponential_powers,
        )
        power = one(BigFloat)
        for degree in 0:numerator_degree
            matrix[row, degree + 1] = power
            power *= x
        end
        power = one(BigFloat)
        for degree in 0:(denominator_degree - 1)
            matrix[row, numerator_degree + degree + 2] = -target * power
            power *= x
        end
        right_hand_side[row] = target * power
    end
    return _legacy_simq!(matrix, right_hand_side)
end

function _legacy_initial_points(lower::BigFloat, upper::BigFloat, equation_count::Int)
    width = upper - lower
    extrema = Vector{BigFloat}(undef, equation_count + 1)
    zeros = Vector{BigFloat}(undef, equation_count + 2)

    extrema[1] = lower
    for index in 1:(equation_count - 1)
        coordinate = 0.5 * (1 - cos(pi * index / Float64(equation_count)))
        coordinate = (exp(coordinate) - 1) / (exp(1.0) - 1)
        extrema[index + 1] = lower + BigFloat(coordinate) * width
    end
    extrema[end] = upper

    twice_count = 2.0 * equation_count
    for index in 0:equation_count
        coordinate = 0.5 * (1 - cos(pi * (2index + 1) / twice_count))
        coordinate = (exp(coordinate) - 1) / (exp(1.0) - 1)
        zeros[index + 1] = lower + BigFloat(coordinate) * width
    end
    zeros[end] = upper
    return zeros, extrema
end

function _legacy_search!(
    zeros::Vector{BigFloat},
    extrema::Vector{BigFloat},
    steps::Vector{BigFloat},
    parameters::Vector{BigFloat},
    numerator_degree::Int,
    denominator_degree::Int,
    lower::BigFloat,
    upper::BigFloat,
    power_num::Int,
    power_den::Int,
    exponential_coefficients::Vector{BigFloat},
    exponential_powers::Vector{Int},
    relaxation::BigFloat,
    previous_spread::BigFloat,
)
    equation_count = numerator_degree + denominator_degree + 1
    magnitudes = Vector{BigFloat}(undef, equation_count + 1)
    closest = BigFloat(Float64(1.0e30))
    farthest = zero(BigFloat)
    left_boundary = lower

    magnitude(x) = abs(
        _legacy_signed_error(
            x,
            parameters,
            numerator_degree,
            denominator_degree,
            power_num,
            power_den,
            exponential_coefficients,
            exponential_powers,
        ),
    )

    for index in 1:(equation_count + 1)
        right_boundary = index == equation_count + 1 ? upper : zeros[index]
        current_location = extrema[index]
        current_magnitude = magnitude(current_location)
        step = steps[index]
        next_location = current_location + step

        if next_location < left_boundary || next_location >= right_boundary
            step = -step
            next_location = current_location
            next_magnitude = current_magnitude
        else
            next_magnitude = magnitude(next_location)
            if next_magnitude < current_magnitude
                step = -step
                next_location = current_location
                next_magnitude = current_magnitude
            end
        end

        march_count = 0
        while next_magnitude >= current_magnitude
            march_count += 1
            march_count > 10 && break
            current_magnitude = next_magnitude
            current_location = next_location
            candidate = current_location + step
            if candidate == current_location ||
               candidate <= left_boundary ||
               candidate >= right_boundary
                break
            end
            next_location = candidate
            next_magnitude = magnitude(next_location)
        end

        extrema[index] = current_location
        magnitudes[index] = current_magnitude
        closest = min(closest, current_magnitude)
        farthest = max(farthest, current_magnitude)
        left_boundary = right_boundary
    end

    spread = farthest - closest
    iszero(closest) || (spread /= closest)
    spread >= previous_spread && (relaxation /= 2)

    for index in 1:equation_count
        ratio = iszero(magnitudes[index + 1]) ? BigFloat(0.0625) :
                magnitudes[index] / magnitudes[index + 1] - one(BigFloat)
        ratio = min(ratio, BigFloat(0.25))
        steps[index] = ratio * (extrema[index + 1] - extrema[index]) * relaxation
    end
    steps[end] = steps[end - 1]

    for index in 1:equation_count
        candidate = zeros[index] - steps[index]
        candidate <= lower && continue
        candidate >= upper && continue
        candidate <= extrema[index] &&
            (candidate = (extrema[index] + zeros[index]) / 2)
        candidate >= extrema[index + 1] &&
            (candidate = (extrema[index + 1] + zeros[index]) / 2)
        zeros[index] = candidate
    end
    return relaxation, spread, magnitudes
end

function _legacy_remez(
    y::Integer,
    z::Integer,
    n::Integer,
    d::Integer,
    lower::Union{Real,AbstractString},
    upper::Union{Real,AbstractString};
    precision::Integer=256,
    tolerance::Union{Nothing,Real,AbstractString}=nothing,
    max_iterations::Integer=10_000,
    exponential_coefficients=nothing,
    exponential_powers=nothing,
)
    y >= 0 || throw(ArgumentError("the exponent numerator must be nonnegative"))
    z > 0 || throw(ArgumentError("the exponent denominator must be positive"))
    n >= 0 || throw(ArgumentError("the numerator degree must be nonnegative"))
    d >= 0 || throw(ArgumentError("the denominator degree must be nonnegative"))
    precision >= 64 || throw(ArgumentError("precision must be at least 64 bits"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))

    return setprecision(BigFloat, Int(precision)) do
        # The C++ constructor and command line both receive IEEE Float64 bounds.
        lower_big = BigFloat(Float64(lower isa AbstractString ? parse(Float64, lower) : lower))
        upper_big = BigFloat(Float64(upper isa AbstractString ? parse(Float64, upper) : upper))
        0 < lower_big < upper_big || throw(ArgumentError("expected 0 < lower < upper"))

        coefficients_big = if exponential_coefficients === nothing
            BigFloat[]
        else
            BigFloat[
                BigFloat(
                    coefficient isa AbstractString ? parse(Float64, coefficient) :
                    Float64(coefficient),
                ) for coefficient in exponential_coefficients
            ]
        end
        powers_int = exponential_powers === nothing ? Int[] :
                     Int.(collect(exponential_powers))
        length(coefficients_big) == length(powers_int) || throw(
            DimensionMismatch(
                "exponential_coefficients and exponential_powers must have equal length",
            ),
        )
        length(coefficients_big) <= 10 ||
            throw(ArgumentError("the legacy algorithm accepts at most 10 exponential terms"))
        all(isfinite, coefficients_big) ||
            throw(ArgumentError("exponential coefficients must be finite"))

        target_tolerance = tolerance === nothing ? BigFloat(Float64(1.0e-15)) :
                           _as_bigfloat(tolerance)
        0 < target_tolerance < 1 ||
            throw(ArgumentError("tolerance must be between zero and one"))

        numerator_degree = Int(n)
        denominator_degree = Int(d)
        power_num = Int(y)
        power_den = Int(z)
        equation_count = numerator_degree + denominator_degree + 1
        zeros, extrema = _legacy_initial_points(lower_big, upper_big, equation_count)
        steps = Vector{BigFloat}(undef, equation_count + 1)
        steps[1] = zeros[1] - lower_big
        for index in 2:equation_count
            steps[index] = zeros[index] - zeros[index - 1]
        end
        steps[end] = steps[end - 1]

        relaxation = BigFloat(0.25)
        spread = BigFloat(Float64(1.0e37))
        parameters = BigFloat[]
        solved_zeros = copy(zeros[1:equation_count])
        magnitudes = fill(BigFloat(Inf), equation_count + 1)
        iterations = 0

        while spread > target_tolerance
            iterations += 1
            iterations <= max_iterations ||
                throw(RemezConvergenceError(Int(max_iterations), spread))
            solved_zeros = copy(zeros[1:equation_count])
            parameters = _legacy_equations(
                solved_zeros,
                numerator_degree,
                denominator_degree,
                power_num,
                power_den,
                coefficients_big,
                powers_int,
            )
            relaxation < target_tolerance &&
                throw(RemezConvergenceError(iterations, spread))
            relaxation, spread, magnitudes = _legacy_search!(
                zeros,
                extrema,
                steps,
                parameters,
                numerator_degree,
                denominator_degree,
                lower_big,
                upper_big,
                power_num,
                power_den,
                coefficients_big,
                powers_int,
                relaxation,
                spread,
            )
        end

        numerator = copy(parameters[1:(numerator_degree + 1)])
        denominator = vcat(
            parameters[(numerator_degree + 2):end],
            one(BigFloat),
        )
        extrema_errors = BigFloat[
            _legacy_signed_error(
                x,
                parameters,
                numerator_degree,
                denominator_degree,
                power_num,
                power_den,
                coefficients_big,
                powers_int,
            ) for x in extrema
        ]
        return RationalApproximation(
            numerator,
            denominator,
            power_num,
            power_den,
            copy(coefficients_big),
            copy(powers_int),
            lower_big,
            upper_big,
            magnitudes[1],
            iterations,
            solved_zeros,
            copy(extrema),
            extrema_errors,
            Int(precision),
        )
    end
end

"""
    remez(y, z, n, d, lower, upper; algorithm=:legacy, kwargs...)

Construct a relative-error minimax rational approximation. The default
`:legacy` algorithm follows the original C++ exchange calculation, including
its initial-point and extrema-search rules. Use `algorithm=:refined` for the
independent golden-section implementation. `precision` is measured in bits.
"""
function remez(
    y::Integer,
    z::Integer,
    n::Integer,
    d::Integer,
    lower::Union{Real,AbstractString},
    upper::Union{Real,AbstractString};
    algorithm::Symbol=:legacy,
    precision::Integer=256,
    tolerance::Union{Nothing,Real,AbstractString}=nothing,
    max_iterations::Integer=10_000,
    extrema_iterations::Integer=160,
    exponential_coefficients=nothing,
    exponential_powers=nothing,
)
    if algorithm === :legacy
        return _legacy_remez(
            y,
            z,
            n,
            d,
            lower,
            upper;
            precision,
            tolerance,
            max_iterations,
            exponential_coefficients,
            exponential_powers,
        )
    elseif algorithm === :refined
        return _refined_remez(
            y,
            z,
            n,
            d,
            lower,
            upper;
            precision,
            tolerance,
            max_iterations,
            extrema_iterations,
            exponential_coefficients,
            exponential_powers,
        )
    end
    throw(ArgumentError("algorithm must be :legacy or :refined"))
end

remez(y::Integer, z::Integer, degree::Integer, lower, upper; kwargs...) =
    remez(y, z, degree, degree, lower, upper; kwargs...)

function _polynomial_derivative_coefficients(coefficients::Vector{BigFloat})
    degree = length(coefficients) - 1
    degree == 0 && return BigFloat[]
    return BigFloat[i * coefficients[i + 1] for i in 1:degree]
end

function _bisect_polynomial_root(
    coefficients::Vector{BigFloat},
    left::BigFloat,
    right::BigFloat,
    tolerance::BigFloat,
    maximum_iterations::Int,
)
    left_value = _evaluate_polynomial(left, coefficients)
    right_value = _evaluate_polynomial(right, coefficients)
    iszero(left_value) && return left
    iszero(right_value) && return right
    signbit(left_value) != signbit(right_value) ||
        throw(ArgumentError("a real polynomial root was not bracketed"))

    for _ in 1:maximum_iterations
        midpoint = (left + right) / 2
        midpoint_value = _evaluate_polynomial(midpoint, coefficients)
        iszero(midpoint_value) && return midpoint
        if signbit(left_value) == signbit(midpoint_value)
            left = midpoint
            left_value = midpoint_value
        else
            right = midpoint
            right_value = midpoint_value
        end
        right - left <= tolerance * max(one(BigFloat), abs(midpoint)) &&
            return (left + right) / 2
    end
    throw(ArgumentError("a real polynomial root did not converge"))
end

"""
Find all distinct real roots using derivative-root isolation and BigFloat
bisection. If a polynomial has degree `m`, its derivative roots divide the real
axis into intervals containing at most one root. A Cauchy bound makes those
intervals finite. Requiring exactly `m` sign-changing intervals both detects
complex roots and rejects repeated roots, which cannot have a simple PFE.
"""
function _polynomial_roots(coefficients::Vector{BigFloat}, precision::Int)
    degree = length(coefficients) - 1
    degree == 0 && return BigFloat[]
    iszero(coefficients[end]) && throw(ArgumentError("the polynomial degree is defective"))
    degree == 1 && return BigFloat[-coefficients[1] / coefficients[2]]

    derivative = _polynomial_derivative_coefficients(coefficients)
    critical_points = _polynomial_roots(derivative, precision)
    leading = abs(coefficients[end])
    bound = 2 * (one(BigFloat) + maximum(abs(coefficient) / leading for coefficient in coefficients[1:end-1]))
    boundaries = vcat(-bound, critical_points, bound)
    tolerance = BigFloat(2)^(-precision + 16)
    maximum_iterations = max(256, 4precision)
    roots = BigFloat[]

    for i in 1:(length(boundaries) - 1)
        left = boundaries[i]
        right = boundaries[i + 1]
        left_value = _evaluate_polynomial(left, coefficients)
        right_value = _evaluate_polynomial(right, coefficients)
        if iszero(left_value) || iszero(right_value)
            throw(ArgumentError("partial fractions require distinct simple real roots"))
        end
        if signbit(left_value) != signbit(right_value)
            push!(
                roots,
                _bisect_polynomial_root(
                    coefficients,
                    left,
                    right,
                    tolerance,
                    maximum_iterations,
                ),
            )
        end
    end

    length(roots) == degree ||
        throw(ArgumentError("partial fractions require distinct real roots"))
    return roots
end

function _polynomial_derivative_value(x::BigFloat, coefficients::Vector{BigFloat})
    degree = length(coefficients) - 1
    degree == 0 && return zero(BigFloat)
    value = BigFloat(degree) * coefficients[end]
    for exponent in (degree - 1):-1:1
        value = muladd(value, x, BigFloat(exponent) * coefficients[exponent + 1])
    end
    return value
end

"""
    _refined_partial_fraction(approximation; inverse=false)

Convert an equal-degree rational approximation to
`constant + sum(residue / (x + shift))`. With `inverse=true`, expand the
reciprocal using derivative-root isolation.
"""
function _refined_partial_fraction(
    approximation::RationalApproximation{BigFloat};
    inverse::Bool=false,
)
    numerator_degree = length(approximation.numerator) - 1
    denominator_degree = length(approximation.denominator) - 1
    numerator_degree == denominator_degree ||
        throw(ArgumentError("partial fractions currently require equal degrees"))

    return setprecision(BigFloat, approximation.precision) do
        numerator = inverse ? approximation.denominator : approximation.numerator
        denominator = inverse ? approximation.numerator : approximation.denominator
        poles = _polynomial_roots(denominator, approximation.precision)
        constant = numerator[end] / denominator[end]
        residues = BigFloat[]
        for pole in poles
            derivative = _polynomial_derivative_value(pole, denominator)
            iszero(derivative) && throw(ArgumentError("partial fractions require simple poles"))
            push!(residues, _evaluate_polynomial(pole, numerator) / derivative)
        end
        shifts = .-poles
        order = sortperm(shifts)
        return PartialFractionApproximation(constant, residues[order], shifts[order])
    end
end

function _legacy_polynomial_value(
    x::BigFloat,
    coefficients::Vector{BigFloat},
    degree::Int,
)
    value = coefficients[degree + 1]
    for index in degree:-1:1
        value = value * x + coefficients[index]
    end
    return value
end

function _legacy_polynomial_derivative_value(
    x::BigFloat,
    coefficients::Vector{BigFloat},
    degree::Int,
)
    value = BigFloat(degree) * coefficients[degree + 1]
    for index in (degree - 1):-1:1
        value = value * x + BigFloat(index) * coefficients[index + 1]
    end
    return value
end

function _legacy_newton_root(coefficients::Vector{BigFloat}, degree::Int)
    root = BigFloat(0.5) * (BigFloat(-100_000) + one(BigFloat))
    tolerance = BigFloat(Float64(1.0e-20))
    for _ in 1:10_000
        value = _legacy_polynomial_value(root, coefficients, degree)
        derivative = _legacy_polynomial_derivative_value(root, coefficients, degree)
        step = value / derivative
        root -= step
        abs(step) < tolerance && return root
    end
    throw(ArgumentError("legacy Newton root finder did not converge"))
end

function _legacy_polynomial_roots(coefficients::Vector{BigFloat})
    degree = length(coefficients) - 1
    iszero(degree) && return BigFloat[]
    polynomial = copy(coefficients)
    roots = Vector{BigFloat}(undef, degree)
    for current_degree in degree:-1:1
        root = _legacy_newton_root(polynomial, current_degree)
        iszero(root) && throw(ArgumentError("legacy Newton root finder returned zero"))
        roots[current_degree] = root
        polynomial[1] = -polynomial[1] / root
        for index in 2:current_degree
            polynomial[index] = (polynomial[index - 1] - polynomial[index]) / root
        end
    end
    return roots
end

function _legacy_pfe(
    numerator_roots::Vector{BigFloat},
    denominator_roots::Vector{BigFloat},
    normalization::BigFloat,
)
    degree = length(numerator_roots)
    length(denominator_roots) == degree ||
        throw(DimensionMismatch("legacy partial fractions require equal degrees"))
    iszero(degree) &&
        return PartialFractionApproximation(normalization, BigFloat[], BigFloat[])

    numerator = zeros(BigFloat, degree)
    denominator = zeros(BigFloat, degree)
    numerator[1] = one(BigFloat)
    denominator[1] = one(BigFloat)

    for root_index in 1:degree
        for coefficient_index in degree:-1:1
            numerator[coefficient_index] *= -numerator_roots[root_index]
            denominator[coefficient_index] *= -denominator_roots[root_index]
            if coefficient_index > 1
                numerator[coefficient_index] += numerator[coefficient_index - 1]
                denominator[coefficient_index] += denominator[coefficient_index - 1]
            end
        end
    end
    numerator .-= denominator

    residues = Vector{BigFloat}(undef, degree)
    for pole_index in 1:degree
        residue = zero(BigFloat)
        for coefficient_index in degree:-1:1
            residue = denominator_roots[pole_index] * residue +
                      numerator[coefficient_index]
        end
        for other_index in degree:-1:1
            pole_index == other_index && continue
            residue /= denominator_roots[pole_index] - denominator_roots[other_index]
        end
        residues[pole_index] = residue * normalization
    end

    shifts = .-denominator_roots
    for index in 1:degree
        smallest = index
        for candidate in (index + 1):degree
            shifts[candidate] < shifts[smallest] && (smallest = candidate)
        end
        if smallest != index
            shifts[index], shifts[smallest] = shifts[smallest], shifts[index]
            residues[index], residues[smallest] = residues[smallest], residues[index]
        end
    end
    return PartialFractionApproximation(normalization, residues, shifts)
end

function _legacy_partial_fractions(approximation::RationalApproximation{BigFloat})
    numerator_degree = length(approximation.numerator) - 1
    denominator_degree = length(approximation.denominator) - 1
    numerator_degree == denominator_degree ||
        throw(ArgumentError("legacy partial fractions require equal degrees"))

    return setprecision(BigFloat, approximation.precision) do
        if iszero(numerator_degree)
            normalization = approximation.numerator[1] / approximation.denominator[1]
            positive = PartialFractionApproximation(
                normalization,
                BigFloat[],
                BigFloat[],
            )
            negative = PartialFractionApproximation(
                inv(normalization),
                BigFloat[],
                BigFloat[],
            )
            return positive, negative
        end

        numerator_roots = _legacy_polynomial_roots(approximation.numerator)
        denominator_roots = _legacy_polynomial_roots(approximation.denominator)
        normalization = approximation.numerator[end]
        positive = _legacy_pfe(numerator_roots, denominator_roots, normalization)
        negative = _legacy_pfe(
            denominator_roots,
            numerator_roots,
            inv(normalization),
        )
        return positive, negative
    end
end

"""
    partial_fraction(approximation; inverse=false, algorithm=:legacy)

Convert an equal-degree rational approximation to partial-fraction form. The
default reproduces the C++ Newton-deflation and residue calculation. Select
`:refined` for derivative-root isolation.
"""
function partial_fraction(
    approximation::RationalApproximation{BigFloat};
    inverse::Bool=false,
    algorithm::Symbol=:legacy,
)
    if algorithm === :legacy
        positive, negative = _legacy_partial_fractions(approximation)
        return inverse ? negative : positive
    elseif algorithm === :refined
        return _refined_partial_fraction(approximation; inverse)
    end
    throw(ArgumentError("algorithm must be :legacy or :refined"))
end

function _legacy_precision_bits(decimal_precision::Integer)
    decimal_precision > 0 || throw(ArgumentError("precision must be positive"))
    requested_bits = floor(Int, 3.321928094 * decimal_precision)
    limb_bits = Sys.WORD_SIZE
    return max(limb_bits, cld(requested_bits, limb_bits) * limb_bits)
end

"""
    calc_coefficients(y, z, n, lower, upper; precision=42, kwargs...)

Drop-in replacement for the historical Julia wrapper around `AlgRemez_jll`.
It computes an `(n,n)` approximation natively and returns Float64
`(coeff_plus, coeff_minus)` partial-fraction coefficients for `x^(y/z)` and
its reciprocal. Here `precision` is measured in decimal digits for compatibility;
the lower-level [`remez`](@ref) API measures precision in bits. The requested
precision is rounded to a GMP limb boundary exactly as in the original C++
implementation (42 decimal digits therefore use 192 bits on a 64-bit host).
"""
function calc_coefficients(
    y::Integer,
    z::Integer,
    n::Integer,
    lower::Union{Real,AbstractString},
    upper::Union{Real,AbstractString};
    precision::Integer=42,
    tolerance::Union{Real,AbstractString}=1.0e-15,
    kwargs...,
)
    precision > 0 || throw(ArgumentError("precision must be positive"))
    bits = _legacy_precision_bits(precision)
    approximation = remez(
        y,
        z,
        n,
        n,
        lower,
        upper;
        precision=bits,
        tolerance=tolerance,
        algorithm=:legacy,
        kwargs...,
    )
    positive_pfe, negative_pfe = _legacy_partial_fractions(approximation)
    positive = AlgRemez_coeffs(positive_pfe)
    negative = AlgRemez_coeffs(negative_pfe)
    return positive, negative
end

"""
    write_legacy_files(approximation; directory=pwd(), ratio=1.01)

Write `approx.dat` and `error.dat` in the format produced by the original C++
executable. The coefficient file contains the partial-fraction expansions of
both the approximation and its reciprocal. As in the original implementation,
this operation requires equal numerator and denominator degrees.
"""
function write_legacy_files(
    approximation::RationalApproximation{BigFloat};
    directory::AbstractString=pwd(),
    ratio::Real=1.01,
)
    ratio > 1 || throw(ArgumentError("ratio must be greater than one"))
    numerator_degree = length(approximation.numerator) - 1
    denominator_degree = length(approximation.denominator) - 1
    numerator_degree == denominator_degree ||
        throw(ArgumentError("legacy partial-fraction output requires equal degrees"))

    mkpath(directory)
    approximation_path = joinpath(directory, "approx.dat")
    error_path = joinpath(directory, "error.dat")
    positive, negative = _legacy_partial_fractions(approximation)

    open(approximation_path, "w") do io
        @printf(
            io,
            "Approximation to f(x) = x^(%d/%d)\n\n",
            approximation.power_num,
            approximation.power_den,
        )
        @printf(io, "alpha[0] = %18.16e\n", Float64(positive.constant))
        for i in 1:numerator_degree
            @printf(
                io,
                "alpha[%d] = %18.16e, beta[%d] = %18.16e\n",
                i,
                Float64(positive.residues[i]),
                i,
                Float64(positive.shifts[i]),
            )
        end

        @printf(
            io,
            "\nApproximation to f(x) = x^(-%d/%d)\n\n",
            approximation.power_num,
            approximation.power_den,
        )
        @printf(io, "alpha[0] = %18.16e\n", Float64(negative.constant))
        for i in 1:numerator_degree
            @printf(
                io,
                "alpha[%d] = %18.16e, beta[%d] = %18.16e\n",
                i,
                Float64(negative.residues[i]),
                i,
                Float64(negative.shifts[i]),
            )
        end
    end

    setprecision(BigFloat, approximation.precision) do
        x = Float64(approximation.lower)
        upper = Float64(approximation.upper)
        step_ratio = Float64(ratio)
        parameters = vcat(approximation.numerator, approximation.denominator[1:end-1])
        open(error_path, "w") do io
            while x < upper
                x_big = BigFloat(x)
                target = Float64(
                    _legacy_target(
                        x_big,
                        approximation.power_num,
                        approximation.power_den,
                        approximation.exponential_coefficients,
                        approximation.exponential_powers,
                    ),
                )
                rational = Float64(
                    _legacy_approximation_value(
                        x_big,
                        parameters,
                        numerator_degree,
                        denominator_degree,
                    ),
                )
                @printf(
                    io,
                    "%e %e\n",
                    x,
                    (rational - target) / target,
                )
                next_x = x * step_ratio
                next_x > x || throw(ArgumentError("ratio does not advance the error grid"))
                x = next_x
            end
        end
    end

    return (approximation=approximation_path, error=error_path)
end

"""
    legacy_main(args=ARGS; directory=pwd())

Command-line-compatible entry point. `args` must be
`y z n d lower upper decimal_precision`. It runs the native Julia algorithm
and writes the two legacy data files.
"""
function legacy_main(args=ARGS; directory::AbstractString=pwd())
    length(args) == 7 || throw(
        ArgumentError(
            "usage: algremez y z numerator_degree denominator_degree lower upper decimal_precision",
        ),
    )
    y = parse(Int, args[1])
    z = parse(Int, args[2])
    n = parse(Int, args[3])
    d = parse(Int, args[4])
    decimal_precision = parse(Int, args[7])
    decimal_precision > 0 || throw(ArgumentError("decimal precision must be positive"))
    bits = _legacy_precision_bits(decimal_precision)

    println("Approximation bounds are [$(args[5]),$(args[6])]")
    println("Precision of arithmetic is $decimal_precision")
    println("Degree of the approximation is ($n,$d)")
    println("Approximating the function x^($y/$z)")
    approximation = remez(
        y,
        z,
        n,
        d,
        args[5],
        args[6];
        precision=bits,
        tolerance=1.0e-15,
        algorithm=:legacy,
    )
    @printf(
        "Converged at %d iterations, error = %e\n",
        approximation.iterations,
        Float64(approximation.maximum_relative_error),
    )
    write_legacy_files(approximation; directory=directory)
    return approximation
end

end
