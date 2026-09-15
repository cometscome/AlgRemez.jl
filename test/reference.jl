module AlgRemezReference

using AlgRemez_jll

export ErrorSample,
       PartialFraction,
       ReferenceResult,
       jll_reference,
       log_grid,
       sign_change_count,
       ulp_distance

const _FLOAT_PATTERN = raw"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
const _CONSTANT_PATTERN = Regex(raw"alpha\[0\]\s*=\s*(" * _FLOAT_PATTERN * ")")
const _TERM_PATTERN = Regex(
    raw"alpha\[(\d+)\]\s*=\s*(" * _FLOAT_PATTERN *
    raw"),\s*beta\[(\d+)\]\s*=\s*(" * _FLOAT_PATTERN * ")",
)
const _CONVERGENCE_PATTERN = Regex(
    raw"Converged at\s+(\d+)\s+iterations,\s+error\s*=\s*(" *
    _FLOAT_PATTERN * ")",
)

"""Partial-fraction representation `constant + sum(residues ./ (x .+ shifts))`."""
struct PartialFraction{T<:AbstractFloat}
    constant::T
    residues::Vector{T}
    shifts::Vector{T}

    function PartialFraction(constant::T, residues::Vector{T}, shifts::Vector{T}) where {T<:AbstractFloat}
        length(residues) == length(shifts) ||
            throw(DimensionMismatch("residue and shift counts differ"))
        return new{T}(constant, residues, shifts)
    end
end

Base.length(p::PartialFraction) = length(p.residues)

function (p::PartialFraction{T})(x::Number) where {T}
    S = promote_type(T, typeof(x))
    value = convert(S, p.constant)
    xx = convert(S, x)
    for i in eachindex(p.residues, p.shifts)
        value += convert(S, p.residues[i]) / (xx + convert(S, p.shifts[i]))
    end
    return value
end

function _ordered_float_bits(x::Float64)
    bits = reinterpret(UInt64, x)
    sign_mask = UInt64(1) << 63
    return (bits & sign_mask) == 0 ? bits | sign_mask : ~bits
end

"""Return the number of representable Float64 values between two values."""
function ulp_distance(x::Float64, y::Float64)
    isfinite(x) && isfinite(y) || throw(ArgumentError("ULP inputs must be finite"))
    x_bits = _ordered_float_bits(x)
    y_bits = _ordered_float_bits(y)
    return x_bits >= y_bits ? x_bits - y_bits : y_bits - x_bits
end

struct ErrorSample{T<:AbstractFloat}
    x::T
    error::T
end

"""All externally observable results produced by one `AlgRemez_jll` run."""
struct ReferenceResult{T<:AbstractFloat}
    positive::PartialFraction{T}
    negative::PartialFraction{T}
    error_samples::Vector{ErrorSample{T}}
    iterations::Int
    reported_error::T
    log::String
end

function _parse_partial_fraction(section::AbstractString, degree::Integer)
    constant_match = match(_CONSTANT_PATTERN, section)
    constant_match === nothing && error("alpha[0] was not found in approx.dat")

    residues = fill(NaN, degree)
    shifts = fill(NaN, degree)
    seen = falses(degree)
    for term_match in eachmatch(_TERM_PATTERN, section)
        alpha_index = parse(Int, term_match.captures[1])
        beta_index = parse(Int, term_match.captures[3])
        alpha_index == beta_index ||
            error("alpha[$alpha_index] is paired with beta[$beta_index]")
        1 <= alpha_index <= degree ||
            error("coefficient index $alpha_index is outside 1:$degree")
        seen[alpha_index] && error("coefficient index $alpha_index occurs twice")
        residues[alpha_index] = parse(Float64, term_match.captures[2])
        shifts[alpha_index] = parse(Float64, term_match.captures[4])
        seen[alpha_index] = true
    end
    all(seen) || error("approx.dat contains $(count(seen)) of $degree expected terms")

    constant = parse(Float64, constant_match.captures[1])
    return PartialFraction(constant, residues, shifts)
end

function _parse_approximation(text::AbstractString, degree::Integer)
    sections = split(text, r"Approximation to f\(x\) =")
    length(sections) == 3 ||
        error("expected positive- and negative-power sections in approx.dat")
    positive = _parse_partial_fraction(sections[2], degree)
    negative = _parse_partial_fraction(sections[3], degree)
    return positive, negative
end

function _parse_error_samples(text::AbstractString)
    samples = ErrorSample{Float64}[]
    for (line_number, line) in enumerate(eachline(IOBuffer(text)))
        isempty(strip(line)) && continue
        fields = split(line)
        length(fields) == 2 ||
            error("error.dat line $line_number does not contain two fields")
        push!(samples, ErrorSample(parse(Float64, fields[1]), parse(Float64, fields[2])))
    end
    isempty(samples) && error("error.dat is empty")
    issorted(sample.x for sample in samples) || error("error.dat abscissae are not sorted")
    return samples
end

function _parse_convergence(log::AbstractString)
    convergence_match = match(_CONVERGENCE_PATTERN, log)
    convergence_match === nothing && error("convergence record was not found in stdout")
    iterations = parse(Int, convergence_match.captures[1])
    reported_error = parse(Float64, convergence_match.captures[2])
    return iterations, reported_error
end

"""
    jll_reference(y, z, n, d, lower, upper; precision=42)

Run the packaged legacy C++ executable in an isolated temporary directory and
return its partial fractions, sampled relative errors, and convergence record.

The legacy partial-fraction routine is only safe for equal positive degrees.
In particular, the executable has an out-of-bounds access for `(n, d) = (0, 0)`.
"""
function jll_reference(
    y::Integer,
    z::Integer,
    n::Integer,
    d::Integer,
    lower::Real,
    upper::Real;
    precision::Integer=42,
)
    y > 0 || throw(ArgumentError("the exponent numerator must be positive"))
    z > 0 || throw(ArgumentError("the exponent denominator must be positive"))
    n == d || throw(ArgumentError("the legacy partial-fraction output requires n == d"))
    n > 0 || throw(ArgumentError("the legacy executable is unsafe for zero degree"))
    0 < lower < upper || throw(ArgumentError("expected 0 < lower < upper"))
    precision > 0 || throw(ArgumentError("precision must be positive"))

    return mktempdir(prefix="algremez-jll-") do directory
        output = IOBuffer()
        executable = AlgRemez_jll.algremez()
        command = `$executable $y $z $n $d $lower $upper $precision`
        run(pipeline(Cmd(command; dir=directory); stdout=output, stderr=output))
        log = String(take!(output))

        approx_path = joinpath(directory, "approx.dat")
        error_path = joinpath(directory, "error.dat")
        isfile(approx_path) || error("legacy executable did not create approx.dat\n$log")
        isfile(error_path) || error("legacy executable did not create error.dat\n$log")

        positive, negative = _parse_approximation(read(approx_path, String), n)
        samples = _parse_error_samples(read(error_path, String))
        iterations, reported_error = _parse_convergence(log)
        return ReferenceResult(positive, negative, samples, iterations, reported_error, log)
    end
end

"""Return geometrically spaced `Float64` points, including both endpoints."""
function log_grid(lower::Real, upper::Real, length::Integer=257)
    0 < lower < upper || throw(ArgumentError("expected 0 < lower < upper"))
    length >= 2 || throw(ArgumentError("the grid needs at least two points"))
    return exp.(range(log(Float64(lower)), log(Float64(upper)); length=length))
end

"""Count sign changes after ignoring samples that are exactly zero."""
function sign_change_count(samples::AbstractVector{<:ErrorSample})
    signs = [sign(sample.error) for sample in samples if !iszero(sample.error)]
    return count(pair -> pair[1] != pair[2], zip(signs, Iterators.drop(signs, 1)))
end

end
