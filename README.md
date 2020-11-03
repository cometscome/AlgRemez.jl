# AlgRemez.jl

[![Build Status](https://travis-ci.org/cometscome/AlgRemez.jl.svg?branch=master)](https://travis-ci.com/cometscome/AlgRemez.jl)
[![Coverage](https://codecov.io/gh/cometscome/AlgRemez.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cometscome/AlgRemez.jl)
[![Coverage](https://coveralls.io/repos/github/cometscome/AlgRemez.jl/badge.svg?branch=master)](https://coveralls.io/github/cometscome/AlgRemez.jl?branch=master) 


AlgRemez.jl is a Julia wrapper for AlgRemez written in c++.  

    Please see,

    https://github.com/maddyscientist/AlgRemez

    
>The archive downloadable here contains an implementation of the Remez algorithm which calculates optimal rational (and polynomial) approximations to the nth root over a given spectral range.  The Remez algorithm, although in principle is extremely straightforward to
program, is quite difficult to get completely correct, e.g., the Maple implementation of the algorithm does not always converge to the 
correct answer.

# How to use
## Install: 
Add ] and install as follows
```
add https://github.com/cometscome/AlgRemez.jl
```

## how to use

There are only two functions, calc_coefficients and fittedfunction. 
You can get coefficients as follows.

```julia
y = 1
z = 2
n = 5
lambda_low = 0.0004
lambda_high = 64
coeff_plus,coeff_minus = AlgRemez.calc_coefficients(y,z,n,lambda_low,lambda_high)
println(coeff_plus.α0)
println(coeff_plus.α)
println(coeff_plus.β)

println(coeff_minus.α0)
println(coeff_minus.α)
println(coeff_minus.β)

funcplus = AlgRemez.fittedfunction(coeff_plus)
println("1: sqrt(x)")
println("Approximated value at x=0.1: ",funcplus(0.1))
println("Exact value at x=0.1: ", sqrt(0.1))

funcminus = AlgRemez.fittedfunction(coeff_minus)
println("2: 1/sqrt(x)")
println("Approximated value at x=0.1: ",funcminus(0.1))
println("Exact value at x=0.1: ", 1/sqrt(0.1))

y = 1
z = 4
n = 5
lambda_low = 0.0004
lambda_high = 64
println("The order is $n")
coeff_plus,coeff_minus = AlgRemez.calc_coefficients(y,z,n,lambda_low,lambda_high)
println(coeff_plus.α)
println(coeff_plus.β)


funcplus = AlgRemez.fittedfunction(coeff_plus)
println("3: x^($y/$z)")
println("Approximated value at x=0.1: ",funcplus(0.1))
println("Exact value at x=0.1: ", 0.1^(y/z))

funcminus = AlgRemez.fittedfunction(coeff_minus)
println("4: x^(-$y/$z)")
println("Approximated value at x=0.1: ",funcminus(0.1))
println("Exact value at x=0.1: ", 0.1^(-y/z))
```

## two functions

```julia
    """
Output: coeff_plus::AlgRemez_coeffs,coeff_minus::AlgRemez_coeffs
struct AlgRemez_coeffs
    α0::Float64
    α::Array{Float64,1}
    β::Array{Float64,1}
    n::Int64
end
f(x) = x^(y/z) = coeff_plus.α0 + sum_i^n coeff_plus.α[i]/(x + coeff_plus.β[i])
f(x) = x^(-y/z) = coeff_minus.α0 + sum_i^n coeff_minus.α[i]/(x + coeff_minus.β[i])

The function to be approximated f(x) = x^(y/z) and f(x) = x^(-y/z), 
with degree (n,n) over the spectral range [lambda_low,lambda_high],
using precision digits of precision in the arithmetic. 
The parameters y and z must be positive, the approximation to f(x) = x^(-y/z) is simply
the inverse of the approximation to f(x) = x^(y/z).
The default value of precision is 42. 
    """
    function calc_coefficients(y,z,n,lambda_low,lambda_high;precision=42)
```

```julia
fittedfunction(coeff::AlgRemez_coeffs)
```