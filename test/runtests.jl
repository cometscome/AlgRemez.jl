using AlgRemez
using Test

@testset "AlgRemez.jl" begin
    # Write your tests here.
    y = 1
    z = 2
    n = 5
    lambda_low = 0.0004
    lambda_high = 64
    println("The order is $n")
    coeff_plus,coeff_minus = AlgRemez.calc_coefficients(y,z,n,lambda_low,lambda_high)
    println(coeff_plus.α)
    println(coeff_plus.β)


    funcplus = AlgRemez.fittedfunction(coeff_plus)
    println("1: sqrt(x)")
    println("Approximated value at x=0.1: ",funcplus(0.1))
    println("Exact value at x=0.1: ", sqrt(0.1))

    funcminus = AlgRemez.fittedfunction(coeff_minus)
    println("2: 1/sqrt(x)")
    println("Approximated value at x=0.1: ",funcminus(0.1))
    println("Exact value at x=0.1: ", 1/sqrt(0.1))
    

    n = 15
    println("The order is $n")
    coeff_plus,coeff_minus = AlgRemez.calc_coefficients(y,z,n,lambda_low,lambda_high)
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
    

    n = 15
    println("The order is $n")
    coeff_plus,coeff_minus = AlgRemez.calc_coefficients(y,z,n,lambda_low,lambda_high)
    funcplus = AlgRemez.fittedfunction(coeff_plus)
    println("3: x^($y/$z)")
    println("Approximated value at x=0.1: ",funcplus(0.1))
    println("Exact value at x=0.1: ", 0.1^(y/z))

    funcminus = AlgRemez.fittedfunction(coeff_minus)
    println("4: x^(-$y/$z)")
    println("Approximated value at x=0.1: ",funcminus(0.1))
    println("Exact value at x=0.1: ", 0.1^(-y/z))
end
