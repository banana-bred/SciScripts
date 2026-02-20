#!/usr/bin/env julia
# Generate coefficients in the Lanczos approximation
#  Γ(z+1) = sqrt(2π) (z+g+1/2)^(z+1/2)  exp(-z-g-1/2) * A(z)

using LinearAlgebra
using Random
using SpecialFunctions

function main()
end # main

g = 27.0
N = 32
NDIGITS   = 40
PREC_BITS = 256

if abspath(PROGRAM_FILE) == @__FILE__
  isempty(ARGS) && exit()
  nargs = length(ARGS)
  NDIGITS = parse(Int, ARGS[1])
  PREC_BITS = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 256
  cs = spouge_coeffs(a, bits)
  for ci in cs
    println(ci)
  end
end

main()
