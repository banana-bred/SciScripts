## Generate Γ expansion coefficients for a given numerical precision

# -- set precision
setprec(bits :: Integer)= (setprecision(BigFloat, bits); nothing)

# -- Spouge error bound
spouge_bound(a :: Integer) = (2BigFloat(pi))^(-a-BigFloat(0.5)) / sqrt(a)

"""
    spouge_coeffs(a :: Int, bits = 256)

Compute Spouge coefficients c0..c(a-1) at given precision (bits)
c0 = sqrt(2π)
ck = (-1)^{k-1}/(k-1)! * exp(a-k)*(a-k)^(k-1/2), k=1..a-1
"""
function spouge_coeffs(a :: Int, bits :: Int=256)
  @assert a ≥ 3
  setprec(bits)
  c = Vector{BigFloat}(undef, a)
  c[begin + 0]= sqrt(2BigFloat(pi)) # c0
  for k in 1:a-1
    kk = BigFloat(k)
    ak = BigFloat(a) - kk
    term = (-one(BigFloat))^(k-1) / factorial(big(k-1))
    term *= exp(ak)
    term *= ak^(kk - BigFloat(0.5))
    c[begin + k] = term
  end
  return c
end

if abspath(PROGRAM_FILE) == @__FILE__
  isempty(ARGS) && exit()
  nargs = length(ARGS)
  a         = parse(Int,  ARGS[1])
  bits      = parse(Int,  ARGS[2])
  logCoeffs = parse(Bool, ARGS[3])
  # bits = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 256
  cs = spouge_coeffs(a, bits)
  if(logCoeffs == true) ; cs = log.(cs); end
  for ci in cs
    println(ci)
  end
end
