using NumericalIntegration: integrate
using DelimitedFiles: writedlm, readdlm

##################################################################
# -- convolve raw data with a function (gaussian, lorentzian, etc)
##################################################################

PROGRAM = basename(@__FILE__)

ln(x) = log(x)

logrange(a,b,n) = (10^x for x in range(log10(a), log10(b), length=n))

function usage(::Int)
  println()
  println(PROGRAM,": Convolves data from a datafile with a function ")
  println()
  println("    usage: ",PROGRAM," [FILE]")
  println()
  return true
end

function gauss(x,x0,γ)
  @. return exp(- (x-x0)^2 / (2*γ^2))
end

function cauchy(x,x0,γ)
  @. return γ / π  / ((x-x0)^2 + γ^2)
end

function convolve(x,x0,γ,conv_type::String)
  if(conv_type == "gauss")
    return gauss(x,x0,γ)
  elseif(conv_type == "cauchy")
    return cauchy(x,x0,γ)
  end
end

function main()

  narg=6

  if(size(ARGS,1) != narg)
    println("STOP: ",PROGRAM," requires ",narg," arguments, but ",string(size(ARGS,1))," were supplied.")
    usage(1)
    exit()
  end

  datafile = ARGS[1]
  outfile  = ARGS[2]
  γ        = parse(Float64,ARGS[3]) # <-- width in x-units
  nx       = parse(Int,ARGS[4]) # <------ number of convolution grid points
  logx     = parse(Bool,ARGS[5]) # <----- T: logarithmic convolution grid
                                 #        F: linear      convolution grid
  conv_type = ARGS[6] # <---------------- type of convolution function
                     #                    - gauss
                     #                    - cauchy

  # -- check if input datafile exists
  if( ! isfile(datafile) )
    usage(1)
    println("STOP: Need a data file to convole")
    exit()
  end

  # -- assumes x y format, one y column (for now)
  data=readdlm(datafile)
  xdata = data[:,1]
  # ny = size(data,2) - 1
  # if(ny < 1)
  #   println("STOP: ",outfile," does not have y values to convolve")
  #   usage(1)
  #   exit()
  # end
  # ydata = data[:,2:ny+1]
  ydata = data[:,2]

  # -- sort data based on x
  indices = sortperm(xdata)
  xdata = xdata[indices]
  ydata = ydata[indices]

  if(logx)
    xconv = [ x for x in logrange(xdata[firstindex(xdata)],xdata[lastindex(xdata)],nx) ]
  else
    xconv = [ x for x in range(xdata[firstindex(xdata)],xdata[lastindex(xdata)],nx) ]
  end

  yconv = [ integrate(xdata,ydata.*convolve(xdata,x,γ,conv_type)) for x in xconv ] ./ [ integrate(xdata,convolve(xdata,x,γ,conv_type)) for x in xconv ]

  writedlm(outfile, [xconv yconv])

end

main()
