##################################################################
# perform a B-spline interpolation of a set of points
##################################################################

using LsqFit
using FStrings
using BSplineKit
using DelimitedFiles: writedlm, readdlm

global const PROGRAM = basename(@__FILE__)

##################################################################
######################### FUNCTIONS ##############################
##################################################################

function usage()
  println()
  println(PROGRAM,": Fit a potential energy curve to a Lennard-Jones potential")
  println()
  println("  usage: ", PROGRAM, "[FILE] [OPTIONS]")
  println()
  return true
end # function usage

function die(message :: String)
    println("@@@@@@@@@@@@@@@@@@@@@@@@@")
    println("ERROR: ",message)
    println("@@@@@@@@@@@@@@@@@@@@@@@@@")
    usage()
    exit()
end # function die

function logrange(xi :: Number, xf :: Number, nx :: Integer)
    return (10^x for x in range(log10(xi), log10(xf), length = nx))
end

function main(
        infile       :: String,
        outfile      :: String,
        kx           :: Integer,
        xi           :: Number,
        xf           :: Number,
        nx           :: Integer,
        logx         :: Bool,
        comment_char :: Char
    )

    # -- read data into x and y arrays
    data   = readdlm(infile, comments = true, comment_char = comment_char)
    _, ny = size(data)
    xdata  = data[:,1]
    ydata  = data[:,2:ny]
    ny    -= 1
    data   = 0 # <-- clear data (it shouldn't be that big anyway)

    dx = (xf - xi) / (nx - 1)
    # -- grid for the interpolated data
    xrange = logx ? logrange(xi, xf, nx) : xi:dx:xf;
    xgrid = [x for x in xrange]

    # -- interpolate the y data
    interpolantArray = [ interpolate(xdata, ydata[:, i], BSplineOrder(kx)) for i in 1:ny ]
    yfit = Array{Float64}(undef, nx, ny)
    for j in 1:ny
        yfit[:, j] .= interpolantArray[j].(xgrid)
    end

    # -- write interpolated data to disk
    writedlm(outfile, [xgrid yfit])

end # fuction main

##################################################################
######################### EXECUTION ##############################
##################################################################

const nargs = 8

if(size(ARGS, 1) != nargs)
    println(f"STOP: {PROGRAM} requires {nargs} arguments, but {string(size(ARGS, 1))} were supplied.")
    usage(1)
    exit()
end

infile       = ARGS[1]                 # <-- : name of the input file containing the potential data
outfile      = ARGS[2]                 # <-- : name of the output file
kx           = parse(Int, ARGS[3])     # <-- : the order of the interpolating splines
xi           = parse(Float64, ARGS[4]) # <-- : the smallest value for the grid over which data will be interpolated
xf           = parse(Float64, ARGS[5]) # <-- : the largest value for the grid over which data will be interpolated
nx           = parse(Int, ARGS[6])     # <-- : the number of points in the grid over which data will be interpolated
logx         = parse(Bool, ARGS[7])    # <-- : is the interpolation grid logaritmic ?
comment_char = only(ARGS[8])           # <-- : the comment character for the readdlm operation

main(infile, outfile, kx, xi, xf, nx, logx, comment_char)
