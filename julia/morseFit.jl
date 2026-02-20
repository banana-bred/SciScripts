##################################################################
# fit a potential energy curve to a Morse potential
##################################################################

using DelimitedFiles: writedlm, readdlm
using LsqFit

global PROGRAM = basename(@__FILE__)

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

function die(message::String)
    usage()
    error(message)
end # function die

function Morse(r :: Union{Number, AbstractArray{<:Number}}, c :: AbstractArray{<:Number})
    a  = c[1]
    r0 = c[2] # equilibrium distance
    D  = c[3] # depth
    diss = c[4] # dissociation limit
    0 in r ? die("At least one of the supplied internuclear distances is 0") :
        @. return D * (exp(-2 * a * (r - r0)) - 2 * exp(-a * (r - r0))) + diss
end # function Morse

function main(
        infile      :: String,
        outfile     :: String,
        headOnly    :: Bool,
        tailOnly    :: Bool,
        rCutoff     :: Float64,
        rMin        :: Float64,
        rMax        :: Float64,
        rStep       :: Float64,
        a           :: Number,
        r0          :: Number,
        D           :: Number,
        diss        :: Number,
        printParams :: Bool
    )

    if(headOnly && tailOnly)
        error("tailOnly and headOnly cannot both be true")
    end

    # -- read data into x and y arrays
    data = readdlm(infile)
    xdata = data[:,1]
    ydata = data[:,2]
    data = 0 # <-- clear data (it shouldn't be that big anyway)

    if tailOnly
        y = ydata[xdata .> rCutoff]
        x = xdata[xdata .> rCutoff]
    elseif headOnly
        y = ydata[xdata .< rCutoff]
        x = xdata[xdata .< rCutoff]
    else
        x = xdata
        y = ydata
    end

    c = [a, r0, D, diss]

    fit = curve_fit(Morse, x, y, c)

    if tailOnly
        rMin = last(x)
    elseif headOnly
        rMax = x[1]
    end

    xfit = [r for r in rMin:rStep:rMax]
    yfit = Morse(xfit, fit.param)

    if tailOnly
        xfit = [xdata; xfit[2:end]]
        yfit = [ydata; yfit[2:end]]
    elseif headOnly
        xfit = [xfit[1:end - 1]; xdata]
        yfit = [yfit[1:end - 1]; ydata]
    end

    writedlm(outfile, [xfit yfit])

    if printParams
        println("a, r0, D, dissociation limit  = ",  fit.param)
    end

end # fuction main

nargs = 13

if(size(ARGS,1) != nargs)
    println("STOP: ", PROGRAM, " requires ", nargs, " arguments, but ", string(size(ARGS, 1)), " were supplied.")
    usage(1)
    exit()
end

infile      = ARGS[1]                  # -- name of the input file containing the potential data
outfile     = ARGS[2]                  # -- name of the output file
headOnly    = parse(Bool, ARGS[3])     # -- (T) only fit the tail of the potential (F) fit the whole thing
tailOnly    = parse(Bool, ARGS[4])     # -- (T) only fit the tail of the potential (F) fit the whole thing
rCutoff     = parse(Float64, ARGS[5])  # -- cutoff distance for fitting the tail. Values smaller than this are ignored
rMin        = parse(Float64, ARGS[6])  # -- smallest distance for theARGS[5] fit. Ignored if tailOnly is true
rMax        = parse(Float64, ARGS[7])  # -- largest  distance for theARGS[6] fit
rStep       = parse(Float64, ARGS[8])  # -- linear step size for the ARGS[7]fit
a           = parse(Float64, ARGS[9])  # -- potential width
r0          = parse(Float64, ARGS[10])  # -- well depth
D           = parse(Float64, ARGS[11]) # -- well minimum
diss        = parse(Float64, ARGS[12]) # -- dissociation limit
printParams = parse(Bool, ARGS[13])    # -- (T) print the optimized potential parameters (F) don't

main(infile, outfile, headOnly, tailOnly, rCutoff, rMin, rMax, rStep, a, r0, D, diss, printParams)
