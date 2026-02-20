# using NumericalIntegration: integrate
using DelimitedFiles: writedlm, readdlm
using FStrings

###############################################################################
#                                                                             #
#    scattering cross sections (area) to a rate coefficient (area/time)       #
#                                                                             #
#    reads a data file formatted as :                                         #
#                                                                             #
#      x  y                                                                   #
#      .  .                                                                   #
#      .  .                                                                   #
#      .  .                                                                   #
#                                                                             #
#    where x is an energy and y is a cross section.                           #
#    Should be used with the wrapper script by the same name.                 #
#                                                                             #
#    Can also be used to convert rate coefficients to cross sections          #
#                                                                             #
###############################################################################

PROGRAM = basename(@__FILE__)

au2cm = 0.5291772083e-8 # <---- atomic units -> centimeters (length)
au2s  = 2.4188843e-17 # <------ atomic units -> seconds (time)
au2rate = au2cm^3 / au2s # <-- atomic units -> centimeters^3 per second (rate coefficient)
au2K = 3.1577465e5 # <--------- atomic units -> Kelvin (temperature)

function usage(::Int)
    println()
    println(f"{PROGRAM}: Convolves data from a datafile with a function ")
    println()
    println(f"    usage: {PROGRAM} [ARGSLIST]")
    println()
    return true
end

function main()

    narg = 9

    if(size(ARGS,1) != narg)
        println(f"STOP: {PROGRAM} requires {narg} arguments, but {string(size(ARGS,1))} were supplied.")
        usage(1)
        exit()
    end

    datafile             = ARGS[1]
    outfile              = ARGS[2]
    xsConversion         = parse(Float64, ARGS[3]) # <-- : multiplicative conversion factor supplied by wrapper to get cross section length in bohr
    rateConversionLength = parse(Float64, ARGS[4]) # <-- : multiplicative conversion factor supplied by wrapper to get the rate coefficient length in bohr
    rateConversionTime   = parse(Float64, ARGS[5]) # <-- : multiplicative conversion factor supplied by wrapper to get the rate coefficient time in seconds
    EConversion          = parse(Float64, ARGS[6]) # <-- : multiplicative conversion factor supplied by wrapper to get energy in hartree atomic units
    comment_char         = only(ARGS[7]) # <-- : convert string to character
    ignore_comments_ARG  = ARGS[8]
    input_type           = ARGS[9]

    ignore_comments = ignore_comments_ARG == "false" ? false : true

    # -- check if input datafile exists
    if( ! isfile(datafile) )
        usage(1)
        println("STOP: Need a data file to convole")
        exit()
    end

    # -- assumes x y format, one y column (for now)
    data  = readdlm(datafile, comments = ignore_comments, comment_char = comment_char)
    xdata = data[:, 1] * EConversion

    # -- sort electron energies by size if not done so already
    indices = sortperm(xdata)
    EelGrid = xdata[indices]

    if input_type == "XS"

        # -- cross sections -> rate coefficient
        output = data[:, 2] * xsConversion ^ 2 .* sqrt.(2 .* EelGrid) .* au2rate

    elseif input_type == "RATE"

        # -- rate coefficient -> cross sections
        output = data[:, 2] * rateConversionLength ^ 3 / rateConversionTime  * au2s ./ sqrt.(2 .* EelGrid) ./ xsConversion^2

    else

        usage(1)
        println("STOP: input_type does not match either XS or RATE")
        exit()

    end

    writedlm(outfile, [data[:, 1] output])

end

main()
