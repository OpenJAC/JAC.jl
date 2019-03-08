
"""
`module  JAC.DoubleAuger`  ... a submodel of JAC that contains all methods for computing double Auger and autoionization amplitudes 
                               and rates; it is using JAC, JAC.ManyElectron, JAC.ImpactExcitation, JAC.Auger.
"""
module DoubleAuger

    using JAC, JAC.ManyElectron, JAC.ImpactExcitation, JAC.Auger


    """
    `struct  DoubleAuger.Settings`  ... defines a type for the settings in estimating double-Auger and autoionization rates.

        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).
        + minAugerEnergy          ::Float64                      ... Minimum energy of free (Auger) electrons to be included.
        + maxAugerEnergy          ::Float64                      ... Maximum energy of free (Auger) electrons to be included.
        + maxKappa                ::Int6464                      ... Maximum kappa value of partial waves to be included.
    """
    struct Settings
        printBeforeComputation    ::Bool
        selectLines               ::Bool  
        selectedLines             ::Array{Tuple{Int64,Int64},1}
        minAugerEnergy            ::Float64
        maxAugerEnergy            ::Float64
        maxKappa                  ::Int64
    end 


    """
    `JAC.DoubleAuger.Settings()`  ... constructor for the default values of DoubleAuger line computations
    """
    function Settings()
        Settings(false, false, Array{Tuple{Int64,Int64},1}[], 0., 10e5, 100)
    end


    """
    `Base.show(io::IO, settings::DoubleAuger.Settings)`  ... prepares a proper printout of the variable settings::DoubleAuger.Settings.
    """
    function Base.show(io::IO, settings::DoubleAuger.Settings) 
        println(io, "printBeforeComputation:       $(settings.printBeforeComputation)  ")
        println(io, "selectLines:                  $(settings.selectLines)  ")
        println(io, "selectedLines:                $(settings.selectedLines)  ")
        println(io, "minAugerEnergy:               $(settings.minAugerEnergy)  ")
        println(io, "maxAugerEnergy:               $(settings.maxAugerEnergy)  ")
        println(io, "maxKappa:                     $(settings.maxKappa)  ")
    end


    """
    `struct  Channel`  ... defines a type for a DoubleAuger channel to help characterize a scattering (continuum) state of many 
                           electron-states with two free electrons.

        + kappa1         ::Int64                ... partial-wave of the free electron_1
        + kappa2         ::Int64                ... partial-wave of the free electron_2
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + phase1         ::Float64              ... phase of the partial wave_1
        + phase2         ::Float64              ... phase of the partial wave_2
        + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
    """
    struct  Channel
        kappa1           ::Int64
        kappa2           ::Int64
        symmetry         ::LevelSymmetry
        phase1           ::Float64
        phase2           ::Float64
        amplitude        ::Complex{Float64}
    end


    """
    `struct  Line`  ... defines a type for a Auger line that may include the definition of sublines and their corresponding amplitudes.

        + initialLevel   ::Level         ... initial-(state) level
        + finalLevel     ::Level         ... final-(state) level
        + electronEnergy ::Float64       ... Energy of the (incoming free) electron.
        + totalRate      ::Float64       ... Total rate of this line.
        + hasChannels    ::Bool          ... Determines whether the individual scattering (sub-) channels are defined in terms of their free-
                                             electron energies, kappas and the total angular momentum/parity as well as the amplitude, or not.
        + channels       ::Array{DoubleAuger.Channel,1}  ... List of DoubleAuger channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        electronEnergy   ::Float64
        totalRate        ::Float64
        hasChannels      ::Bool
        channels         ::Array{DoubleAuger.Channel,1}
    end 


    """
    `Base.show(io::IO, line::DoubleAuger.Line)`  ... prepares a proper printout of the variable line::DoubleAuger.Line.
    """
    function Base.show(io::IO, line::DoubleAuger.Line) 
        println(io, "initialLevel:           $(line.initialLevel)  ")
        println(io, "finalLevel:             $(line.finalLevel)  ")
        println(io, "electronEnergy:         $(line.electronEnergy)  ")
        println(io, "totalRate:              $(line.totalRate)  ")
        println(io, "hasChannels:            $(line.hasChannels)  ")
        println(io, "channels:               $(line.channels)  ")
    end

end # module
