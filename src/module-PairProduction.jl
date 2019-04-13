
"""
`module  JAC.PairProduction`  ... a submodel of JAC that contains all methods for computing positron-bound-electron pair production (PEPP) 
                                  by single-photon absorption cross sections and rates; photon + |i(N)> --> |f(N+1)> + e^+;
                                  it is using JAC, JAC.ManyElectron. 
"""
module PairProduction 

    using JAC, JAC.ManyElectron


    """
    `struct  PairProduction.Settings`  ... defines a type for the details and parameters in computing positron-bound-electron pair production 
                                           (PEPP) with single-photon absorption (positron) lines; photon + |i(N)> --> |f(N+1)> + e^+.

        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + photonEnergies          ::Array{Float64,1}             ... List of photon energies.
        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).
    """
    struct Settings 
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        photonEnergies            ::Array{Float64,1} 
        printBeforeComputation    ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1} 
    end 


    """
    `JAC.PairProduction.Settings()`  ... constructor for the default values of pair-production positron line computations
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], false, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::PairProduction.Settings)`  ... prepares a proper printout of the variable settings::PairProduction.Settings.
    function Base.show(io::IO, settings::PairProduction.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "photonEnergies:           $(settings.photonEnergies)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
    end


    """
    `struct  PairProduction.Channel`  ... defines a type for a positron-bound-electron pair production (PEPP) with single-photon absorption
                                          (positron) channel that specifies all quantum numbers, phases and amplitudes.

        + multipole      ::EmMultipole          ... Multipole of the photon absorption.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + kappa          ::Int64                ... partial-wave of the free electron
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + phase          ::Float64              ... phase of the partial wave
        + amplitude      ::Complex{Float64}     ... Photoionization amplitude associated with the given channel.
   """
    struct  Channel
        multipole        ::EmMultipole
        gauge            ::EmGauge
        kappa            ::Int64
        symmetry         ::LevelSymmetry
        phase            ::Float64
        amplitude        ::Complex{Float64}
    end


    """
    `struct  PairProduction.Line`  ... defines a type for a pair-production (positron) line that may include the definition of channels.

        + initialLevel   ::Level                  ... initial-(state) level
        + finalLevel     ::Level                  ... final-(state) level
        + positronEnergy ::Float64                ... Energy of the (outgoing free) positron.
        + photonEnergy   ::Float64                ... Energy of the absorbed photon.
        + crossSection   ::EmProperty             ... Cross section for this pair-production (positron) line.
        + hasChannels    ::Bool                   ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                      free-electron energy, kappa, multipole, etc., or not.
        + channels       ::Array{PairProduction.Channel,1}  ... List of PairProduction.Channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        positronEnergy   ::Float64
        photonEnergy     ::Float64
        crossSection     ::EmProperty
        hasChannels      ::Bool
        channels         ::Array{PairProduction.Channel,1}
    end


    """
    `JAC.PairProduction.Line()`  ... 'empty' constructor for a pair-production (positron) line between a specified initial and final level.
    """
    function Line()
        Line(Level(), Level(), 0., 0., EmProperty(0., 0.), false, PairProduction[] )
    end


    # `Base.show(io::IO, line::PairProduction.Line)`  ... prepares a proper printout of the variable line::PairProduction.Line.
    function Base.show(io::IO, line::PairProduction.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "positronEnergy:    $(line.positronEnergy)  ")
        println(io, "photonEnergy:      $(line.photonEnergy)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
        println(io, "hasChannels:       $(line.hasChannels)  ")
        println(io, "channels:          $(line.channels)  ")
    end

end # module
