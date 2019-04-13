
"""
`module  JAC.PairAnnihilation1Photon`  ... a submodel of JAC that contains all methods for computing positron-bound-electron pair annihilation 
                                           (PEPA) with two-photon emission cross sections and rates; e^+ + |i(N)> --> |f(N-1)> + photon + photon;
                                           it is using JAC, JAC.ManyElectron. 
"""
module PairAnnihilation2Photon 

    using JAC, JAC.ManyElectron


    """
    `struct  PairAnnihilation2Photon.Settings`  ... defines a type for the details and parameters in computing positron-bound-electron pair 
             annihilation (PEPA) with single-photon emission cross sections and rates; e^+ + |i(N)> --> |f(N-1)> + photon1 + photon2.

        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + positronEnergies        ::Array{Float64,1}             ... List of positron energies.
        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).
    """
    struct Settings 
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        positronEnergies          ::Array{Float64,1} 
        printBeforeComputation    ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1} 
    end 


    """
    `JAC.PairAnnihilation2Photon.Settings()`  ... constructor for the default values of pair-annihilation photon line computations
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], false, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::PairAnnihilation2Photon.Settings)`  
    #		... prepares a proper printout of the variable settings::PairAnnihilation2Photon.Settings.
    function Base.show(io::IO, settings::PairAnnihilation2Photon.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "positronEnergies:         $(settings.positronEnergies)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
    end



    """
    `struct  PairAnnihilation2Photon.Channel`  ... defines a type for a positron-bound-electron pair annihilation (PEPA) with two-photon 
                                                   emission channel that specifies all quantum numbers, phases and amplitudes.

        + multipole1     ::EmMultipole          ... Multipole of the emitted photon1.
        + multipole2     ::EmMultipole          ... Multipole of the emitted photon2.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + kappa          ::Int64                ... partial-wave of the incoming free positron
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + phase          ::Float64              ... phase of the partial wave
        + amplitude      ::Complex{Float64}     ... PairAnnihilation2PhotonChannel amplitude associated with the given channel.
   """
    struct  Channel
        multipole1       ::EmMultipole
        multipole2       ::EmMultipole
        gauge            ::EmGauge
        kappa            ::Int64
        symmetry         ::LevelSymmetry
        phase            ::Float64
        amplitude        ::Complex{Float64}
    end


    """
    `struct  PairAnnihilation2Photon.Line`  ... defines a type for a positron-bound-electron pair-annihilation (two-photon) line that may 
                                                include the definition of channels.

        + initialLevel   ::Level                 ... initial-(state) level
        + finalLevel     ::Level                 ... final-(state) level
        + positronEnergy ::Float64                ... Energy of the (incoming free) positron.
        + photon1Energy  ::Float64                ... Energy of the emitted photon1.
        + photon2Energy  ::Float64                ... Energy of the emitted photon2.
        + crossSection   ::EmProperty             ... Cross section for this pair-annihilation (photon) line.
        + hasChannels    ::Bool                   ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                      free-positron energy, kappa, multipole, etc., or not.
        + channels       ::Array{PairAnnihilation2Photon.Channel,1}  ... List of PairAnnihilation2Photon.Channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        positronEnergy   ::Float64
        photon1Energy    ::Float64
        photon2Energy    ::Float64
        crossSection     ::EmProperty
        hasChannels      ::Bool
        channels         ::Array{PairAnnihilation2Photon.Channel,1}
    end



    """
    `JAC.PairAnnihilation2Photon.Line()`  ... 'empty' constructor for a pair-annihilation (photon) line between a specified initial and 
                                              final level.
    """
    function Line()
        Line(Level(), Level(), 0., 0., 0., EmProperty(0., 0.), false, PairAnnihilation2Photon[] )
    end


    # `Base.show(io::IO, line::PairAnnihilation2Photon.Line)`  ... prepares a proper printout of the variable line::PairAnnihilation2Photon.Line.
    function Base.show(io::IO, line::PairAnnihilation1Photon.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "positronEnergy:    $(line.positronEnergy)  ")
        println(io, "photon1Energy:     $(line.photon1Energy)  ")
        println(io, "photon2Energy:     $(line.photon2Energy)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
        println(io, "hasChannels:       $(line.hasChannels)  ")
        println(io, "channels:          $(line.channels)  ")
    end

end # module
