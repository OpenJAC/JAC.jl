
"""
`module  JAC.ResonantTwoColorPhotonIon`  
    ... a submodel of JAC that contains all methods for computing resonant two-color (two-photon & single-electron) ionization cross 
        sections.
"""
module ResonantTwoColorPhotonIon

    using ..ManyElectron


    """
    `struct  ResonantTwoColorPhotonIon.Settings`  
        ... defines a type for the settings in estimating resonant two-color (two-photon & single-electron) ionization cross sections.

        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + photonEnergies          ::Array{Float64,1}             ... List of photon energies.
        + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.
    """
    struct Settings 
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        photonEnergies            ::Array{Float64,1} 
        printBefore               ::Bool 
        lineSelection             ::LineSelection 
   end 


    """
    `ResonantTwoColorPhotonIon.Settings()`  ... constructor for the default values of resonant two-color (two-photon & single-electron) 
         ionization estimates.
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], false, LineSelection() )
    end


    # `Base.show(io::IO, settings::ResonantTwoColorPhotonIon.Settings)`  
    #		... prepares a proper printout of the variable settings::ResonantTwoColorPhotonIon.Settings.
    function Base.show(io::IO, settings::ResonantTwoColorPhotonIon.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "photonEnergies:           $(settings.photonEnergies)  ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "lineSelection:            $(settings.lineSelection)  ")
    end


    """
    `struct  ResonantTwoColorPhotonIon.Channel`  ... defines a type for a two-color channel to help characterize resonant two-color (two-photon & 
             single-electron) ionization with well-defined multipolarities.

        + multipoles     ::Array{EmMultipole,1}   ... Multipoles of all N incoming photons.
        + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
        + kappa          ::Int64                  ... partial-wave of the free electron
        + symmetry       ::LevelSymmetry          ... total angular momentum and parity of the scattering state
        + phase          ::Float64                ... phase of the partial wave
        + amplitude      ::Complex{Float64}       ... multi-photon ionization amplitude associated with the given channel.
   """
    struct  Channel
        multipoles       ::Array{EmMultipole,1}
        gauge            ::EmGauge 
        kappa            ::Int64
        symmetry         ::LevelSymmetry 
        phase            ::Float64
        amplitude        ::Complex{Float64}
    end



    """
    `struct  ResonantTwoColorPhotonIon.Line`  ... defines a type for a multi-photon ionization line that may include the definition of channels.

        + initialLevel     ::Level           ... initial-(state) level
        + finalLevel       ::Level           ... final-(state) level
        + photonEnergy     ::Float64         ... Energy of the incoming photons; all photons are assumed to have equal energy.
        + crossSection     ::EmProperty      ... Cross section for this multi-photon ionization.
        + channels         ::Array{ResonantTwoColorPhotonIon.Channel,1}  ... List of ResonantTwoColorPhotonIon.Channels of this line.
    """
    struct  Line
        initialLevel       ::Level
        finalLevel         ::Level
        photonEnergy       ::Float64
        crossSection       ::EmProperty
        channels           ::Array{ResonantTwoColorPhotonIon.Channel,1}
    end


    # `Base.show(io::IO, line::ResonantTwoColorPhotonIon.Line)`  ... prepares a proper printout of the variable line::ResonantTwoColorPhotonIon.Line.
    function Base.show(io::IO, line::ResonantTwoColorPhotonIon.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "photonEnergy:      $(line.photonEnergy)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
       println(io, "channels:          $(line.channels)  ")
    end

end # module
