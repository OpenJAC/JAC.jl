
"""
`module  JAC.ImpactIonization`  ... a submodel of JAC that contains all methods for computing electron impact ionization cross sections; 
                                    it is using JAC, JAC.ManyElectron, JAC.ImpactExcitation, JAC.AutoIonization.
"""
module ImpactIonization

    using JAC, JAC.ManyElectron, JAC.ImpactExcitation, JAC.AutoIonization


    """
    `struct  ImpactIonization.Settings`  ... defines a type for the settings in estimating electron-impact ionization cross sections.

        + electronEnergies   ::Array{Float64,1}             ... List of impact-energies of the incoming elecgtrons.
        + calcShellDependent ::Bool                         ... True of shell-dependent cross sections are to be calculated, and false otherwise.
        + shells             ::Array{Shell,1}               ... List of shells for which shell-dependent cross sections are calculated.
    """
    struct Settings
        electronEnergies     ::Array{Float64,1}
        calcShellDependent   ::Bool 
        shells               ::Array{Shell,1} 
    end 


    """
    `JAC.ImpactIonization.Settings()`  ... constructor for the default values of electron-impact ionization cross section estimates.
    """
    function Settings()
       Settings( Float64[], false, Shell[] )
    end


    ## `Base.show(io::IO, settings::ImpactIonization.Settings)`  ... prepares a proper printout of the variable settings::ImpactIonization.Settings.
    function Base.show(io::IO, settings::ImpactIonization.Settings) 
        println(io, "electronEnergies:       $(settings.electronEnergies)  ")
        println(io, "calcShellDependent:     $(settings.calcShellDependent)  ")
        println(io, "shells:                 $(settings.shells)  ")
    end


    """
    `struct  ImpactIonization.Channel`  ... defines a type for a electron-impact ionization channel to help characterize the incoming and 
                                            two outgoing (continuum) states of many electron-states with (two) free electron

        + initialEpsilon   ::Float64            ... energy of the incoming free-electron
        + finalEpsilon     ::Float64            ... energy of the outgoing free-electron
        + releasedEpsilon  ::Float64            ... energy of the released free-electron
        + initialKappa     ::Int64              ... partial-wave of the incoming free electron
        + finalKappa       ::Int64              ... partial-wave of the outgoing free electron
        + releasedKappa    ::Int64              ... partial-wave of the released free electron
        + symmetry         ::LevelSymmetry      ... total angular momentum and parity of the scattering state
        + initialPhase     ::Float64            ... phase of the incoming partial wave
        + finalPhase       ::Float64            ... phase of the outgoing partial wave
        + releasedPhase    ::Float64            ... phase of the released partial wave
        + amplitude        ::Complex{Float64}   ... Collision amplitude associated with the given channel.
    """
    struct  Channel
        initialEpsilon     ::Float64
        finalEpsilon       ::Float64
        releasedEpsilon    ::Float64
        initialKappa       ::Int64 
        finalKappa         ::Int64 
        releasedKappa      ::Int64 
        symmetry           ::LevelSymmetry
        initialPhase       ::Float64
        finalPhase         ::Float64
        releasedPhase      ::Float64
        amplitude          ::Complex{Float64}
    end


    # `Base.show(io::IO, channel::ImpactIonization.Channel)`  ... prepares a proper printout of the variable channel::ImpactIonization.Channel.
    function Base.show(io::IO, channel::ImpactIonization.Channel) 
        println(io, "initialEpsilon:     $(channel.initialEpsilon)  ")
        println(io, "finalEpsilon:       $(channel.finalEpsilon)  ")
        println(io, "releasedEpsilon:    $(channel.releasedEpsilon)  ")
        println(io, "initialKappa:       $(channel.initialKappa)  ")
        println(io, "finalKappa:         $(channel.finalKappa)  ")
        println(io, "releasedKappa:      $(channel.releasedKappa)  ")
        println(io, "symmetry:           $(channel.symmetry)  ")
        println(io, "initialPhase:       $(channel.initialPhase)  ")
        println(io, "finalPhase:         $(channel.finalPhase)  ")
        println(io, "releasedPhase:      $(channel.releasedPhase)  ")
        println(io, "amplitude:          $(channel.amplitude)  ")
    end


    """
    `struct  ImpactIonization.Line`  ... defines a type for a electron-impact ionization line that may include the definition of channels and 
                                         their corresponding amplitudes.

        + initialLevel   ::Level          ... initial-(state) level
        + finalLevel     ::Level          ... final-(state) level
        + crossSection   ::Float64        ... total cross section of this line
        + hasChannels    ::Bool           ... Determines whether the individual scattering (sub-) channels are defined in terms of their free-
                                              electron energies, kappa's, phases and the total angular momentum/parity as well as the amplitude, 
                                              or not.
        + channels       ::Array{ImpactIonization.Line,1}  ... List of Eimex channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        crossSection     ::Float64 
        hasChannels      ::Bool
        channels         ::Array{ImpactIonization.Channel,1}
    end 


    # `Base.show(io::IO, line::ImpactIonization.Line)`  ... prepares a proper printout of the variable line::ImpactIonization.Line.
    function Base.show(io::IO, line::ImpactIonization.Line) 
        println(io, "initialLevel:     $(line.initialLevel)  ")
        println(io, "finalLevel:       $(line.finalLevel)  ")
        println(io, "crossSection:     $(line.crossSection)  ")
        println(io, "hasChannels:      $(line.hasChannels)  ")
        println(io, "channels:         $(line.channels)  ")
    end

end # module
