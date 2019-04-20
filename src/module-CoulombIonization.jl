
"""
`module  JAC.CoulombIonization`  ... a submodel of JAC that contains all methods for computing Coulomb ionization cross sections and alignment
                                     parameters for the excitation of target or projectile electrons by fast ion impact; it is using JAC, 
                                     JAC.ManyElectron, JAC.Radial.
"""
module CoulombIonization

    using JAC.BasicTypes, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0


    """
    `struct  CoulombIonization.Settings`  ... defines a type for the details and parameters of computing Coulomb-ionization lines.

        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + energies                ::Array{Float64,1}             ... List of ... energies.
        + calcAlignment           ::Bool                         ... True, if alignment parameters to be calculated and false otherwise.
        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).
    """
    struct Settings 
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        photonEnergies            ::Array{Float64,1} 
        calcAlignment             ::Bool 
        printBeforeComputation    ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1} 
    end 


    """
    `JAC.CoulombIonization.Settings()`  ... constructor for the default values of Coulomb-excitation line computations.
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], false, false, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::CoulombIonization.Settings)`  ... prepares a proper printout of the variable settings::CoulombIonization.Settings.
    function Base.show(io::IO, settings::CoulombIonization.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "photonEnergies:           $(settings.photonEnergies)  ")
        println(io, "calcAlignment:            $(settings.calcAlignment)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
    end


    """
    `struct  CoulombIonization.Channel`  ... defines a type for a Coulomb-excitation channel to help characterize a single magnetic substate.
                                             !!! This need to adapted !!!

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
    `struct  Line`  ... defines a type for a Coulomb-excitation line that may include the definition of channels.

        + initialLevel   ::Level                  ... initial-(state) level
        + finalLevel     ::Level                  ... final-(state) level
        + crossSection   ::EmProperty             ... Cross section for this photoionization.
        + alignmentA2    ::EmProperty             ... Alignment A_2 parameter.
        + hasChannels    ::Bool                   ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                      magnetic substates, etc., or not.
        + channels       ::Array{CoulombIonization.Channel,1}  ... List of CoulombIonization.Channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        crossSection     ::EmProperty
        alignmentA2      ::EmProperty
        hasChannels      ::Bool
        channels         ::Array{CoulombIonization.Channel,1}
    end


    # `Base.show(io::IO, line::CoulombIonization.Line)`  ... prepares a proper printout of the variable line::CoulombIonization.Line.
    function Base.show(io::IO, line::CoulombIonization.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
    end

end # module
