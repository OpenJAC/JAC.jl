
"""
`module  JAC.CoulombExcitation`  ... a submodel of JAC that contains all methods for computing Coulomb excitation cross sections and alignment
                                     parameters for the excitation of target or projectile electrons by fast ion impact; it is using JAC.Basics. 
                                     JAC.ManyElectron, JAC.Radial.
"""
module CoulombExcitation

    using JAC.Basics, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0


    """
    `struct  CoulombExcitation.Settings`  ... defines a type for the details and parameters of computing Coulomb-excitation lines.

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
    `JAC.CoulombExcitation.Settings()`  ... constructor for the default values of Coulomb-excitation line computations.
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], false, false, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::CoulombExcitation.Settings)`  ... prepares a proper printout of the variable settings::CoulombExcitation.Settings.
    function Base.show(io::IO, settings::CoulombExcitation.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "photonEnergies:           $(settings.photonEnergies)  ")
        println(io, "calcAlignment:            $(settings.calcAlignment)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
    end


    """
    `struct  CoulombExcitation.Channel`  ... defines a type for a Coulomb-excitation channel to help characterize a single magnetic substate.
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
        + channels       ::Array{CoulombExcitation.Channel,1}  ... List of CoulombExcitation.Channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        crossSection     ::EmProperty
        alignmentA2      ::EmProperty
        hasChannels      ::Bool
        channels         ::Array{CoulombExcitation.Channel,1}
    end


    # `Base.show(io::IO, line::CoulombExcitation.Line)`  ... prepares a proper printout of the variable line::CoulombExcitation.Line.
    function Base.show(io::IO, line::CoulombExcitation.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
    end



    """
    `JAC.CoulombExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::CoulombExcitation.Settings; output=true)`  
               ... to compute the Coulomb excitation amplitudes and all properties as requested by the given settings. A list of 
                   lines::Array{CoulombExcitation.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::CoulombExcitation.Settings; output=true)
        println("")
        printstyled("JAC.CoulombExcitation.computeLines(): The computation of Coulomb excitation cross sections starts now ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        println("Not yet implemented: Data structures and properties still need to be worked out.")
        #
        lines = "Not yet implemented !"
        if    output    return( lines )
        else            return( nothing )
        end
    end

end # module
