
"""
`module  JAC.RayleighCompton`  ... a submodel of JAC that contains all methods for computing elastic Rayleigh and inelastic Compton 
                                    photon scattering cross sections; it is using JAC, JAC.ManyElectron.
"""
module RayleighCompton

    using JAC.Basics, JAC.ManyElectron, JAC.Radial


    """
    `struct  RayleighCompton.Settings`  ... defines a type for the settings in estimating lastic Rayleigh and inelastic Compton 
                                            photon scattering cross sections

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
    `JAC.RayleighCompton.Settings()`  ... constructor for the default values of Rayleigh-Compton photon-scattering estimates.
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], false, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::RayleighCompton.Settings)`  ... prepares a proper printout of the variable settings::RayleighCompton.Settings.
    function Base.show(io::IO, settings::RayleighCompton.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "photonEnergies:           $(settings.photonEnergies)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
    end


    """
    `struct  RayleighCompton.Channel`  ... defines a type for a Rayleigh-Compton channel to help characterize the scattering of a single photon
                                           with well-defined multipolarities.

        + inMultipole    ::EmMultipole          ... Multipole of the incoming photon.
        + outMultipole   ::EmMultipole          ... Multipole of the outgoing, scattered photon.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + amplitude      ::Complex{Float64}     ... Photoionization amplitude associated with the given channel.
   """
    struct  Channel
        inMultipole      ::EmMultipole
        outMultipole     ::EmMultipole
        gauge            ::EmGauge
        amplitude        ::Complex{Float64}
    end


    """
    `struct  RayleighCompton.Line`  ... defines a type for a Rayleigh-Compton scattering line that may include the definition of channels.

        + initialLevel     ::Level           ... initial-(state) level
        + finalLevel       ::Level           ... final-(state) level
        + inPhotonEnergy   ::Float64         ... Energy of the incoming photon.
        + outPhotonEnergy  ::Float64         ... Energy of the outgoing, scattered photon.
        + crossSection     ::EmProperty      ... Cross section for this photoionization.
        + hasChannels      ::Bool            ... Determines whether the individual (sub-) channels are defined in terms of their multipolarities, 
                                                 etc., or not.
        + channels         ::Array{RayleighCompton.Channel,1}  ... List of RayleighCompton.Channels of this line.
    """
    struct  Line
        initialLevel       ::Level
        finalLevel         ::Level
        inPhotonEnergy     ::Float64
        outPhotonEnergy    ::Float64
        crossSection       ::EmProperty
        hasChannels        ::Bool
        channels           ::Array{RayleighCompton.Channel,1}
    end


    # `Base.show(io::IO, line::RayleighCompton.Line)`  ... prepares a proper printout of the variable line::RayleighCompton.Line.
    function Base.show(io::IO, line::RayleighCompton.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "inPhotonEnergy:    $(line.inPhotonEnergy)  ")
        println(io, "outPhotonEnergy:   $(line.outPhotonEnergy)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
        println(io, "hasChannels:       $(line.hasChannels)  ")
        println(io, "channels:          $(line.channels)  ")
    end



    """
    `JAC.RayleighCompton.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::RayleighCompton.Settings; output=true)`  
               ... to compute the Rayleigh-Compton transition amplitudes and all properties as requested by the given settings. A list of 
                   lines::Array{RayleighCompton.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::RayleighCompton.Settings; output=true)
        println("")
        printstyled("JAC.RayleighCompton.computeLines(): The computation of Rayleigh-ComptonAuger cross sections starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------------- \n", color=:light_green)
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
