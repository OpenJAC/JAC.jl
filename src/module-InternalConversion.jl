
"""
`module  JAC.InternalConversion`  ... a submodel of JAC that contains all methods for computing photo-excitation-autoionization cross 
                                          sections and rates; it is using JAC, JAC.ManyElectron, JAC.Radial, JAC.PhotoEmission, JAC.AutoIonization.
"""
module InternalConversion 

    using JAC.Basics, JAC.ManyElectron, JAC.Radial, JAC.PhotoEmission, JAC.AutoIonization
    global JAC_counter = 0


    """
    `struct  InternalConversion.Settings`  ... defines a type for the details and parameters of computing photon-impact 
                                                   excitation-autoionization pathways |i(N)>  --> |m(N)>  --> |f(N-1)>.

        + multipoles              ::Array{EmMultipole,1}               ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge,1}                  ... Specifies the gauges to be included into the computations.
        + printBeforeComputation  ::Bool                               ... True, if all energies and lines are printed before their evaluation.
        + selectPathways          ::Bool                               ... True if particular pathways are selected for the computations.
        + selectedPathways        ::Array{Tuple{Int64,Int64,Int64},1}  ... List of list of pathways, given by tupels (inital, inmediate, final).
        + maxKappa                ::Int64                              ... Maximum kappa value of partial waves to be included.
    """
    struct Settings
        multipoles                ::Array{EmMultipole,1}
        gauges                    ::Array{UseGauge,1} 
        printBeforeComputation    ::Bool
        selectPathways            ::Bool
        selectedPathways          ::Array{Tuple{Int64,Int64,Int64},1}
        maxKappa                  ::Int64 
    end 


    """
    `JAC.InternalConversion.Settings()`  ... constructor for the default values of photon-impact excitation-autoionizaton settings.
    """
    function Settings()
        Settings( JAC.EmMultipole[], UseGauge[], false,  false, Tuple{Int64,Int64,Int64}[], 0)
    end


    # `Base.show(io::IO, settings::InternalConversion.Settings)`  
    #   ... prepares a proper printout of the variable settings::InternalConversion.Settings.  
    function Base.show(io::IO, settings::InternalConversion.Settings) 
        println(io, "multipoles:              $(settings.multipoles)  ")
        println(io, "gauges:                  $(settings.gauges)  ")
        println(io, "printBeforeComputation:  $(settings.printBeforeComputation)  ")
        println(io, "selectPathways:          $(settings.selectPathways)  ")
        println(io, "selectedPathways:        $(settings.selectedPathways)  ")
        println(io, "maxKappa:                $(settings.maxKappa)  ")
    end



    """
    `struct  JAC.InternalConversion.Channel`  ... defines a type for a photon-impact excitaton & autoionization channel that specifies 
                                                      all quantum numbers, phases and amplitudes.

        + excitationChannel  ::PhotoEmission.Channel       ... Channel that describes the photon-impact excitation process.
        + augerChannel       ::AutoIonization.Channel  ... Channel that describes the subsequent Auger/autoionization process.
    """
    struct  Channel
        excitationChannel    ::PhotoEmission.Channel
        augerChannel         ::AutoIonization.Channel
    end 


    """
    `struct  JAC.InternalConversion.Pathway`  ... defines a type for a photon-impact excitation pathway that may include the definition 
                                                      of different excitation and autoionization channels and their corresponding amplitudes.

        + initialLevel        ::Level                  ... initial-(state) level
        + intermediateLevel   ::Level                  ... intermediate-(state) level
        + finalLevel          ::Level                  ... final-(state) level
        + photonEnergy        ::Float64                 ... energy of the (incoming) electron
        + electronEnergy      ::Float64                 ... energy of the (finally outgoing, scattered) electron
        + crossSection        ::EmProperty              ... total cross section of this pathway
        + hasChannels         ::Bool                    ... Determines whether the individual excitation and autoionization channels are defined 
                                                            in terms of their multipole, gauge, free-electron kappa, phases and the total 
                                                            angular momentum/parity as well as the amplitude, or not.
        + channels            ::Array{InternalConversion.Channel,1}  ... List of channels of this pathway.
    """
    struct  Pathway
        initialLevel          ::Level
        intermediateLevel     ::Level
        finalLevel            ::Level
        photonEnergy          ::Float64
        electronEnergy        ::Float64
        crossSection          ::EmProperty
        hasChannels           ::Bool
        channels              ::Array{InternalConversion.Channel,1}
    end 


    """
    `JAC.InternalConversion.Pathway()`  ... 'empty' constructor for an photon-impact excitation-autoionization pathway between a specified 
                                                initial, intermediate and final level.
    """
    function Pathway()
        Pathway(Level(), Level(), Level(), 0., 0., 0., false, InternalConversion.Channel[] )
    end


    # `Base.show(io::IO, pathway::InternalConversion.Pathway)`  
    # ... prepares a proper printout of the variable pathway::InternalConversion.Pathway.
    function Base.show(io::IO, pathway::InternalConversion.Pathway) 
        println(io, "initialLevel:               $(pathway.initialLevel)  ")
        println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
        println(io, "finalLevel:                 $(pathway.finalLevel)  ")
        println(io, "photonEnergy                $(pathway.photonEnergy)  ") 
        println(io, "electronEnergy              $(pathway.electronEnergy)  ")
        println(io, "crossSection:               $(pathway.crossSection)  ")
        println(io, "hasChannels:                $(pathway.hasChannels)  ")
        println(io, "channels:                   $(pathway.channels)  ")
    end



    """
    `JAC.InternalConversion.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                              settings::InternalConversion.Settings; output=true)` ... to compute the multiphoton transition 
         amplitudes and all properties as requested by the given settings. A list of lines::Array{InternalConversion.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                           settings::InternalConversion.Settings; output=true)
        println("")
        printstyled("JAC.InternalConversion.computeLines(): The computation of internal conversion amplitudes starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------- ---------------------------------------------- \n", color=:light_green)
        println("")
        #
        #
        pathways = "Not yet implemented !"
        #
        if    output    return( pathways )
        else            return( nothing )
        end
    end

end # module
