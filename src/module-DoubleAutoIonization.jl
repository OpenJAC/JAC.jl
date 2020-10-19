
"""
`module  JAC.DoubleAutoIonization`  
    ... a submodel of JAC that contains all methods for computing double Auger and autoionization amplitudes and rates.
"""
module DoubleAutoIonization

    using JAC, ..AngularMomentum, ..AtomicState, ..AutoIonization, ..Basics, ..ManyElectron, ..TableStrings


    """
    `struct  DoubleAutoIonization.Settings`  ... defines a type for the settings in estimating double-Auger and autoionization rates.

        + green                   ::Array{GreenChannel,1}   ... Precalculated and user-specified Green function of the ion.
        + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
        + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
        + minAugerEnergy          ::Float64                 ... Minimum energy of free (Auger) electrons to be included.
        + maxAugerEnergy          ::Float64                 ... Maximum energy of free (Auger) electrons to be included.
        + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
        + operator                ::AbstractEeInteraction   ... Auger operator that is to be used for evaluating the Auger amplitudes: 
                                                                allowed values are: CoulombInteraction(), BreitInteraction(), ...
        + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    """
    struct Settings
        green                     ::Array{GreenChannel,1}
        NoEnergySharings          ::Int64 
        printBefore               ::Bool 
        minAugerEnergy            ::Float64
        maxAugerEnergy            ::Float64
        maxKappa                  ::Int64  
        operator                  ::AbstractEeInteraction
        lineSelection             ::LineSelection 
    end 


    """
    `DoubleAutoIonization.Settings()`  ... constructor for the default values of DoubleAutoIonization line computations
    """
    function Settings()
        Settings(GreenChannel[], 0, false, 0., 10e5, 100, CoulombInteraction(), LineSelection())
    end


    # `Base.show(io::IO, settings::DoubleAutoIonization.Settings)`  ... prepares a proper printout of the variable settings::DoubleAutoIonization.Settings.
    function Base.show(io::IO, settings::DoubleAutoIonization.Settings) 
        println(io, "green:                         (settings.green)  ")
        println(io, "NoEnergySharings:             $(settings.NoEnergySharings)  ")
        println(io, "printBefore:                  $(settings.printBefore)  ")
        println(io, "minAugerEnergy:               $(settings.minAugerEnergy)  ")
        println(io, "maxAugerEnergy:               $(settings.maxAugerEnergy)  ")
        println(io, "maxKappa:                     $(settings.maxKappa)  ")
        println(io, "operator:                     $(settings.operator)  ")
        println(io, "lineSelection:                $(settings.lineSelection)  ")
    end


    """
    `struct  DoubleAutoIonization.ReducedChannel`  
        ... defines a type for a DoubleAutoIonization (reduced) channel to help characterize a scattering (continuum) state of many 
            electron-states with two free electrons.

        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + kappa1         ::Int64                ... partial-wave of the free electron_1
        + energy1        ::Float64              ... energy of the partial wave_1
        + phase1         ::Float64              ... phase of the partial wave_1
        + kappa2         ::Int64                ... partial-wave of the free electron_2
        + energy2        ::Float64              ... energy of the partial wave_2
        + phase2         ::Float64              ... phase of the partial wave_2
        + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
    """
    struct  ReducedChannel
        symmetry         ::LevelSymmetry
        kappa1           ::Int64 
        energy1          ::Float64 
        phase1           ::Float64
        kappa2           ::Int64   
        energy2          ::Float64 
        phase2           ::Float64 
        amplitude        ::Complex{Float64}
    end


    # `Base.show(io::IO, channel::DoubleAutoIonization.ReducedChannel)`  
    #       ... prepares a proper printout of the variable cannel::DoubleAutoIonization.ReducedChannel.
    function Base.show(io::IO, channel::DoubleAutoIonization.ReducedChannel) 
        sa = "reduced DA channel for J^P=(channel.symmetry)  with " * 
             "kappa1=$(channel.kappa1), energy1=$(channel.energy1), phase1=$(channel.phase1), " *
             "kappa2=$(channel.kappa2), energy2=$(channel.energy2), phase2=$(channel.phase2), amp=$(channel.amplitude) "
        println(io, sa)
    end


    """
    `struct  DoubleAutoIonization.Sharing`  
        ... defines a type for a DoubleAutoIonization sharing to help characterize energy sharing between the two emitted electrons.

        + epsilon1       ::Float64         ... Energy of (free) electron 1.
        + epsilon2       ::Float64         ... Energy of (free) electron 2.
        + differentialCs ::EmProperty      ... differential cross section of this energy sharing.
        + channels       ::Array{DoubleAutoIonization.ReducedChannel,1}  ... List of DoubleAutoIonization (reduced) channels of this line.
    """
    struct  Sharing
        epsilon1         ::Float64
        epsilon2         ::Float64
        differentialCs   ::EmProperty
        channels         ::Array{DoubleAutoIonization.ReducedChannel,1}
    end


    # `Base.show(io::IO, sharing::DoubleAutoIonization.Sharing)`  ... prepares a proper printout of the variable sharing::DoubleAutoIonization.Sharing.
    function Base.show(io::IO, sharing::DoubleAutoIonization.Sharing) 
        println(io, "epsilon1:               $(sharing.epsilon1)  ")
        println(io, "epsilon2:               $(sharing.epsilon2)  ")
        println(io, "differentialCs:         $(sharing.differentialCs)  ")
        println(io, "channels:               $(sharing.channels)  ")
    end


    """
    `struct  DoubleAutoIonization.Line`  ... defines a type for a double Auger line that includes sharings and their reduced amplitudes.

        + initialLevel   ::Level            ... initial-(state) level
        + finalLevel     ::Level            ... final-(state) level
        + totalRate      ::EmProperty       ... Total rate of this line.
        + sharings       ::Array{DoubleAutoIonization.Sharing,1}  ... List of DoubleAutoIonization sharings of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        totalRate        ::EmProperty
        sharings         ::Array{DoubleAutoIonization.Sharing,1}
    end 


    # `Base.show(io::IO, line::DoubleAutoIonization.Line)`  ... prepares a proper printout of the variable line::DoubleAutoIonization.Line.
     function Base.show(io::IO, line::DoubleAutoIonization.Line) 
        println(io, "initialLevel:           $(line.initialLevel)  ")
        println(io, "finalLevel:             $(line.finalLevel)  ")
        println(io, "totalRate:              $(line.totalRate)  ")
        println(io, "sharings:               $(line.sharings)  ")
    end


    """
    `DoubleAutoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                       settings::DoubleAutoIonization.Settings; output=true)`  
        ... to compute the double-Auger transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{DoubleAutoIonization.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                           settings::DoubleAutoIonization.Settings; output=true)
        println("")
        printstyled("DoubleAutoIonization.computeLines(): The computation of radiative Auger rates starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = DoubleAutoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    DoubleAutoIonization.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = DoubleAutoIonization.Line[]
        for  line in lines
            newLine = DoubleAutoIonization.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        DoubleAutoIonization.displayTotalRates(lines)
        DoubleAutoIonization.displayDifferentialRates(lines)
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `DoubleAutoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::DoubleAutoIonization.Settings)` 
        ... to determine a list of DoubleAutoIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
            and  by taking into account the particular selections and settings for this computation; an Array{DoubleAutoIonization.Line,1} is 
            returned. Apart from the level specification and sharing, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::DoubleAutoIonization.Settings)
        lines = DoubleAutoIonization.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    energy     = iLevel.energy - fLevel.energy
                    if  energy < settings.minAugerEnergy  ||  energy > settings.maxAugerEnergy    continue   end  
                    channels = DoubleAutoIonization.determineSharingsAndChannels(fLevel, iLevel, energy, settings) 
                    push!( lines, DoubleAutoIonization.Line(iLevel, fLevel, EmProperty(0., 0.,), channels) )
                end
            end
        end
        return( lines )
    end

end # module
