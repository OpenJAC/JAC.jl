
"""
`module  JAC.PhotoDoubleIonization`  
    ... a submodel of JAC that contains all methods for computing double Auger and autoionization amplitudes and rates.
"""
module PhotoDoubleIonization

    using JAC, ..AngularMomentum, ..AtomicState, ..Basics, ..ManyElectron, ..PhotoIonization, ..TableStrings


    """
    `struct  PhotoDoubleIonization.Settings`  ... defines a type for the settings in estimating single-photon double-ionization cross sections.

        + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
        + photonEnergies          ::Array{Float64,1}        ... List of photon energies.  
        + green                   ::Array{GreenChannel,1}   ... Precalculated and user-specified Green function of the ion.
        + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
        + calcAnisotropy          ::Bool                    ... True, if the beta anisotropy parameters are to be calculated and false o/w. 
        + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
        + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    """
    struct Settings
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge} 
        photonEnergies            ::Array{Float64,1}
        green                     ::Array{GreenChannel,1}
        NoEnergySharings          ::Int64 
        calcAnisotropy            ::Bool 
        printBefore               ::Bool
        lineSelection             ::LineSelection 
    end 


    """
    `PhotoDoubleIonization()`  ... constructor for the default values of PhotoDoubleIonization line computations
    """
    function Settings()
        Settings(EmMultipole[E1], UseGauge[UseCoulomb, UseBabushkin], Float64[], GreenChannel[], 0, false, false, LineSelection())
    end


    # `Base.show(io::IO, settings::PhotoDoubleIonization.Settings)`  ... prepares a proper printout of the variable settings::PhotoDoubleIonization.Settings.
    function Base.show(io::IO, settings::PhotoDoubleIonization.Settings) 
        println(io, "multipoles:                   $(settings.multipoles)  ")
        println(io, "gauges:                       $(settings.gauges)  ")
        println(io, "photonEnergies:               $(settings.photonEnergies)  ")
        println(io, "green:                         (settings.green)  ")
        println(io, "NoEnergySharings:             $(settings.NoEnergySharings)  ")
        println(io, "calcAnisotropy:               $(settings.calcAnisotropy)  ")
        println(io, "printBefore:                  $(settings.printBefore)  ")
        println(io, "lineSelection:                $(settings.lineSelection)  ")
    end


    """
    `struct  PhotoDoubleIonization.ReducedChannel`  
        ... defines a type for a PhotoDoubleIonization (reduced) channel to help characterize a scattering (continuum) state of many 
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


    # `Base.show(io::IO, channel::PhotoDoubleIonization.ReducedChannel)`  
    #       ... prepares a proper printout of the variable cannel::PhotoDoubleIonization.ReducedChannel.
    function Base.show(io::IO, channel::PhotoDoubleIonization.ReducedChannel) 
        sa = "reduced DA channel for J^P=(channel.symmetry)  with " * 
             "kappa1=$(channel.kappa1), energy1=$(channel.energy1), phase1=$(channel.phase1), " *
             "kappa2=$(channel.kappa2), energy2=$(channel.energy2), phase2=$(channel.phase2), amp=$(channel.amplitude) "
        println(io, sa)
    end


    """
    `struct  PhotoDoubleIonization.Sharing`  
        ... defines a type for a PhotoDoubleIonization sharing to help characterize energy sharing between the two emitted electrons.

        + epsilon1       ::Float64         ... Energy of (free) electron 1.
        + epsilon2       ::Float64         ... Energy of (free) electron 2.
        + differentialCs ::EmProperty      ... differential cross section of this energy sharing.
        + channels       ::Array{PhotoDoubleIonization.ReducedChannel,1}  ... List of PhotoDoubleIonization (reduced) channels of this line.
    """
    struct  Sharing
        epsilon1         ::Float64
        epsilon2         ::Float64
        differentialCs   ::EmProperty
        channels         ::Array{PhotoDoubleIonization.ReducedChannel,1}
    end


    # `Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing)`  ... prepares a proper printout of the variable sharing::PhotoDoubleIonization.Sharing.
    function Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing) 
        println(io, "epsilon1:               $(sharing.epsilon1)  ")
        println(io, "epsilon2:               $(sharing.epsilon2)  ")
        println(io, "differentialCs:         $(sharing.differentialCs)  ")
        println(io, "channels:               $(sharing.channels)  ")
    end


    """
    `struct  PhotoDoubleIonization.Line`  ... defines a type for a double Auger line that includes sharings and their reduced amplitudes.

        + initialLevel   ::Level            ... initial-(state) level
        + finalLevel     ::Level            ... final-(state) level
        + totalRate      ::EmProperty       ... Total rate of this line.
        + sharings       ::Array{PhotoDoubleIonization.Sharing,1}  ... List of PhotoDoubleIonization sharings of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        totalRate        ::EmProperty
        sharings         ::Array{PhotoDoubleIonization.Sharing,1}
    end 


    # `Base.show(io::IO, line::PhotoDoubleIonization.Line)`  ... prepares a proper printout of the variable line::PhotoDoubleIonization.Line.
     function Base.show(io::IO, line::PhotoDoubleIonization.Line) 
        println(io, "initialLevel:           $(line.initialLevel)  ")
        println(io, "finalLevel:             $(line.finalLevel)  ")
        println(io, "totalRate:              $(line.totalRate)  ")
        println(io, "sharings:               $(line.sharings)  ")
    end


    """
    `PhotoDoubleIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                       settings::PhotoDoubleIonization.Settings; output=true)`  
        ... to compute the double-Auger transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{PhotoDoubleIonization.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                           settings::PhotoDoubleIonization.Settings; output=true)
        println("")
        printstyled("PhotoDoubleIonization.computeLines(): The computation of single-photon double-ionization cs starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = PhotoDoubleIonization.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    PhotoDoubleIonization.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = PhotoDoubleIonization.Line[]
        for  line in lines
            newLine = PhotoDoubleIonization.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        PhotoDoubleIonization.displayTotalRates(lines)
        PhotoDoubleIonization.displayDifferentialRates(lines)
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `PhotoDoubleIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoDoubleIonization.Settings)` 
        ... to determine a list of PhotoDoubleIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
            and  by taking into account the particular selections and settings for this computation; an Array{PhotoDoubleIonization.Line,1} is 
            returned. Apart from the level specification and sharing, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoDoubleIonization.Settings)
        lines = PhotoDoubleIonization.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    energy   = iLevel.energy - fLevel.energy
                    channels = PhotoDoubleIonization.determineSharingsAndChannels(fLevel, iLevel, energy, settings) 
                    push!( lines, PhotoDoubleIonization.Line(iLevel, fLevel, EmProperty(0., 0.,), channels) )
                end
            end
        end
        return( lines )
    end

end # module
