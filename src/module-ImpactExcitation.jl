
"""
`module  JAC.ImpactExcitation`  
    ... a submodel of JAC that contains all methods for computing electron impact excitation cross sections and collision strengths.
"""
module ImpactExcitation 

    using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..Nuclear, ..TableStrings

    """
    `struct  ImpactExcitationSettings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing electron-impact excitation lines.

        + electronEnergies        ::Array{Float64,1}             ... List of impact-energies of the incoming elecgtrons.
        + includeBreit            ::Bool                         ... True if the Breit interaction is to be included, and false otherwise.
        + calcCollisionStrength   ::Bool                         ... True, if collision strength need to be calculated, and false otherwise.
        + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.
        + maxKappa                ::Int64                        ... Maximum kappa value of partial waves to be included.
        + energyShift             ::Float64                      ... An overall energy shift for all transitions |i> --> |f>.
    """
    struct Settings  <:  AbstractProcessSettings
        electronEnergies          ::Array{Float64,1}
        includeBreit              ::Bool 
        calcCollisionStrength     ::Bool
        printBefore               ::Bool 
        lineSelection             ::LineSelection  
        maxKappa                  ::Int64
        energyShift               ::Float64    
    end 


    """
    `ImpactExcitation.Settings()`  ... constructor for the default values of electron-impact excitation line computations.
    """
    function Settings()
       Settings( Float64[], false, false, false, LineSelection(), 0, 0.)
    end


    # `Base.show(io::IO, settings::ImpactExcitation.Settings)`  ... prepares a proper printout of the variable settings::ImpactExcitation.Settings.
    function Base.show(io::IO, settings::ImpactExcitation.Settings) 
        println(io, "electronEnergies:           $(settings.electronEnergies)  ")
        println(io, "includeBreit:               $(settings.includeBreit)  ")
        println(io, "calcCollisionStrength:      $(settings.calcCollisionStrength)  ")
        println(io, "printBefore:                $(settings.printBefore)  ")
        println(io, "lineSelection:              $(settings.lineSelection)  ")
        println(io, "maxKappa:                   $(settings.maxKappa)  ")
        println(io, "energyShift:                $(settings.energyShift)  ")
    end


    """
    `struct  ImpactExcitation.Channel`  
        ... defines a type for a electron-impact excitaiton channel to help characterize the incoming and outgoing (continuum) states of 
            many electron-states with a single free electron

        + initialKappa     ::Int64              ... partial-wave of the incoming free electron
        + finalKappa       ::Int64              ... partial-wave of the outgoing free electron
        + symmetry         ::LevelSymmetry      ... total angular momentum and parity of the scattering state
        + initialPhase     ::Float64            ... phase of the incoming partial wave
        + finalPhase       ::Float64            ... phase of the outgoing partial wave
        + amplitude        ::Complex{Float64}   ... Collision amplitude associated with the given channel.
    """
    struct  Channel
        initialKappa       ::Int64 
        finalKappa         ::Int64 
        symmetry           ::LevelSymmetry
        initialPhase       ::Float64
        finalPhase         ::Float64
        amplitude          ::Complex{Float64}
    end


    # `Base.show(io::IO, channel::ImpactExcitation.Channel)`  ... prepares a proper printout of the variable channel::ImpactExcitation.Channel.
    function Base.show(io::IO, channel::ImpactExcitation.Channel) 
        println(io, "initialKappa:       $(channel.initialKappa)  ")
        println(io, "finalKappa:         $(channel.finalKappa)  ")
        println(io, "symmetry:           $(channel.symmetry)  ")
        println(io, "initialPhase:       $(channel.initialPhase)  ")
        println(io, "finalPhase:         $(channel.finalPhase)  ")
        println(io, "amplitude:          $(channel.amplitude)  ")
    end


    """
    `struct  ImpactExcitation.Line`  
        ... defines a type for a electron-impact excitation line that may include the definition of channels and their corresponding amplitudes.

        + initialLevel           ::Level         ... initial- (bound-state) level
        + finalLevel             ::Level         ... final- (bound-state) level
        + initialElectronEnergy  ::Float64       ... energy of the incoming (initial-state) free-electron
        + finalElectronEnergy    ::Float64       ... energy of the outgoing (final-state) free-electron
        + crossSection           ::Float64       ... total cross section of this line
        + hasChannels            ::Bool          ... Determines whether the individual scattering (sub-) channels are defined in terms of their 
                                                     kappa's, phases and the total angular momentum/parity as well as the amplitude, or not.
        + channels               ::Array{ImpactExcitation.Channel,1}  ... List of ImpactExcitation channels of this line.
    """
    struct  Line
        initialLevel             ::Level
        finalLevel               ::Level
        initialElectronEnergy    ::Float64
        finalElectronEnergy      ::Float64
        crossSection             ::Float64 
        hasChannels              ::Bool
        channels                 ::Array{ImpactExcitation.Channel,1}
    end 


    """
    `ImpactExcitation.Line()`  ... 'empty' constructor for an electron-impact excitation line between a specified initial and final level.
    """
    function Line()
        Line(Level(), Level(), 0., 0., 0., false, ImpactExcitation.Channel[] )
    end


    """
    `ImpactExcitation.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`  
        ... constructor for an electron-impact excitation line between a specified initial and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)
        Line(initialLevel, finalLevel, 0., 0., crossSection, false, ImpactExcitation.Channel[] )
    end


    # `Base.show(io::IO, line::ImpactExcitation.Line)`  ... prepares a proper printout of the variable line::ImpactExcitation.Line.
    function Base.show(io::IO, line::ImpactExcitation.Line) 
        println(io, "initialLevel:            $(line.initialLevel)  ")
        println(io, "finalLevel:              $(line.finalLevel)  ")
        println(io, "initialElectronEnergy:   $(line.initialElectronEnergy)  ")
        println(io, "finalElectronEnergy:     $(line.finalElectronEnergy)  ")
        println(io, "crossSection:            $(line.crossSection)  ")
        println(io, "hasChannels:             $(line.hasChannels)  ")
        println(io, "channels:                $(line.channels)  ")
    end


    """
    `ImpactExcitation.computeAmplitudesProperties(line::ImpactExcitation.Line, grid::Radial.Grid, settings::ImpactExcitation.Settings)`  
        ... to compute all amplitudes and properties of the given line; a line::ImpactExcitation.Line is returned for which the amplitudes and
            properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::ImpactExcitation.Line, grid::Radial.Grid, settings::ImpactExcitation.Settings)
        newChannels = ImpactExcitation.Channel[]
        for channel in line.channels
            # Generate two continuum orbitals
            println("ImpactExcitation.computeAmplitudesProperties-aa: warning ... no coninuum orbitals are generated.")
            initialPhase = 0.;    finalPhase = 0.
            # Define a proper continuum basis from the initialLevel.basis and finalLevel.basis with the two continuum orbital
            println("ImpactExcitation.computeAmplitudesProperties-ab: warning ... no continuum bases are generated.")
            # Compute the transition matrix for the two constructed continuum bases
            println("ImpactExcitation.computeAmplitudesProperties-ac: warning ... no transition matrix is computed.")
            # matrix    = ImpactExcitation.computeMatrix(channel.multipole, channel.gauge, line.omega, line.finalLevel.basis, 
            #                                                line.initialLevel.basis, grid, settings)
            # amplitude = line.finalLevel.mc * matrix * line.initialLevel.mc 
            amplitude = 1.0 
            push!( newChannels, ImpactExcitation.Channel( channel.initialKappa, channel.finalKappa, channel.symmetry, 
                                                          initialPhase, finalPhase, amplitude) )
        end
        # Calculate the electron-impact excitation strength and cross section
        println("ImpactExcitation.computeAmplitudesProperties-ba: warning ... cs and strength set to -1.")
        crossSection = -1.
        line = ImpactExcitation.Line( line.initialLevel, line.finalLevel, line.initialElectronEnergy, line.finalElectronEnergy, 
                                      crossSection, true, newChannels)
        return( line )
    end


    """
    `ImpactExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                   settings::ImpactExcitation.Settings; output=true)`  
        ... to compute the electron-impact excitation transition amplitudes and all properties as requested by the given settings. 
            A list of lines::Array{ImpactExcitation.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactExcitation.Settings; 
                           output=true)
        println("")
        printstyled("ImpactExcitation.computePathways(): The computation of electron-impact excitation cross sections starts now ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = ImpactExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    ImpactExcitation.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = ImpactExcitation.Line[]
        for  line in lines
            newLine = ImpactExcitation.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        ImpactExcitation.displayResults(lines)
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `ImpactExcitation.computeMatrix(finalBasis::Basis, initialBasis::Basis, settings::ImpactExcitation.Settings)`  
        ... to compute the transition matrix  (<finalContinuumCSF_r|| V(e-e) ||initialContinuumCSF_s>)  between the CSF_r from the 
            finalContinuumBasis and the CSF_s from the initialContinuumBasis. A (non-quadratic) matrix::Array{Float64,2} with 
            dimensions [length(finalContinuumBasis.csfs) x length(initialContinuumBasis.csfs)] is returned. Note that this transition 
            matrix is typically specific to just one Eimex channel due to the different energies, partial waves and overall symmetry 
            of the scattering states. **Not yet implemented !**
    """
    function computeMatrix(finalBasis::Basis, initialBasis::Basis, settings::ImpactExcitation.Settings)   
        error("Not yet implemented.")
    end



    """
    `ImpactExcitation.determineChannels(finalLevel::Level, initialLevel::Level, settings::ImpactExcitation.Settings)`  
        ... to determine a list of electron-impact excitation Channels for a transitions from the initial to the final level and by 
            taking into account the particular settings of for this computation; an Array{ImpactExcitation.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::ImpactExcitation.Settings)
        channels = ImpactExcitation.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        #
        for  inKappa = -(settings.maxKappa+1):settings.maxKappa
            if  inKappa == 0    continue    end
            symtList = AngularMomentum.allowedTotalSymmetries(symi, inKappa)
            for  symt in symtList
                outKappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
                for  outKappa in outKappaList
                    push!(channels, ImpactExcitation.Channel( inKappa, outKappa, symt, 0., 0.,Complex(0.)) )
                end
            end
        end
        return( channels )  
    end


    """
    `ImpactExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ImpactExcitation.Settings)`  
        ... to determine a list of ImpactExcitation.Line's for transitions between levels from the initial- and final-state multiplets, 
            and by taking into account the particular selections and settings for this computation; an Array{ImpactExcitation.Line,1} is 
            returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ImpactExcitation.Settings)
        lines = ImpactExcitation.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    for  en in settings.electronEnergies
                        initialElectronEnergy  = en
                        finalElectronEnergy    = en - (iLevel.energy - fLevel.energy) + settings.energyShift
                        if  finalElectronEnergy < 0    continue   end  
                        channels = ImpactExcitation.determineChannels(fLevel, iLevel, settings) 
                        push!( lines, ImpactExcitation.Line(iLevel, fLevel, initialElectronEnergy, finalElectronEnergy, 0., true, channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `ImpactExcitation.displayLines(lines::Array{ImpactExcitation.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all 
            selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{ImpactExcitation.Line,1})
        nx = 180
        println(" ")
        println("  Selected electron-impact ionization lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(12, "Energy e_in"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(12, "Energy e_out"; na=4);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.flushleft(57, "List of partial waves and total symmetries"; na=4)  
        sb = sb * TableStrings.flushleft(57, "partial-in [total J^P] partial-out        "; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #
        NoLines = 0   
        for  line in lines
            NoLines = NoLines + 1
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "    "
            kappaInOutSymmetryList = Tuple{Int64,Int64,LevelSymmetry}[]
            for  i in 1:length(line.channels)
                push!( kappaInOutSymmetryList, (line.channels[i].initialKappa, line.channels[i].finalKappa, line.channels[i].symmetry) ) 
            end
            wa = TableStrings.kappaKappaSymmetryTupels(95, kappaInOutSymmetryList)
            sb = sa * wa[1];    println( sb )  
            for  i = 2:length(wa)
                NoLines = NoLines + 1
                sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
            # Avoid long table printouts
            if NoLines > 100  println("\n  ImpactExcitation.displayLines():  A maximum of 100 lines are printed in this table. \n")
               break   
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `ImpactExcitation.displayResults(lines::Array{ImpactExcitation.Line,1})`  
        ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is 
            returned otherwise.
    """
    function  displayResults(lines::Array{ImpactExcitation.Line,1})
        nx = 113
        println(" ")
        println("  Electron-impact excitation cross sections:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(12, "Energy e_in"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(12, "Energy e_out"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(15, "Cross section"; na=3)      
        sb = sb * TableStrings.center(15, TableStrings.inUnits("cross section"); na=3)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "     "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection))    * "    "
            println(sa)
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

end # module
