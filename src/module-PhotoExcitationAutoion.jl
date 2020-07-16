
"""
`module  PhotoExcitationAutoion`  
    ... a submodel of JAC that contains all methods for computing photo-excitation-autoionization cross sections and rates.
"""
module PhotoExcitationAutoion 

    using Printf, ..AngularMomentum, ..AutoIonization, ..Basics, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear, ..PhotoEmission, 
                  ..Radial, ..TableStrings

    """
    `struct  PhotoExcitationAutoion.Settings`  
        ... defines a type for the details and parameters of computing photon-impact excitation-autoionization pathways 
            |i(N)>  --> |m(N)>  --> |f(N-1)>.

        + multipoles              ::Array{Basics.EmMultipole,1}    ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{Basics.UseGauge,1}       ... Specifies the gauges to be included into the computations.
        + printBefore             ::Bool                           ... True, if all energies and lines are printed before their evaluation.
        + pathwaySelection        ::PathwaySelection               ... Specifies the selected levels/pathways, if any.
        + maxKappa                ::Int64                          ... Maximum kappa value of partial waves to be included.
    """
    struct Settings
        multipoles                ::Array{Basics.EmMultipole,1}
        gauges                    ::Array{Basics.UseGauge,1} 
        printBefore               ::Bool 
        pathwaySelection          ::PathwaySelection 
        maxKappa                  ::Int64 
    end 


    """
    `PhotoExcitationAutoion.Settings()`  
        ... constructor for the default values of photon-impact excitation-autoionizaton settings.
    """
    function Settings()
        Settings( Basics.EmMultipole[], UseGauge[], false,  PathwaySelection(), 0)
    end


    # `Base.show(io::IO, settings::PhotoExcitationAutoion.Settings)`  
    # 	... prepares a proper printout of the variable settings::PhotoExcitationAutoion.Settings.  
    function Base.show(io::IO, settings::PhotoExcitationAutoion.Settings) 
        println(io, "multipoles:              $(settings.multipoles)  ")
        println(io, "gauges:                  $(settings.gauges)  ")
        println(io, "printBefore:             $(settings.printBefore)  ")
        println(io, "pathwaySelection:        $(settings.pathwaySelection)  ")
        println(io, "maxKappa:                $(settings.maxKappa)  ")
    end


    """
    `struct  PhotoExcitationAutoion.Pathway`  
        ... defines a type for a photon-impact excitation pathway that may include the definition of different 
            excitation and autoionization channels and their corresponding amplitudes.

        + initialLevel        ::Level           ... initial-(state) level
        + intermediateLevel   ::Level           ... intermediate-(state) level
        + finalLevel          ::Level           ... final-(state) level
        + excitEnergy         ::Float64         ... photon excitation energy of this pathway
        + electronEnergy      ::Float64         ... energy of the (finally outgoing, scattered) electron
        + crossSection        ::EmProperty      ... total cross section of this pathway
        + hasChannels         ::Bool            
            ... Determines whether the individual excitation and autoionization channels are defined in terms of their multipole, 
                gauge, free-electron kappa, phases and the total angular momentum/parity as well as the amplitude, or not.
        + excitChannels       ::Array{PhotoEmission.Channel,1}  ... List of excitation channels of this pathway.
        + augerChannels       ::Array{AutoIonization.Channel,1} ... List of Auger channels of this pathway.
    """
    struct  Pathway
        initialLevel          ::Level
        intermediateLevel     ::Level
        finalLevel            ::Level
        excitEnergy           ::Float64
        electronEnergy        ::Float64
        crossSection          ::EmProperty
        hasChannels           ::Bool
        excitChannels         ::Array{PhotoEmission.Channel,1}  
        augerChannels         ::Array{AutoIonization.Channel,1}
    end 


    """
    `PhotoExcitationAutoion.Pathway()`  
        ... 'empty' constructor for an photon-impact excitation-autoionization pathway between a specified initial, intermediate 
            and final level.
    """
    function Pathway()
        Pathway(Level(), Level(), Level(), 0., 0., EmProperty(0., 0.), false, PhotoEmission.Channel[], AutoIonization.Channel[] )
    end


    # `Base.show(io::IO, pathway::PhotoExcitationAutoion.Pathway)`  
    #		... prepares a proper printout of the variable pathway::PhotoExcitationAutoion.Pathway.
    function Base.show(io::IO, pathway::PhotoExcitationAutoion.Pathway) 
        println(io, "initialLevel:               $(pathway.initialLevel)  ")
        println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
        println(io, "finalLevel:                 $(pathway.finalLevel)  ")
        println(io, "excitEnergy                 $(pathway.excitEnergy)  ") 
        println(io, "electronEnergy              $(pathway.electronEnergy)  ")
        println(io, "crossSection:               $(pathway.crossSection)  ")
        println(io, "hasChannels:                $(pathway.hasChannels)  ")
        println(io, "excitEnergy:                $(pathway.excitEnergy)  ")
        println(io, "augerEnergy:                $(pathway.augerEnergy)  ")
    end



    """
    `PhotoExcitationAutoion.computeAmplitudesProperties(pathway::PhotoExcitationAutoion.Pathway, nm::Nuclear.Model, grid::Radial.Grid, 
                                                        nrContinuum::Int64, settings::PhotoExcitationAutoion.Settings)` 
        ... to compute all amplitudes and properties of the given pathway; a line::PhotoExcitationAutoion.Pathway is returned for which 
            the amplitudes and properties have now been evaluated.
    """
    function  computeAmplitudesProperties(pathway::PhotoExcitationAutoion.Pathway, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                          settings::PhotoExcitationAutoion.Settings)
        # Compute all excitation channels
        neweChannels = PhotoEmission.Channel[]
        for eChannel in pathway.excitChannels
            amplitude   = PhotoEmission.amplitude("absorption", eChannel.multipole, eChannel.gauge, pathway.excitEnergy, 
                                                  pathway.intermediateLevel, pathway.initialLevel, grid)
             push!( neweChannels, PhotoEmission.Channel( eChannel.multipole, eChannel.gauge, amplitude))
        end
        # Compute all AutoIonization decay channels
        newaChannels = AutoIonization.Channel[];   contSettings = Continuum.Settings(false, nrContinuum)
        for aChannel in pathway.augerChannels
            newnLevel   = Basics.generateLevelWithSymmetryReducedBasis(pathway.intermediateLevel)
            newnLevel   = Basics.generateLevelWithExtraSubshell(Subshell(101, aChannel.kappa), newnLevel)
            newfLevel   = Basics.generateLevelWithSymmetryReducedBasis(pathway.finalLevel)
            cOrbital, phase  = Continuum.generateOrbitalForLevel(pathway.electronEnergy, Subshell(101, aChannel.kappa), newfLevel, nm, grid, contSettings)
            newcLevel   = Basics.generateLevelWithExtraElectron(cOrbital, aChannel.symmetry, newfLevel)
            newcChannel = AutoIonization.Channel( aChannel.kappa, aChannel.symmetry, phase, Complex(0.))
            amplitude = 1.0
            ## amplitude   = AutoIonization.amplitude("Coulomb", aChannel, newnLevel, newcLevel, grid)
            push!( newaChannels, AutoIonization.Channel( aChannel.kappa, aChannel.symmetry, phase, amplitude))
        end
        #
        crossSection = EmProperty(-1., -1.)
        pathway = PhotoExcitationAutoion.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, 
                                                  pathway.excitEnergy, pathway.electronEnergy, crossSection, true, neweChannels, newaChannels)
        return( pathway )
    end



    """
    `PhotoExcitationAutoion.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                            nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoExcitation.Settings; output=true)`  
        ... to compute the photo-excitation-autoionization amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{PhotoExcitationAutoion.Lines} is returned.
    """
    function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, 
                              grid::Radial.Grid, settings::PhotoExcitationAutoion.Settings; output=true)
        println("")
        printstyled("PhotoExcitationAutoion.computePathways(): The computation of photo-excitation-autoionization amplitudes starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        pathways = PhotoExcitationAutoion.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    PhotoExcitationAutoion.displayPathways(pathways)    end
        # Determine maximum (electron) energy and check for consistency of the grid
        maxEnergy = 0.;   for  pathway in pathways   maxEnergy = max(maxEnergy, pathway.electronEnergy)   end
        nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newPathways = PhotoExcitationAutoion.Pathway[]
        for  pathway in pathways
            push!( newPathways, PhotoExcitationAutoion.computeAmplitudesProperties(pathway, nm, grid, nrContinuum, settings) )
        end
        # Print all results to screen
        PhotoExcitationAutoion.displayResults(stdout, newPathways)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    PhotoExcitationAutoion.displayResults(iostream, newPathways)   end
        #
        if    output    return( newPathways )
        else            return( nothing )
        end
    end


    """
    `PhotoExcitationAutoion.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                              settings::PhotoExcitationAutoion.Settings)`  
        ... to determine a list of photoexcitation-autoionization pathways between the levels from the given initial-, intermediate- and 
            final-state multiplets and by taking into account the particular selections and settings for this computation; an 
            Array{PhotoExcitationAutoion.Line,1} is returned. Apart from the level specification, all physical properties are set to zero 
            during the initialization process.  
    """
    function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                settings::PhotoExcitationAutoion.Settings)
        pathways = PhotoExcitationAutoion.Pathway[]
        for  iLevel  in  initialMultiplet.levels
            for  nLevel  in  intermediateMultiplet.levels
                for  fLevel  in  finalMultiplet.levels
                    if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)
                        eEnergy = nLevel.energy - iLevel.energy
                        aEnergy = nLevel.energy - fLevel.energy
                        if  eEnergy < 0.   ||   aEnergy < 0    continue    end
                        rSettings = PhotoEmission.Settings( settings.multipoles, settings.gauges, false, false, LineSelection(), 0., 0., 0.)
                        eChannels = PhotoEmission.determineChannels(nLevel, iLevel, rSettings) 
                        aSettings = AutoIonization.Settings( false, false, LineSelection(), 0., 0., settings.maxKappa, CoulombInteraction())
                        aChannels = AutoIonization.determineChannels(fLevel, nLevel, aSettings) 
                        push!( pathways, PhotoExcitationAutoion.Pathway(iLevel, nLevel, fLevel, eEnergy, aEnergy, EmProperty(0., 0.), 
                                                                        true, eChannels, aChannels) )
                    end
                end
            end
        end
        return( pathways )
    end


    """
    `PhotoExcitationAutoion.displayPathways(pathways::Array{PhotoExcitationAutoion.Line,1})`  
        ... to display a list of pathways and channels that have been selected due to the prior settings. A neat table of 
            all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayPathways(pathways::Array{PhotoExcitationAutoion.Pathway,1})
        nx = 170
        println(" ")
        println("  Selected photo-excitation-autoionization pathways:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "     ";   sb = "     "
        sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);          
        sa = sa * TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=3);
        sa = sa * TableStrings.center(14, "Energy m-i"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(14, "Energy m-f"; na=3);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
        sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  pathway in pathways
            sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=2)
            sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))   * "    "
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "   "
            kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
            for  ech in pathway.excitChannels
                for  ach in pathway.augerChannels
                    ##x eChannel = pathway.channels[i].excitationChannel;    aChannel = pathway.channels[i].augerChannel;  
                    push!( kappaMultipoleSymmetryList, (ach.kappa, ech.multipole, ech.gauge, ach.symmetry) )
                end
            end
            wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
            if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
            for  i = 2:length(wa)
                sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `PhotoExcitationAutoion.displayResults(stream::IO, pathways::Array{PhotoExcitationAutoion.Line,1})`  
        ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but 
            nothing is returned otherwise.
    """
    function  displayResults(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1})
        nx = 135
        println(stream, " ")
        println(stream, "  Partial photo-excitation & autoionization cross sections:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "     ";   sb = "     "
        sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);          
        sa = sa * TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=3);
        sa = sa * TableStrings.center(14, "Energy m-i"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(14, "Energy m-f"; na=3);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(30, "Cou -- Cross sections -- Bab"; na=2);       
        sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section")*"          "*
                                              TableStrings.inUnits("cross section"); na=2)
        println(stream, sa);    println(stream, sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  pathway in pathways
            sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=2)
            sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))    * "    "
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  ech in pathway.excitChannels
                multipoles = push!( multipoles, ech.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", pathway.crossSection.Coulomb))     * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", pathway.crossSection.Babushkin))   * "    "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

end # module
