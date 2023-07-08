
"""
`module  JAC.PhotoDoubleIonization`  
    ... a submodel of JAC that contains all methods for computing double Auger and autoionization amplitudes and rates.
"""
module PhotoDoubleIonization

    using Printf, JAC, ..AngularMomentum, ..AtomicState, ..Basics, ..ManyElectron, ..PhotoIonization, ..TableStrings


    """
    `struct  PhotoDoubleIonization.Settings  <:  AbstractProcessSettings`  
        ... defines a type for the settings in estimating single-photon double-ionization cross sections.

        + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
        + photonEnergies          ::Array{Float64,1}        ... List of photon energies.  
        + activeShells            ::Array{Shell,1}          ... Non-relativistic shells that emulate parts of the two-electron continuum.
        + csScaling               ::Float64                 ... Scaling factor to modify all differential and total cross sections.
        + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
        + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
        + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
        + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    """
    struct Settings  <:  AbstractProcessSettings
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge} 
        photonEnergies            ::Array{Float64,1}
        activeShells              ::Array{Shell,1} 
        csScaling                 ::Float64    
        NoEnergySharings          ::Int64 
        printBefore               ::Bool
        maxKappa                  ::Int64 
        lineSelection             ::LineSelection 
    end 


    """
    `PhotoDoubleIonization()`  ... constructor for the default values of PhotoDoubleIonization line computations
    """
    function Settings()
        Settings(EmMultipole[E1], UseGauge[UseCoulomb, UseBabushkin], Float64[], Shell[], 0., 0, false, false, 0, LineSelection())
    end


    # `Base.show(io::IO, settings::PhotoDoubleIonization.Settings)`  ... prepares a proper printout of settings::PhotoDoubleIonization.Settings.
    function Base.show(io::IO, settings::PhotoDoubleIonization.Settings) 
        println(io, "multipoles:                   $(settings.multipoles)  ")
        println(io, "gauges:                       $(settings.gauges)  ")
        println(io, "photonEnergies:               $(settings.photonEnergies)  ")
        println(io, "activeShells:                 $(settings.activeShells)  ")
        println(io, "csScaling:                    $(settings.csScaling)  ")
        println(io, "NoEnergySharings:             $(settings.NoEnergySharings)  ")
        println(io, "printBefore:                  $(settings.printBefore)  ")
        println(io, "lineSelection:                $(settings.lineSelection)  ")
    end


    """
    `struct  PhotoDoubleIonization.Sharing`  
        ... defines a type for a PhotoDoubleIonization sharing to help characterize energy sharing between the two emitted electrons.

        + omega          ::Float64         ... Energy of the incident photon
        + epsilon1       ::Float64         ... Energy of (free) electron 1.
        + epsilon2       ::Float64         ... Energy of (free) electron 2.
        + weight         ::Float64         ... Gauss-Lengendre weight of this sharing for energy-integrated quantities.
        + differentialCs ::EmProperty      ... differential cross section of this energy sharing.
    """
    struct  Sharing
        omega            ::Float64
        epsilon1         ::Float64
        epsilon2         ::Float64
        weight           ::Float64
        differentialCs   ::EmProperty
    end


    # `Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing)`  ... prepares a proper printout of sharing::PhotoDoubleIonization.Sharing.
    function Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing) 
        println(io, "omega:                  $(sharing.omega)  ")
        println(io, "epsilon1:               $(sharing.epsilon1)  ")
        println(io, "epsilon2:               $(sharing.epsilon2)  ")
        println(io, "weight:                 $(sharing.weight)  ")
        println(io, "differentialCs:         $(sharing.differentialCs)  ")
    end


    """
    `struct  PhotoDoubleIonization.Line`  ... defines a type for a double Auger line that includes sharings and their reduced amplitudes.

        + initialLevel   ::Level            ... initial-(state) level
        + finalLevel     ::Level            ... final-(state) level
        + omega          ::Float64          ... photon energy of this line
        + totalCs        ::EmProperty       ... Total rate of this line.
        + sharings       ::Array{PhotoDoubleIonization.Sharing,1}  ... List of PhotoDoubleIonization sharings of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        omega            ::Float64
        totalCs          ::EmProperty
        sharings         ::Array{PhotoDoubleIonization.Sharing,1}
    end 


    # `Base.show(io::IO, line::PhotoDoubleIonization.Line)`  ... prepares a proper printout of line::PhotoDoubleIonization.Line.
     function Base.show(io::IO, line::PhotoDoubleIonization.Line) 
        println(io, "initialLevel:           $(line.initialLevel)  ")
        println(io, "finalLevel:             $(line.finalLevel)  ")
        println(io, "omega:                  $(line.omega)  ")
        println(io, "totalCs:                $(line.totalCs)  ")
        println(io, "sharings:               $(line.sharings)  ")
    end


    """
    `PhotoDoubleIonization.computeContinuumOrbitals(finalLevel::Level, kappas::Array{Int64,1}, epsilons::Array{Float6464,1}, 
                                                    nm::Nuclear.Model, grid::Radial.Grid)` 
        ... computes in turn all the necessary continuum orbitals for the photo-double ionization computations; a
            cOrbitals::Dict{String, Orbital} is returned, where the string is formed by string(kappa, " ", epsilon).
    """
    function computeContinuumOrbitals(finalLevel::Level, kappas::Array{Int64,1}, epsilons::Array{Float64,1}, nm::Nuclear.Model, grid::Radial.Grid)
        cOrbitals = Dict{String, Orbital}()
        #        
        # Generate potential for continuum orbitals for this step
        nrContinuum  = Continuum.gridConsistency(maximum(epsilons), grid)
        contSettings = Continuum.Settings(false, nrContinuum);   
        npot         = Nuclear.nuclearPotential(nm, grid)
        wp           = Basics.computePotentialDFS(grid, finalLevel)
        pot          = Basics.add(npot, wp)
        #  Generate continuum orbitals
        for  (ie, epsilon)  in  enumerate(epsilons)
            for  kappa in kappas    
                sh             = Subshell(101, kappa)
                cOrbital, phase, normF  = Continuum.generateOrbitalLocalPotential(epsilon, sh, pot, contSettings)
                key            = string(kappa, "  ", epsilon)
                cOrbitals[key] = cOrbital
                println(">> New continum orbital generated for $sh and energy $epsilon ")
            end
        end
        
        return( cOrbitals )
    end
     


    """
    `PhotoDoubleIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                        settings::PhotoDoubleIonization.Settings; output=true)`  
        ... to compute the -- differential and total --- double-photoionization cross sections due to the given settings. 
            A list of lines::Array{PhotoDoubleIonization.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
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
            newLine = PhotoDoubleIonization.computeLineCrossSections(line, nm, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        PhotoDoubleIonization.displayTotalCrossSections(stdout, newLines, settings)
        PhotoDoubleIonization.displayDifferentialCrossSections(stdout, newLines, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   PhotoDoubleIonization.displayTotalCrossSections(iostream, newLines, settings)
                           PhotoDoubleIonization.displayDifferentialCrossSections(iostream, newLines, settings)     end
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end


    """
    `PhotoDoubleIonization.computeLineCrossSections(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                    settings::PhotoDoubleIonization.Settings)` 
        ... to compute the differential and total cross section for the line and all its sharings. The procedure applies an 
            effective approach to treat the double continuum by means of selected active electrons with shifted resonance
            conditions. The procedure "emulates" the modulus-square of the transition amplitudes
            
                | <(alpha_f J_f, epsilon kappa) J_t || A || alpha_i J_i> |^2     with second-order operator
        
                A = V^(e-e)  O^(photoionization) + O^(photoionization)  V^(e-e)
                  = Sum_n  [ V^(e-e)  |n >< n|  O^(photoionization)  +  O^(photoionization)  |n >< n|  V^(e-e) ] / (E_i + w - E_n)
                  
            due to the electron-photon interaction V^(e-e) and the photoionization amplitude for the given final and initial level.
            
            A newLine::PhotoDoubleIonization.Line is returned for which all missing (differential) cross sections are 
            calculated.         
    """
    function  computeLineCrossSections(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                       settings::PhotoDoubleIonization.Settings)
        totalCs   = Basics.EmProperty(0.);    newSharings   = PhotoDoubleIonization.Sharing[];    asfSettings = AsfSettings()   
        conf1List = Configuration[];          multiplet1List = Multiplet[]
        conf2List = Configuration[];          multiplet2List = Multiplet[]
        #
        # Generate the list of configuration that are used in these computations
        finalConfs      = Basics.extractNonrelativisticConfigurations(line.finalLevel.basis)
        for  sha  in  settings.activeShells
            addConfs = Basics.generateConfigurationsWithAdditionalElectrons(finalConfs, [sha])
            append!(conf1List, addConfs)
        end
        conf1List = unique(conf1List)
        #
        for  sha  in  settings.activeShells
            for  shb  in  settings.activeShells
                addConfs = Basics.generateConfigurationsWithAdditionalElectrons(finalConfs, [sha, shb])
                append!(conf2List, addConfs)
            end
        end
        conf2List = unique(conf2List)
        
        # Generate all associated multiplets with these "excited" configurations
        for  conf  in conf1List
            sa = "  Multiplet computations for $(string(conf)[1:end])   with $(conf.NoElectrons) electrons ... ";   print(sa)
            basis     = Basics.performSCF([conf], nm, grid, asfSettings; printout=false)
            multiplet = Basics.performCI(basis,   nm, grid, asfSettings; printout=false)
            push!(multiplet1List, multiplet)
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
        for  conf  in conf2List
            sa = "  Multiplet computations for $(string(conf)[1:end])   with $(conf.NoElectrons) electrons ... ";   print(sa)
            basis     = Basics.performSCF([conf], nm, grid, asfSettings; printout=false)
            multiplet = Basics.performCI(basis,   nm, grid, asfSettings; printout=false)
            push!(multiplet2List, multiplet)
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
        
        # Generate all continuum orbitals in the potential of the final level; use cOrbitals with string(kappa, "  ", epsilon2)
        kappas = Int64[];   epsilons = Float64[]
        for  kappa = 1:settings.maxKappa    push!(kappas, -kappa);       push!(kappas, kappa)   end
        for  sharing  in  line.sharings     push!(epsilons, sharing.epsilon2)                   end
        epsilons  = unique(epsilons)   
        cOrbitals = PhotoDoubleIonization.computeContinuumOrbitals(line.finalLevel, kappas, epsilons, nm, grid)        
        
        # Compute the energy-differential cross sections for all sharings by adding the contributions from the various
        # configurations with two and one additional electrons; the iteration is done over all pairs of multiplets, while
        # the proper selection rules are applied based on the given configurations
        for  sharing  in  line.sharings
            differentialCs  = Basics.EmProperty(0.)
            for  multiplet2 in multiplet2List
                for  multiplet1 in multiplet1List
                    for (k,corb)  in  cOrbitals
                        differentialCs = differentialCs + PhotoDoubleIonization.computeProductOfMatrixElements(line.finalLevel, 
                                                                                                      multiplet2, multiplet1, corb, grid)
                    end 
                end 
            end 
            #
            totalCs    = totalCs + sharing.weight * differentialCs
            newSharing = PhotoDoubleIonization.Sharing(sharing.omega, sharing.epsilon1, sharing.epsilon2, sharing.weight, differentialCs)
            push!(newSharings, newSharing)
        end
  
        newLine = PhotoDoubleIonization.Line(line.initialLevel, line.finalLevel, line.omega, totalCs, newSharings )
        
        return( newLine )
    end


    """
    `PhotoDoubleIonization.computeProductOfMatrixElements(finalLevel::Level, multiplet2::Multiplet, multiplet1::Multiplet, 
                                                          cOrbital::Orbital, grid::Radial.Grid)` 
        ... to compute a single term to the differential cross section by considering just one multiplet with two additional
            electrons (compared to the finalLevel) as well as with one additional bound electron and the continuum orbital
            (cOrbital). This term is non-zero for all levels with:
            
                1) a non-zero V^(e-e) interaction matrix element to some level from multiplet2 and
                2) a non-zero O^(E1) dipole matrix element between levels of multiplet2  and  (multiplet1, cOrbital).
            
            A product::EmProperty is returned that is formed by Sum |V^(e-e)|^2 * |O^(E1)|^2, and where the summation runs formally 
            over all levels of multiplet2  and  (multiplet1, cOrbital), and where a proper coupling of the continuum orbital is 
            considered.
            
            A second term due to the reversed treatment of multiplet2  and  (multiplet1, cOrbital) is neglected as a direct transition
            from finalLevel --> (multiplet1, cOrbital) is not possible.
            
            To recognize the "zeros" of the product terms due to the shell occupation of multiplet2  and  (multiplet1, cOrbital),
            the associated configurations are extracted the compared with each other. A zero value is returned if this can be 
            extracted from the configurations involved.
    """
    function  computeProductOfMatrixElements(finalLevel::Level, multiplet2::Multiplet, multiplet1::Multiplet, 
                                             cOrbital::Orbital, grid::Radial.Grid)
        product = Basics.EmProperty(0.)
        # Analyze the configurations of finalLevel, multiplet2 and multiplet1 and compare for shell occupation and parity selection
        # (if combined with cOrbital)
        finalConfs = Basics.extractNonrelativisticConfigurations(finalLevel.basis)
        mult2Confs = Basics.extractNonrelativisticConfigurations(multiplet2.levels[1].basis)
        mult1Confs = Basics.extractNonrelativisticConfigurations(multiplet1.levels[1].basis)
        if  length(finalConfs) != 1   error("stop a: More than single configuration $finalConfs")   end
        if  length(mult2Confs) != 1   error("stop a: More than single configuration $mult2Confs")   end
        if  length(mult1Confs) != 1   error("stop a: More than single configuration $mult2Confs")   end
        #
        if      Basics.determineParity(finalConfs[1]) != Basics.determineParity(mult2Confs[1])   return( product )
        elseif  Basics.determineParity(mult2Confs[1]) * Basics.Parity("-") != 
                Basics.determineParity(mult1Confs[1]) * Basics.determineParity(cOrbital)         return( product )
        end
        
        println("Code need to be continued here ... ")
        
        
        return( product )
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
                    energy   = fLevel.energy - iLevel.energy
                    for  omega in settings.photonEnergies
                        channels = PhotoDoubleIonization.determineSharings(fLevel, iLevel, omega, energy, settings) 
                        push!( lines, PhotoDoubleIonization.Line(iLevel, fLevel, omega, EmProperty(0., 0.,), channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `PhotoDoubleIonization.determineSharings(finalLevel::Level, initialLevel::Level, omega::Float64, energy::Float64,
                                                        settings::PhotoDoubleIonization.Settings)`  
        ... to determine a list of PhotoDoubleIonization Sharing's for a transitions from the initial to final level 
            and by taking into account the particular settings of for this computation; 
            an Array{PhotoDoubleIonization.Sharing,1} is returned.
    """
    function determineSharings(finalLevel::Level, initialLevel::Level, omega::Float64, energy::Float64,
                               settings::PhotoDoubleIonization.Settings)
        sharings  = PhotoDoubleIonization.Sharing[]
        eSharings = Basics.determineEnergySharings(omega - energy, settings.NoEnergySharings) 
        for  es in eSharings
            epsilon1  = es[1];    epsilon2 = es[2];    weight = es[3] 
            push!(sharings, PhotoDoubleIonization.Sharing(omega, epsilon1, epsilon2, weight, EmProperty(0., 0.)) )
        end
        return( sharings )  
    end
    

    """
    `PhotoDoubleIonization.displayDifferentialCrossSections(stream::IO, lines::Array{PhotoDoubleIonization.Line,1},
                                                            settings::PhotoDoubleIonization.Settings)`  
        ... to display all differential rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayDifferentialCrossSections(stream::IO, lines::Array{PhotoDoubleIonization.Line,1}, 
                                               settings::PhotoDoubleIonization.Settings)
        #
        # First, print lines and sharings
        nx = 130
        println(stream, " ")
        println(stream, "  Energy-differential cross sections of selected photo-double ionization lines:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.flushleft(38, "Energies (all in " * TableStrings.inUnits("energy") * ")"; na=3);              
        sb = sb * TableStrings.flushleft(38, "  i -- f        epsilon_1     epsilon_2"; na=3)
        sa = sa * TableStrings.center(14, "Weight"; na=0);                            sb = sb * TableStrings.hBlank(13)
        sa = sa * TableStrings.center(34, "Cou -- diff. cross section -- Bab"; na=3)      
        sb = sb * TableStrings.center(34, TableStrings.inUnits("differential cross section") * "     " * 
                                          TableStrings.inUnits("differential cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
            energy = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
            #
            for  (is, sharing)  in  enumerate(line.sharings)
                if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1))                        * "    "
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2))                        * "    "
                sb = sb * @sprintf("%.4e",                                              sharing.weight)                           * "       "
                sb = sb * @sprintf("%.5e", Defaults.convertUnits("energy-diff. cross section: from atomic", sharing.differentialCs.Coulomb))   * "   "
                sb = sb * @sprintf("%.5e", Defaults.convertUnits("energy-diff. cross section: from atomic", sharing.differentialCs.Babushkin)) * "   "
                println(stream,  sb )
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `PhotoDoubleIonization.displayLines(lines::Array{PhotoDoubleIonization.Line,1})`  
        ... to display a list of lines and sharings that have been selected due to the prior settings. A neat table 
            of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{PhotoDoubleIonization.Line,1})
        #
        # First, print lines and sharings
        nx = 94
        println(" ")
        println("  Selected photo-double ionization lines & sharings:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.flushleft(54, "Energies (all in )" * TableStrings.inUnits("energy") * ")"; na=5);              
        sb = sb * TableStrings.flushleft(54, "  i -- f         omega        epsilon_1    epsilon_2"; na=5)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
            energy = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
            #
            for  (is, sharing)  in  enumerate(line.sharings)
                if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega))    * "    "
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1)) * "   "
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2)) * "   "
                println( sb )
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `PhotoDoubleIonization.displayTotalCrossSections(stream::IO, lines::Array{PhotoDoubleIonization.Line,1},
                                                     settings::PhotoDoubleIonization.Settings)`  
        ... to display all total rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayTotalCrossSections(stream::IO, lines::Array{PhotoDoubleIonization.Line,1}, 
                                        settings::PhotoDoubleIonization.Settings)
        #
        # First, print lines and sharings
        nx = 103
        println(stream, " ")
        println(stream, "  Total rates of selected photo-double ionization lines:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=3)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(30, "Cou --   cross section   -- Bab"; na=4)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "        " * 
                                          TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
            energy = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
            #
            sb = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega))                    * "         "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.totalCs.Coulomb))   * "     "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.totalCs.Babushkin)) * "   "
            println(stream,  sb )
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

     
     
    #########################################################################################################################
    #########################################################################################################################
    #########################################################################################################################


    """
    `PhotoDoubleIonization.amplitude(kind::String, channel::PhotoDoubleIonization.ReducedCoupling, 
                                     continuumLevel::Level, green::Array{GreenChannel,1}, initialLevel::Level, grid::Radial.Grid)`  
        ... to compute the kind = (photoionization) amplitude  <(alpha_f J_f, epsilon kappa) J_t || A || alpha_i J_i>  with second-order operator
        
                A = V^(e-e)  O^(photoionization) + O^(photoionization)  V^(e-e)
                  = Sum_n  V^(e-e)  |n >< n|  O^(photoionization)  +  O^(photoionization)  |n >< n|  V^(e-e) / (E_i + w - E_n)
                  
            due to the electron-photon interaction V^(e-e) and the photoionization amplitude for the given final and initial level. 
            Of course, this amplitude depends on the partial wave of the outgoing electron as well as the given multipole and gauge. 
            A value::ComplexF64 is returned.
    """

end # module
