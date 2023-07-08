
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
        + green                   ::Array{GreenChannel,1}   ... Precalculated and user-specified Green function of the ion.
        + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
        + calcAnisotropy          ::Bool                    ... True, if the beta anisotropy parameters are to be calculated and false o/w. 
        + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
        + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
        + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    """
    struct Settings  <:  AbstractProcessSettings
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge} 
        photonEnergies            ::Array{Float64,1}
        green                     ::Array{GreenChannel,1}
        NoEnergySharings          ::Int64 
        calcAnisotropy            ::Bool 
        printBefore               ::Bool
        maxKappa                  ::Int64 
        lineSelection             ::LineSelection 
    end 


    """
    `PhotoDoubleIonization()`  ... constructor for the default values of PhotoDoubleIonization line computations
    """
    function Settings()
        Settings(EmMultipole[E1], UseGauge[UseCoulomb, UseBabushkin], Float64[], GreenChannel[], 0, false, false, 0, LineSelection())
    end


    # `Base.show(io::IO, settings::PhotoDoubleIonization.Settings)`  ... prepares a proper printout of settings::PhotoDoubleIonization.Settings.
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

        + multipole      ::EmMultipole          ... Multipole of the photon absorption.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + omega          ::Float64              ... photon energy
        + energy1        ::Float64              ... energy of the partial wave_1
        + kappa1         ::Int64                ... partial-wave of the free electron_1
        + phase1         ::Float64              ... phase of the partial wave_1
        + energy2        ::Float64              ... energy of the partial wave_2
        + kappa2         ::Int64                ... partial-wave of the free electron_2
        + phase2         ::Float64              ... phase of the partial wave_2
        + tSymmetry      ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
    """
    struct  ReducedChannel
        multipole        ::EmMultipole 
        gauge            ::EmGauge 
        omega            ::Float64 
        energy1          ::Float64 
        kappa1           ::Int64 
        phase1           ::Float64
        energy2          ::Float64 
        kappa2           ::Int64   
        phase2           ::Float64 
        tSymmetry        ::LevelSymmetry
        amplitude        ::Complex{Float64}
    end


    # `Base.show(io::IO, channel::PhotoDoubleIonization.ReducedChannel)`  
    #       ... prepares a proper printout of the variable cannel::PhotoDoubleIonization.ReducedChannel.
    function Base.show(io::IO, channel::PhotoDoubleIonization.ReducedChannel) 
        println(io, "multipole:                    $(channel.multipole)  ")
        println(io, "gauge:                        $(channel.gauge)  ")
        println(io, "omega:                        $(channel.omega)  ")
        println(io, "energy1:                      $(channel.energy1)  ")
        println(io, "kappa1:                       $(channel.kappa1)  ")
        println(io, "phase1:                       $(channel.phase1)  ")
        println(io, "energy2:                      $(channel.energy2)  ")
        println(io, "kappa2:                       $(channel.kappa2)  ")
        println(io, "phase2:                       $(channel.phase2)  ")
        println(io, "tSymmetry:                    $(channel.tSymmetry)  ")
        println(io, "amplitude:                    $(channel.amplitude)  ")
    end


    """
    `struct  PhotoDoubleIonization.ReducedCoupling`  
        ... defines a type for a PhotoDoubleIonization (reduced) coupling to help characterize a scattering (continuum) state of many 
            electron-states with two free electrons.

        + multipole      ::EmMultipole          ... Multipole of the photon absorption.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + omega          ::Float64              ... photon energy
        + kappa1         ::Int64                ... partial-wave of the free electron_1
        + phase1         ::Float64              ... phase of the partial wave_1
        + xSymmetry      ::LevelSymmetry        ... angular momentum and parity after coupling the first electron
        + kappa2         ::Int64                ... partial-wave of the free electron_2
        + phase2         ::Float64              ... phase of the partial wave_2
        + tSymmetry      ::LevelSymmetry        ... total angular momentum and parity of the scattering state
    """
    struct  ReducedCoupling
        multipole        ::EmMultipole 
        gauge            ::EmGauge 
        omega            ::Float64 
        kappa1           ::Int64 
        phase1           ::Float64
        xSymmetry        ::LevelSymmetry
        kappa2           ::Int64   
        phase2           ::Float64 
        tSymmetry        ::LevelSymmetry
    end


    # `Base.show(io::IO, channel::PhotoDoubleIonization.ReducedCoupling)`  
    #       ... prepares a proper printout of the variable cannel::PhotoDoubleIonization.ReducedCoupling.
    function Base.show(io::IO, channel::PhotoDoubleIonization.ReducedCoupling) 
        println(io, "multipole:                    $(channel.multipole)  ")
        println(io, "gauge:                        $(channel.gauge)  ")
        println(io, "omega:                        $(channel.omega)  ")
        println(io, "kappa1:                       $(channel.kappa1)  ")
        println(io, "phase1:                       $(channel.phase1)  ")
        println(io, "xSymmetry:                    $(channel.xSymmetry)  ")
        println(io, "kappa2:                       $(channel.kappa2)  ")
        println(io, "phase2:                       $(channel.phase2)  ")
        println(io, "tSymmetry:                    $(channel.tSymmetry)  ")
    end


    """
    `struct  PhotoDoubleIonization.Sharing`  
        ... defines a type for a PhotoDoubleIonization sharing to help characterize energy sharing between the two emitted electrons.

        + omega          ::Float64         ... Energy of the incident photon
        + epsilon1       ::Float64         ... Energy of (free) electron 1.
        + epsilon2       ::Float64         ... Energy of (free) electron 2.
        + weight         ::Float64         ... Gauss-Lengendre weight of this sharing for energy-integrated quantities.
        + differentialCs ::EmProperty      ... differential cross section of this energy sharing.
        + channels       ::Array{PhotoDoubleIonization.ReducedChannel,1}  ... List of PhotoDoubleIonization (reduced) channels of this line.
    """
    struct  Sharing
        omega            ::Float64
        epsilon1         ::Float64
        epsilon2         ::Float64
        weight           ::Float64
        differentialCs   ::EmProperty
        channels         ::Array{PhotoDoubleIonization.ReducedChannel,1}
    end


    # `Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing)`  ... prepares a proper printout of sharing::PhotoDoubleIonization.Sharing.
    function Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing) 
        println(io, "omega:                  $(sharing.omega)  ")
        println(io, "epsilon1:               $(sharing.epsilon1)  ")
        println(io, "epsilon2:               $(sharing.epsilon2)  ")
        println(io, "weight:                 $(sharing.weight)  ")
        println(io, "differentialCs:         $(sharing.differentialCs)  ")
        println(io, "channels:               $(sharing.channels)  ")
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
    `PhotoDoubleIonization.amplitude(kind::String, channel::PhotoDoubleIonization.ReducedCoupling, 
                                     continuumLevel::Level, green::Array{GreenChannel,1}, initialLevel::Level, grid::Radial.Grid)`  
        ... to compute the kind = (photoionization) amplitude  <(alpha_f J_f, epsilon kappa) J_t || A || alpha_i J_i>  with second-order operator
        
                A = V^(e-e)  O^(photoionization) + O^(photoionization)  V^(e-e)
                  = Sum_n  V^(e-e)  |n >< n|  O^(photoionization)  +  O^(photoionization)  |n >< n|  V^(e-e) / (E_i + w - E_n)
                  
            due to the electron-photon interaction V^(e-e) and the photoionization amplitude for the given final and initial level. 
            Of course, this amplitude depends on the partial wave of the outgoing electron as well as the given multipole and gauge. 
            A value::ComplexF64 is returned.
    """
    function amplitude(kind::String, channel::PhotoDoubleIonization.ReducedCoupling,
                       continuumLevel::Level, green::Array{GreenChannel,1}, initialLevel::Level, grid::Radial.Grid)
        if      kind in [ "PDI"]
        #-----------------------
            symi = LevelSymmetry(initialLevel.J, initialLevel.parity);      symc = LevelSymmetry(continuumLevel.J, continuumLevel.parity)
            eni  = initialLevel.energy + channel.omega
            amplitude = ComplexF64(0.);
            # Compute the terms for <(alpha_f J_f, epsilon kappa) J_t || V^(e-e)  |n >< n|  O^(photoionization) || alpha_i J_i>
            for greenChannel in green
                if  AngularMomentum.isAllowedMultipole(greenChannel.symmetry, channel.multipole, symi)  &&   symc == greenChannel.symmetry
                    println(">>> Apply Green function channel with symmetry $(greenChannel.symmetry)")
                    for  (nlev, level)  in  enumerate(greenChannel. gMultiplet.levels)
                        ## if  nlev > 1  println(">>> Skip nlev > 1 ... need to be corrected !!!!");    break     end
                        piChannel = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa1, channel.xSymmetry, 
                                                            channel.phase1, ComplexF64(0.))
                        am_ni     = PhotoIonization.amplitude("photoionization", piChannel, channel.omega, level, initialLevel, grid::Radial.Grid)
                        auChannel = AutoIonization.Channel(channel.kappa2, channel.tSymmetry, channel.phase2, ComplexF64(0.))
                        am_fn     = AutoIonization.amplitude(CoulombInteraction(), auChannel, continuumLevel, level, grid; printout=false)
                        amplitude = amplitude  +  am_fn * am_ni / (eni - level.energy)
                        ## println(">> PDI first  for level $nlev  with  am_fn(A) = $am_fn  and  am_ni(P) = $am_ni")
                    end 
                else    ## println(">>> Skip Green function channel with symmetry $(greenChannel.symmetry)")
                end
            end
            #
            # Compute the terms for <(alpha_f J_f, epsilon kappa) J_t || O^(photoionization)  |n >< n|  V^(e-e) || alpha_i J_i>
            for greenChannel in green
                if  AngularMomentum.isAllowedMultipole(greenChannel.symmetry, channel.multipole, symc)  &&   symi == greenChannel.symmetry
                    println(">>> Apply Green function channel with symmetry $(greenChannel.symmetry)")
                    for  (nlev, level)  in  enumerate(greenChannel. gMultiplet.levels)
                        ## if  nlev > 1  println(">>> Skip nlev > 1 ... need to be corrected !!!!");    break     end
                        auChannel = AutoIonization.Channel(channel.kappa1, channel.xSymmetry, channel.phase1, ComplexF64(0.))
                        am_ni     = AutoIonization.amplitude(CoulombInteraction(), auChannel, level, initialLevel, grid; printout=false)
                        piChannel = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa2, channel.tSymmetry, 
                                                            channel.phase2, ComplexF64(0.))
                        am_fn     = PhotoIonization.amplitude("photoionization", piChannel, channel.omega, continuumLevel, level, grid::Radial.Grid)
                        amplitude = amplitude  +  am_fn * am_ni / (eni - level.energy)
                        ## println(">> PDI second for level $nlev  with  am_fn(P) = $am_fn  and  am_ni(A) = $am_ni")
                    end 
                else    ## println(">>> Skip Green function channel with symmetry $(greenChannel.symmetry)")
                end
            end
            #
        else    error("stop b")
        end
        
        return( amplitude )
    end


    """
    `PhotoDoubleIonization.computeAmplitudesProperties(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                                       settings::PhotoDoubleIonization.Settings)`  
        ... to compute all amplitudes and properties of the given line; a line::PhotoDoubleIonization.Line is returned 
            for which the amplitudes and properties have now been evaluated.
    """
    function  computeAmplitudesProperties(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                          settings::PhotoDoubleIonization.Settings)
        #
        symf         = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)
        newSharings  = PhotoDoubleIonization.Sharing[];   
        contSettings = Continuum.Settings(false, grid.NoPoints-50);    tcs = EmProperty(0.);     
        #
        for (ish, sharing)  in  enumerate(line.sharings)
            # Calculate only "half" of the sharings since the differential CS is always symmetric
            @show ish, sharing.epsilon1, sharing.epsilon2, sharing.weight
            if  ish > length(line.sharings)/2     continue    end
            newChannels = PhotoDoubleIonization.ReducedChannel[];    dcsC = 0.;    dcsB = 0.
            for (ich, ch)  in  enumerate(sharing.channels)
                ## @show  ich, ch.kappa1, ch.kappa2
                newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
                newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, newiLevel.basis.subshells)
                # Generate two continuum orbitals as associated with this sharing; these orbitals move each in the potential of the final level
                cOrbital1, phase1 = Continuum.generateOrbitalForLevel(ch.energy1, Subshell(101, ch.kappa1), newfLevel, nm, grid, contSettings)
                cOrbital2, phase2 = Continuum.generateOrbitalForLevel(ch.energy2, Subshell(102, ch.kappa2), newfLevel, nm, grid, contSettings)
                # Expand these continuum orbitals in the orbital set as given by the Green function
                greenOrbs  = settings.green[1].gMultiplet.levels[1].basis.orbitals
                bOrbitals1 = Basics.expandOrbital(cOrbital1, greenOrbs, 1.0e-1, grid, printout=true)
                bOrbitals2 = Basics.expandOrbital(cOrbital2, greenOrbs, 1.0e-1, grid, printout=true)
                # Loop through all sub-orbital contribution in the expansion of the continuum orbitals
                amplitude  = ComplexF64(0.)
                for  bOrb1  in  bOrbitals1
                    if   bOrb1.subshell  in  newfLevel.basis.subshells        continue    end
                    for  bOrb2  in  bOrbitals2
                        if   bOrb2.subshell  in  newfLevel.basis.subshells    continue    end
                        println(">>> Omega = $(line.omega), Sharing = $ish,  channel = $ich,  subsh1 = $(bOrb1.subshell), subsh2 = $(bOrb2.subshell)")
                        if      bOrb1.subshell < bOrb2.subshell
                            symxList = AngularMomentum.allowedDoubleKappaSymmetries(symf, bOrb1.subshell.kappa, bOrb2.subshell.kappa, ch.tSymmetry)
                            ##x newcLevel  = Basics.generateLevelWithExtraTwoElectrons(bOrb1, ch.xSymmetry, bOrb2, ch.tSymmetry, newfLevel)
                            for  symx  in  symxList
                                newCoupling = PhotoDoubleIonization.ReducedCoupling(ch.multipole, ch.gauge, ch.omega, 
                                                                                    bOrb1.subshell.kappa, phase1, symx,
                                                                                    bOrb2.subshell.kappa, phase2, ch.tSymmetry)
                                newcLevel   = Basics.generateLevelWithExtraTwoElectrons(bOrb1, symx, bOrb2, ch.tSymmetry, newfLevel)
                                amplitude   = amplitude + 
                                              PhotoDoubleIonization.amplitude("PDI", newCoupling, newcLevel, settings.green, newiLevel, grid)
                            end
                        elseif  bOrb1.subshell > bOrb2.subshell
                            symxList = AngularMomentum.allowedDoubleKappaSymmetries(symf, bOrb2.subshell.kappa, bOrb1.subshell.kappa, ch.tSymmetry)
                            for  symx  in  symxList
                                newCoupling = PhotoDoubleIonization.ReducedCoupling(ch.multipole, ch.gauge, ch.omega, 
                                                                                    bOrb2.subshell.kappa, phase2, symx,
                                                                                    bOrb1.subshell.kappa, phase1, ch.tSymmetry)
                                newcLevel  = Basics.generateLevelWithExtraTwoElectrons(bOrb2, symx, bOrb1, ch.tSymmetry, newfLevel)
                                amplitude  = amplitude - 
                                             PhotoDoubleIonization.amplitude("PDI", newCoupling, newcLevel, settings.green, newiLevel, grid)
                            end
                        elseif   bOrb1.subshell  ==  bOrb2.subshell     error("stop a")
                        end
                    end
                end
                push!( newChannels, ReducedChannel(ch.multipole, ch.gauge, ch.omega, ch.energy1, ch.kappa1, phase1, 
                                                                                     ch.energy2, ch.kappa2, phase2, ch.tSymmetry, amplitude) )
                if       ch.gauge == Basics.Coulomb     dcsC = dcsC + abs(amplitude)^2
                elseif   ch.gauge == Basics.Babushkin   dcsB = dcsB + abs(amplitude)^2
                elseif   ch.gauge == Basics.Magnetic    dcsB = dcsB + abs(amplitude)^2;   dcsC = dcsC + abs(amplitude)^2
                end
                @show ich, ch.gauge, dcsC, dcsB, "   ", amplitude
            end
            # Calculate the differential cross section 
            csFactor = 4 * pi^2 / Defaults.getDefaults("alpha") / line.omega
            diffCs   = EmProperty(csFactor*dcsC, csFactor*dcsB)
            push!( newSharings, PhotoDoubleIonization.Sharing( sharing.omega, sharing.epsilon1, sharing.epsilon2, sharing.weight, diffCs, newChannels) )
            tcs = tcs + sharing.weight * 2 * diffCs  ## Factor 2 because only half of the symmetric differential CS are calculated
        end
        # Calculate the total cross section 
        newLine = PhotoDoubleIonization.Line( line.initialLevel, line.finalLevel, line.omega, tcs, newSharings)
        return( newLine )
    end


    """
    `PhotoDoubleIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                       settings::PhotoDoubleIonization.Settings; output=true)`  
        ... to compute the double-Auger transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{PhotoDoubleIonization.Lines} is returned.
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
            newLine = PhotoDoubleIonization.computeAmplitudesProperties(line, nm, grid, settings) 
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
                        channels = PhotoDoubleIonization.determineSharingsAndChannels(fLevel, iLevel, omega, energy, settings) 
                        push!( lines, PhotoDoubleIonization.Line(iLevel, fLevel, omega, EmProperty(0., 0.,), channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `PhotoDoubleIonization.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, omega::Float64, energy::Float64,
                                                        settings::PhotoDoubleIonization.Settings)`  
        ... to determine a list of PhotoDoubleIonization Sharing's and Channel's for a transitions from the initial to final level 
            and by taking into account the particular settings of for this computation; an Array{PhotoDoubleIonization.Sharing,1} 
            is returned.
    """
    function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, omega::Float64, energy::Float64,
                                          settings::PhotoDoubleIonization.Settings)
        sharings  = PhotoDoubleIonization.Sharing[]
        eSharings = Basics.determineEnergySharings(omega - energy, settings.NoEnergySharings) 
        for  es in eSharings
            epsilon1  = es[1];    epsilon2 = es[2];    weight = es[3] 
            channels  = PhotoDoubleIonization.ReducedChannel[];   
            symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
            for  mp in settings.multipoles
                for  gauge in settings.gauges
                    symList = AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                    for  symt in symList
                        # Use the allowed values kappa_1 kappa_2 and deal with their coupling later
                        kappaPairs = AngularMomentum.allowedDoubleKappas(symf, symt, settings.maxKappa)
                        for  kappas  in  kappaPairs
                            # Include further restrictions if appropriate
                            if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      
                                push!(channels, ReducedChannel(mp, Basics.Coulomb,   omega, epsilon1, kappas[1], 0.,
                                                                                            epsilon2, kappas[2], 0., symt, Complex(0.)) ) 
                            elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                                push!(channels, ReducedChannel(mp, Basics.Babushkin, omega, epsilon1, kappas[1], 0.,
                                                                                            epsilon2, kappas[2], 0., symt, Complex(0.)) ) 
                            elseif string(mp)[1] == 'M'                                
                                push!(channels, ReducedChannel(mp, Basics.Magnetic,  omega, epsilon1, kappas[1], 0.,
                                                                                            epsilon2, kappas[2], 0., symt, Complex(0.)) ) 
                            end 
                        end  
                    end
                end
            end
            push!(sharings, PhotoDoubleIonization.Sharing(omega, epsilon1, epsilon2, weight, EmProperty(0., 0.), channels) )
        end
        return( sharings )  
    end


    """
    `PhotoDoubleIonization.displayDifferentialCrossSections(stream::IO, lines::Array{PhotoDoubleIonization.Line,1}, settings::PhotoDoubleIonization.Settings)`  
        ... to display all differential rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayDifferentialCrossSections(stream::IO, lines::Array{PhotoDoubleIonization.Line,1}, settings::PhotoDoubleIonization.Settings)
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
        ... to display a list of lines, sharings and channels that have been selected due to the prior settings. A neat table 
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
        #
        # Second, print lines and channles
        nx = 130
        println(" ")
        println("  Selected photo-double ionization lines & channels:   ... channels are shown only for the first sharing")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.flushleft(100, "Channels (all energies in " * TableStrings.inUnits("energy") * ")" ; na=5);              
        sb = sb * TableStrings.flushleft(100, "Multipole  Gauge      omega       (epsilon, kappa)_1   (epsilon, kappa)_2   -->    J^P_t"; na=5)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #
        for  line in lines
            sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            for (ic, ch) in enumerate(line.sharings[1].channels)
                if  ic == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
                skappa1 = "  " * string(ch.kappa1);         skappa2 = "  " * string(ch.kappa2);    
                stsym   = "      " * string(ch.tSymmetry)
                sb = sb * string(ch.multipole) * "         " * string(ch.gauge)[1:3] * "      " 
                sb = sb * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", ch.omega)) * "     "
                sb = sb * " (" * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", ch.energy1)) * ", " * skappa1[end-2:end] * ")    "
                sb = sb * " (" * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", ch.energy2)) * ", " * skappa2[end-2:end] * ")    -->"
                sb = sb * stsym[end-8:end]
                println( sb )
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `PhotoDoubleIonization.displayTotalCrossSections(stream::IO, lines::Array{PhotoDoubleIonization.Line,1}, settings::PhotoDoubleIonization.Settings)`  
        ... to display all total rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayTotalCrossSections(stream::IO, lines::Array{PhotoDoubleIonization.Line,1}, settings::PhotoDoubleIonization.Settings)
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

end # module
