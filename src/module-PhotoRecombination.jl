
"""
`module  JAC.PhotoRecombination`  
    ... a submodel of JAC that contains all methods for computing photo-recomnbination, i.e. radiative recombination 
        and radiative electron capture properties between some initial and final-state multiplets.
"""
module PhotoRecombination

    using Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..ManyElectron, ..PhotoEmission, ..Radial, ..Nuclear, ..TableStrings

    """
    `struct  Settings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing photo recombination lines.

        + multipoles          ::Array{EmMultipole}  ... Multipoles of the radiation field that are to be included.
        + gauges              ::Array{UseGauge}     ... Gauges to be included into the computations.
        + electronEnergies    ::Array{Float64,1}    ... List of electron energies [in default units].
        + ionEnergies         ::Array{Float64,1}    ... List of ion energies [in MeV/u].
        + useIonEnergies      ::Bool                ... Make use of ion energies in [MeV/u] to obtain the electron energies.
        + calcAnisotropy      ::Bool                ... True, if the overall anisotropy is to be calculated.
        + calcTensors         ::Bool                ... True, if the statistical tensors are to be calculated and 
                                                        false otherwise.
        + printBefore         ::Bool                ... True, if all energies and lines are printed before their evaluation.
        + lineSelection       ::LineSelection       ... Specifies the selected levels, if any.
    """
    struct Settings  <:  AbstractProcessSettings
        multipoles            ::Array{EmMultipole}
        gauges                ::Array{UseGauge}
        electronEnergies      ::Array{Float64,1} 
        ionEnergies           ::Array{Float64,1}
        useIonEnergies        ::Bool
        calcAnisotropy        ::Bool
        calcTensors           ::Bool 
        printBefore           ::Bool 
        lineSelection         ::LineSelection 
    end 


    """
    `PhotoRecombination.Settings()`  ... constructor for the default values of photo recombination line computations
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], Float64[], false, false, false, false, LineSelection() )
    end


    # `Base.show(io::IO, settings::PhotoRecombination.Settings)`  
    #		... prepares a proper printout of the variable settings::PhotoRecombination.Settings.
    function Base.show(io::IO, settings::PhotoRecombination.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "electronEnergies:         $(settings.electronEnergies)  ")
        println(io, "ionEnergies:              $(settings.ionEnergies)  ")
        println(io, "useIonEnergies:           $(settings.useIonEnergies)  ")
        println(io, "calcAnisotropy:           $(settings.calcAnisotropy)  ")
        println(io, "calcTensors:              $(settings.calcTensors)  ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "lineSelection:            $(settings.lineSelection)  ")
    end


    """
    `struct  PhotoRecombination.Channel`  
        ... defines a type for a photorecombination channel to help characterize a single multipole and scattering 
            (continuum) state of many electron-states with a single free electron.

        + multipole      ::EmMultipole          ... Multipole of the photon emission/absorption.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + kappa          ::Int64                ... partial-wave of the free electron
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + phase          ::Float64              ... phase of the partial wave
        + amplitude      ::Complex{Float64}     ... Rec amplitude associated with the given channel.
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
    `struct  PhotoRecombination.Line`  
        ... defines a type for a Photorecombination line that may include the definition of channels.

        + initialLevel   ::Level                  ... initial-(state) level
        + finalLevel     ::Level                  ... final-(state) level
        + electronEnergy ::Float64                ... Energy of the (incoming free) electron.
        + photonEnergy   ::Float64                ... Energy of the emitted photon.
        + betaGamma2     ::Float64                ... beta^2 * gamma^2.
        + crossSection   ::EmProperty             ... Cross section for this electron capture.
        + hasChannels    ::Bool                   ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                      free-electron energy, kappa, multipole, etc., or not.
        + channels       ::Array{PhotoRecombination.Channel,1}    ... List of photorecombination channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        electronEnergy   ::Float64
        photonEnergy     ::Float64
        betaGamma2       ::Float64 
        crossSection     ::EmProperty
        hasChannels      ::Bool
        channels         ::Array{PhotoRecombination.Channel,1}
    end 


    """
    `PhotoRecombination.Line()`  
        ... constructor for an `empty instance` of a photorecombination line between a specified initial 
            and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level)
        Line(initialLevel::Level, finalLevel::Level, 0., 0., 0., EmProperty(0., 0.), false, Channel[])
    end


    # `Base.show(io::IO, line::PhotoRecombination.Line)`  ... prepares a proper printout of the variable line::PhotoRecombination.Line.
    function Base.show(io::IO, line::PhotoRecombination.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "electronEnergy:    $(line.electronEnergy)  ")
        println(io, "photonEnergy:      $(line.photonEnergy)  ")
        println(io, "betaGamma2:        $(line.betaGamma2)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
        println(io, "hasChannels:       $(line.hasChannels)  ")
        println(io, "channels:          $(line.channels)  ")
    end


    """
    `PhotoRecombination.amplitude(kind::String, channel::PhotoRecombination.Channel, energy::Float64, finalLevel::Level, 
                                  continuumLevel::Level, grid::Radial.Grid)`  
        ... to compute the kind = (photorecombination) amplitude  
            < alpha_f J_f || O^(photorecombination) || (alpha_i J_i, epsilon kappa) J_t>  due to the electron-photon 
            interaction for the given final and continuum level, the partial wave of the outgoing electron as well as 
            the given multipole and gauge. A value::ComplexF64 is returned.
    """
    function amplitude(kind::String, channel::PhotoRecombination.Channel, energy::Float64, finalLevel::Level, 
                       continuumLevel::Level, grid::Radial.Grid)
        if      kind in [ "photorecombination"]
        #-----------------------------------
            amplitude = PhotoEmission.amplitude("emission", channel.multipole, channel.gauge, energy, finalLevel, continuumLevel, grid)
            amplitude = im^Basics.subshell_l(Subshell(101, channel.kappa)) * exp( im*channel.phase ) * amplitude
        else    error("stop b")
        end
        
        return( amplitude )
    end


    """
    `PhotoRecombination.checkConsistentMultiplets(finalMultiplet::Multiplet, initialMultiplet::Multiplet)`  
        ... to check that the given initial- and final-state levels and multiplets are consistent to each other and
            to avoid later problems with the computations. An error message is issued if an inconsistency occurs,
            and nothing is returned otherwise.
    """
    function  checkConsistentMultiplets(finalMultiplet::Multiplet, initialMultiplet::Multiplet)
        initialSubshells      = initialMultiplet.levels[1].basis.subshells;             ni = length(initialSubshells)
        finalSubshells        = finalMultiplet.levels[1].basis.subshells
        
        if initialSubshells[1:end] == finalSubshells[1:ni]
        else
            @show initialSubshells
            @show finalSubshells
            error("\nThe order of subshells must be equal for the initial- and final states. \n" *
                  "However, the initial states can have less subshells; this limitation arises from the angular coefficients.")
        end
            
        return( nothing )
    end


    """
    `PhotoRecombination.computeAmplitudesProperties(line::PhotoRecombination.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                    settings::PhotoRecombination.Settings)`  
        ... to compute all amplitudes and properties of the given line; a line::PhotoRecombination.Line is returned for 
            which the amplitudes and properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::PhotoRecombination.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                          settings::PhotoRecombination.Settings)
        newChannels = PhotoRecombination.Channel[];;   contSettings = Continuum.Settings(false, nrContinuum);    csC = 0.;    csB = 0.
        for channel in line.channels
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, line.finalLevel.basis.subshells)
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, newfLevel.basis.subshells)
            newfLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newfLevel)
            cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newiLevel, nm, grid, contSettings)
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newiLevel)
            newChannel = PhotoRecombination.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
            amplitude  = PhotoRecombination.amplitude("photorecombination", channel, line.photonEnergy, newfLevel, newcLevel, grid)
            push!( newChannels, PhotoRecombination.Channel(newChannel.multipole, newChannel.gauge, newChannel.kappa, newChannel.symmetry, 
                                                           newChannel.phase, amplitude) )
            if       channel.gauge == Basics.Coulomb     csC = csC + abs(amplitude)^2
            elseif   channel.gauge == Basics.Babushkin   csB = csB + abs(amplitude)^2
            elseif   channel.gauge == Basics.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
            end
        end
        Ji2 = Basics.twice(line.initialLevel.J)
        csFactor     = 8 * pi^3 * Defaults.getDefaults("alpha")^3 * line.photonEnergy / (Ji2 + 1)
        crossSection = 1.0 /line.betaGamma2 * EmProperty(csFactor * csC, csFactor * csB)
        newLine      = PhotoRecombination.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, line.betaGamma2, 
                                                crossSection, true, newChannels)
        return( newLine )
    end



    """
    `PhotoRecombination.computeAnisotropyParameter(nu::Int64, gauge::EmGauge, line::PhotoRecombination.Line)`  
        ... to compute the anisotropy parameter of the emitted photons for the photorecombination of an initially unpolarized ion. 
            A value::ComplexF64 is returned.
    """
    function  computeAnisotropyParameter(nu::Int64, gauge::EmGauge, line::PhotoRecombination.Line)
        if  !line.hasChannels   error("No channels are defined for the given PhotoRecombination.line.")         end
        wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J   
        wn = 0.;    
        for  ch in line.channels   
            if  gauge != ch.gauge  &&  gauge != Basics.Magnetic   continue    end
            wn = wn + conj(ch.amplitude) * ch.amplitude   
        end
        
        for  cha  in line.channels
            if  gauge != cha.gauge  &&  gauge != Basics.Magnetic   continue    end
            J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
            j = AngularMomentum.kappa_j(cha.kappa);    l = AngularMomentum.kappa_l(cha.kappa)
            #
            for  chp  in line.channels  
                if  gauge != chp.gauge  &&  gauge != Basics.Magnetic   continue    end
                Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                jp = AngularMomentum.kappa_j(chp.kappa);     lp = AngularMomentum.kappa_l(chp.kappa)
                #
                if  1 + (-1)^(L + p + Lp + pp - nu) == 0    continue    end
                @show 1.0im^(L + p - Lp - pp)
                @show AngularMomentum.phaseFactor([Ji, -1, AngularJ64(1//2), -1, Jf])
                @show sqrt( AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), l, lp, j, jp, J, Jp]) )
                wa = wa + 1.0im^(L + p - Lp - pp) * AngularMomentum.phaseFactor([Ji, -1, AngularJ64(1//2), -1, Jf]) *
                                 sqrt( AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), l, lp, j, jp, J, Jp]) ) *  
                                 AngularMomentum.ClebschGordan( l, AngularM64(0), lp, AngularM64(0),  AngularJ64(nu),  AngularM64(0)) *
                                 AngularMomentum.ClebschGordan( AngularJ64(L), AngularM64(1), AngularJ64(Lp), AngularM64(-1), 
                                                                AngularJ64(nu), AngularM64(0)) *
                                 AngularMomentum.Wigner_6j(J, Jp, AngularJ64(nu), AngularJ64(Lp), AngularJ64(L), Jf) * 
                                 AngularMomentum.Wigner_6j(J, Jp, AngularJ64(nu), jp, j, Ji) * 
                                 AngularMomentum.Wigner_6j(j, jp, AngularJ64(nu), lp, l, AngularJ64(1//2)) * 
                                 conj(cha.amplitude) * chp.amplitude
            end
        end
        
        value = - 0.5 * wa / wn
        return( value )
    end


    """
    `PhotoRecombination.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                         settings::PhotoRecombination.Settings; output::Bool=true)` 
        ... to compute the photo recombination transition amplitudes and all properties as requested by the given settings. 
            A list of lines::Array{PhotoRecombination.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::PhotoRecombination.Settings; output::Bool=true)
        println("")
        printstyled("PhotoRecombination.computeLines(): The computation of photo-recombination properties starts now ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        # Check the consistency of the initial- and final-state multiplets
        PhotoRecombination.checkConsistentMultiplets(finalMultiplet, initialMultiplet)
        #
        lines = PhotoRecombination.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    PhotoRecombination.displayLines(lines)    end
        # Determine maximum energy and check for consistency of the grid
        maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
        nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newLines = PhotoRecombination.Line[]
        for  line in lines
            newLine = PhotoRecombination.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        PhotoRecombination.displayResults(stdout, newLines, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   PhotoRecombination.displayResults(iostream, newLines, settings)    end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end

    

    """
    `PhotoRecombination.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoRecombination.Settings)`  
        ... to determine a list of RecChannel for a transitions from the initial to final level and by taking into account 
            the particular settings of for this computation; an Array{PhotoRecombination.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoRecombination.Settings)
        channels = PhotoRecombination.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp in settings.multipoles
            for  gauge in settings.gauges
                symList = AngularMomentum.allowedMultipoleSymmetries(symf, mp)
                for  symt in symList
                    kappaList = AngularMomentum.allowedKappaSymmetries(symi, symt)
                    for  kappa in kappaList
                        # Include further restrictions if appropriate
                        if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      
                            push!(channels, PhotoRecombination.Channel(mp, Basics.Coulomb,   kappa, symt, 0., Complex(0.)) )
                        elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                            push!(channels, PhotoRecombination.Channel(mp, Basics.Babushkin, kappa, symt, 0., Complex(0.)) )  
                        elseif string(mp)[1] == 'M'                                
                            push!(channels, PhotoRecombination.Channel(mp, Basics.Magnetic,  kappa, symt, 0., Complex(0.)) ) 
                        end 
                    end
                end
            end
        end
        return( channels )  
    end


    """
    `PhotoRecombination.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoRecombination.Settings)`  
        ... to determine a list of PhotoRecombination.Line's for transitions between levels from the initial- and final-state 
            multiplets, and by taking into account the particular selections and settings for this computation; 
            an Array{PhotoRecombination.Line,1} is returned. Apart from the level specification, all physical properties are set 
            to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoRecombination.Settings)
        lines = PhotoRecombination.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    #
                    electronEnergies = Float64[]
                    if  settings.useIonEnergies 
                        for en in settings.ionEnergies
                            # E^(electron) [keV] = E^(projectile) [MeV/u] / 1.8228885
                            en_au = 1000. / 1.8228885 * Defaults.convertUnits("energy: from eV to atomic", en);      push!(electronEnergies, en_au)
                        end
                    else
                        for en in settings.electronEnergies 
                            en_au = Defaults.convertUnits("energy: to atomic", en);      push!(electronEnergies, en_au)
                        end
                        betaGamma2 = 1.0
                    end
                    #
                    for  en in electronEnergies
                        betaGamma2 = 1.0
                        if  true ## settings.useIonEnergies 
                            wc         = Defaults.getDefaults("speed of light: c")
                            Gamma      = 1.0 + en / wc^2
                            beta       = sqrt( 1.0 - 1.0/Gamma^2)
                            betaGamma2 = beta^2 * Gamma^2
                        end
                        #
                        if  en < 0    continue   end 
                        omega    = en + iLevel.energy - fLevel.energy
                        channels = PhotoRecombination.determineChannels(fLevel, iLevel, settings) 
                        push!( lines, PhotoRecombination.Line(iLevel, fLevel, en, omega, betaGamma2, EmProperty(0., 0.), true, channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `PhotoRecombination.displayLines(lines::Array{PhotoRecombination.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{PhotoRecombination.Line,1})
        nx = 181
        println(" ")
        println("  Selected photorecombination lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(10, "Energy_if"; na=2);              
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(12, "Energy e_r"; na=1);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(10, "omega"; na=5);              
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center( 7, "beta^2*"; na=2)
        sb = sb * TableStrings.center( 7, "gamma^2"; na=2)
        sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
        sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            energy = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))              * "   "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "   "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.2e", line.betaGamma2)                                                   * "  "
            kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
            for  i in 1:length(line.channels)
                push!( kappaMultipoleSymmetryList, (line.channels[i].kappa, line.channels[i].multipole, line.channels[i].gauge, 
                                                    line.channels[i].symmetry) )
            end
            wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
            sb = sa * wa[1];    println( sb )  
            for  i = 2:length(wa)
                sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
        end
        println("  ", TableStrings.hLine(nx), "\n")
        #
        return( nothing )
    end


    """
    `PhotoRecombination.displayResults(stream::IO, lines::Array{PhotoRecombination.Line,1}, settings::PhotoRecombination.Settings)`  
        ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is 
            returned otherwise.
    """
    function  displayResults(stream::IO, lines::Array{PhotoRecombination.Line,1}, settings::PhotoRecombination.Settings)
        if  settings.useIonEnergies   nx = 165    else  nx = 141   end
        println(stream, " ")
        println(stream, "  Photorecombination cross sections:")
        println(stream, "  ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(16, "i--J^P--f"   ; na=3);                       sb = sb * TableStrings.hBlank(19)
        sa = sa * TableStrings.center(12, "i--Energy--f"; na=4)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_r"; na=3)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
        if  settings.useIonEnergies
            sa = sa * TableStrings.center(12, "Ion energies"; na=3)             
            sb = sb * TableStrings.center(12, "[MeV/u]";      na=3)
        end
        sa = sa * TableStrings.center( 7, "beta^2*"; na=4)
        sb = sb * TableStrings.center( 7, "gamma^2"; na=4)
        sa = sa * TableStrings.flushleft(16, "Multipoles"; na=1);                      sb = sb * TableStrings.hBlank(17)
        sa = sa * TableStrings.center(30, "Cou -- Cross section -- Bab"; na=3)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " * 
                                              TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = " ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                          fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(17, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(17, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            if  settings.useIonEnergies 
                # E^(electron) [keV] = E^(projectile) [MeV/u] / 1.8228885
                enIon = 1.8228885 * Defaults.convertUnits("energy: from atomic to eV", line.electronEnergy) / 1000.
                sa = sa * @sprintf("%.6e", enIon) * "    "
            end
            sa = sa * @sprintf("%.2e", line.betaGamma2)                                                   * "   "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "                              "
            sa = sa * TableStrings.flushleft(16, mpString[1:16];  na=2)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Coulomb))     * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Babushkin))   * "    "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        #
        if  settings.calcAnisotropy  
            nx = 133
            println(stream, " ")
            println(stream, "  Anisotropy angular parameters beta_nu of the emitted photons for initially unpolarized ions:")
            println(stream, " ")
            println(stream, "  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
            sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(20)
            sa = sa * TableStrings.center(16, "i--J^P--f"   ; na=3);                       sb = sb * TableStrings.hBlank(19)
            sa = sa * TableStrings.center(12, "i--Energy--f"; na=4)               
            sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
            sa = sa * TableStrings.center(12, "omega"     ; na=4)             
            sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
            sa = sa * TableStrings.center(12, "Energy e_r"; na=3)             
            sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
            sa = sa * TableStrings.center(10, "Multipoles"; na=5);                         sb = sb * TableStrings.hBlank(15)
            sa = sa * TableStrings.flushleft(57, "beta's"; na=4)  
            println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
            #   
            for  line in lines
                sa  = " ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                              fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * TableStrings.center(17, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * TableStrings.center(17, TableStrings.symmetries_if(isym, fsym); na=3)
                en = line.initialLevel.energy - line.finalLevel.energy
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
                multipoles = EmMultipole[]
                for  ch in line.channels
                    multipoles = push!( multipoles, ch.multipole)
                end
                multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
                sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
                be = PhotoRecombination.computeAnisotropyParameter(1, Basics.Coulomb,   line);   sa = sa * @sprintf("%.4e", be.re) * " (1, C);  "
                be = PhotoRecombination.computeAnisotropyParameter(1, Basics.Babushkin, line);   sa = sa * @sprintf("%.4e", be.re) * " (1, B);  "
                println(stream, sa);   sa = TableStrings.hBlank(102)
                be = PhotoRecombination.computeAnisotropyParameter(2, Basics.Coulomb,   line);   sa = sa * @sprintf("%.4e", be.re) * " (2, C);  "
                be = PhotoRecombination.computeAnisotropyParameter(2, Basics.Babushkin, line);   sa = sa * @sprintf("%.4e", be.re) * " (2, B);  "
                println(stream, sa);   sa = TableStrings.hBlank(102)
                be = PhotoRecombination.computeAnisotropyParameter(3, Basics.Coulomb,   line);   sa = sa * @sprintf("%.4e", be.re) * " (3, C);  "
                be = PhotoRecombination.computeAnisotropyParameter(3, Basics.Babushkin, line);   sa = sa * @sprintf("%.4e", be.re) * " (3, B);  "
                println(stream, sa);   sa = TableStrings.hBlank(102)
                be = PhotoRecombination.computeAnisotropyParameter(4, Basics.Coulomb,   line);   sa = sa * @sprintf("%.4e", be.re) * " (4, C);  "
                be = PhotoRecombination.computeAnisotropyParameter(4, Basics.Babushkin, line);   sa = sa * @sprintf("%.4e", be.re) * " (4, B);  "
                println(stream, sa);   sa = TableStrings.hBlank(102)
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        #
        if  settings.calcTensors   
            println(stream, " ")
            println(stream, "  Reduced statistical tensors of the recombined ion ... not yet implemented !!")
            println(stream, " ")
        end
        #
        return( nothing )
    end

end # module
