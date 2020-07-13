    #
    # Two-photon emission of initially-unpolarized atoms.
    #
    """
    `struct  MultiPhotonDeExcitation.Channel_2pEmission`  
        ... defines a type for a two-photon emission channel for the emission of photons with well-defined 
            multipolarities.

        + K              ::AngularJ64             ... Rank K of the channel.
        + omega1         ::Float64                ... omega1.
        + omega2         ::Float64                ... omega2.
        + multipole1     ::EmMultipole            ... Multipole M1.
        + multipole2     ::EmMultipole            ... Multipole M2.
        + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
        + Jsym           ::LevelSymmetry          ... Symmetry of the Green function channel used in the summation.
        + amplitude      ::Complex{Float64}       ... reduced two-photon emission amplitude U^(K, 2gamma, emission) (..)
                                                      associated with the given channel.
    """
    struct  Channel_2pEmission
        K                ::AngularJ64 
        omega1           ::Float64
        omega2           ::Float64 
        multipole1       ::EmMultipole
        multipole2       ::EmMultipole
        gauge            ::EmGauge
        Jsym             ::LevelSymmetry
        amplitude        ::Complex{Float64}
    end


    """
    `struct  MultiPhotonDeExcitation.Line_2pEmission`  
        ... defines a type for a two-photon absorption line by monochromatic light that may include the definition of channels.

        + initialLevel     ::Level                ... initial-(state) level
        + finalLevel       ::Level                ... final-(state) level
        + energy           ::Float64              ... Total transition energy.
        + omegas           ::Array{Float64,1}     ... Energy of the emitted photons.
        + weights          ::Array{Float64,1}     ... Associated weights for the given photon energy in a Gauss-Legendre integration.
        + energyDiffCs     ::Array{EmProperty,1}  ... energy-differential cross section for the selected omegas.
        + totalCs          ::EmProperty           ... Total cross section for right-circularly polarized incident light.
        + csUnpolarized    ::EmProperty           ... Total cross section.
        + hasChannels      ::Bool                 ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                      multipolarities, etc., or not.
        + channels         ::Array{MultiPhotonDeExcitation.Channel_2pEmission,1}  
                                                  ... List of MultiPhotonDeExcitation.Channel_2pEmission's of this line.
    """
    struct  Line_2pEmission
        initialLevel       ::Level
        finalLevel         ::Level
        energy             ::Float64
        omegas             ::Array{Float64,1} 
        weights            ::Array{Float64,1}
        energyDiffCs       ::Array{EmProperty,1} 
        totalCs            ::EmProperty 
        hasChannels        ::Bool
        channels           ::Array{MultiPhotonDeExcitation.Channel_2pEmission,1}
    end


    # `Base.show(io::IO, line::MultiPhotonDeExcitation.Line_2pEmission)`  
    #   ... prepares a proper printout of the variable line::MultiPhotonDeExcitation.Line_2pEmission.
    function Base.show(io::IO, line::MultiPhotonDeExcitation.Line_2pEmission) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "energy:            $(line.energy)  ")
        println(io, "omegas:            $(line.omegas)  ")
        println(io, "weights:           $(line.weights)  ")
        println(io, "energyDiffCs:      $(line.energyDiffCs)  ")
        println(io, "totalCs:           $(line.totalCs)  ")
        println(io, "hasChannels:       $(line.hasChannels)  ")
        println(io, "channels:          $(line.channels)  ")
    end


    """
    `MultiPhotonDeExcitation.computeAmplitudesProperties_2pEmission(line::MultiPhotonDeExcitation.Line_2pEmission, 
                             grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings)` 
        ... to compute all amplitudes and properties of the given line; a line::MultiPhotonDeExcitation.Line_2pEmission 
            is returned for which the amplitudes and properties are now evaluated.
    """
    function  computeAmplitudesProperties_2pEmission(line::MultiPhotonDeExcitation.Line_2pEmission, 
                                                     grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings)
        newChannels = MultiPhotonDeExcitation.Channel_2pEmission[]
        for channel in line.channels
            amplitude = MultiPhotonDeExcitation.computeReducedAmplitudeEmission(channel.K, line.finalLevel, channel.omega2, channel.multipole2, 
                                         channel.Jsym, channel.omega1, channel.multipole1, line.initialLevel, channel.gauge, grid, settings.greenChannels)
            push!( newChannels, MultiPhotonDeExcitation.Channel_2pEmission(channel.K, channel.omega1, channel.omega2, 
                                                        channel.multipole1, channel.multipole2, channel.gauge, channel.Jsym, amplitude) )
        end
        # Take amplitudes into accout for the given line
        newLine   = MultiPhotonDeExcitation.Line_2pEmission( line.initialLevel, line.finalLevel, line.energy, line.omegas, line.weights, 
                                                             line.energyDiffCs, line.totalCs, true, newChannels)
        # Calculate the requested cross sections, etc.
        energyDiffCs = EmProperty[]
        for  omega1  in  line.omegas
             omega2      = line.energy - omega1
            enDiffCs_Cou = MultiPhotonDeExcitation.computeEnergyDiffCs(omega1, omega2, newLine, EmGauge("Coulomb"), settings)
            enDiffCs_Bab = MultiPhotonDeExcitation.computeEnergyDiffCs(omega1, omega2, newLine, EmGauge("Babushkin"), settings)
            push!( energyDiffCs, EmProperty(enDiffCs_Cou, enDiffCs_Bab))
        end
        totalCs_Cou = 0.;   totalCs_Bab = 0.
        for  i = 1:length(line.omegas)
            totalCs_Cou = totalCs_Cou + energyDiffCs[i].Coulomb   * line.weights[i]
            totalCs_Bab = totalCs_Bab + energyDiffCs[i].Babushkin * line.weights[i]
        end
        totalCs   = EmProperty(totalCs_Cou, totalCs_Bab)
        newLine   = MultiPhotonDeExcitation.Line_2pEmission( line.initialLevel, line.finalLevel, line.energy, line.omegas, line.weights, energyDiffCs,
                                                             totalCs, true, newChannels)
        return( newLine )
    end



    """
    `MultiPhotonDeExcitation.computeLines(process::TwoPhotonEmission, finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                          grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings; output=true)` 
        ... to compute the multiphoton transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1} is returned.
    """
    function  computeLines(process::TwoPhotonEmission, finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                           grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings; output=true)
        println("")
        printstyled("MultiPhotonDeExcitation.computeLines(::TwoPhotonEmission): The computation of amplitudes starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = MultiPhotonDeExcitation.determineLines_2pEmission(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    MultiPhotonDeExcitation.displayLines_2pEmission(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = MultiPhotonDeExcitation.Line_2pEmission[]
        for  line in lines
            newLine = MultiPhotonDeExcitation.computeAmplitudesProperties_2pEmission(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        MultiPhotonDeExcitation.displayCrossSections_2pEmission(stdout, settings.process.properties, lines)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    MultiPhotonDeExcitation.displayCrossSections_2pEmission(iostream, settings.process.properties, lines)   end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `MultiPhotonDeExcitation.computeEnergyDiffCs(omega1::Float64, omega2::Float64, line::MultiPhotonDeExcitation.Line_2pEmission, 
                                                 settings::MultiPhotonDeExcitation.Settings)`  
        ... to compute the energy differential cross section for the given line and the energy sharing omega1. A dcs::Float64 is returned.
    """
    function computeEnergyDiffCs(omega1::Float64, omega2::Float64, line::MultiPhotonDeExcitation.Line_2pEmission, gauge::EmGauge,
                                 settings::MultiPhotonDeExcitation.Settings)
        dcs = 0.
        for  mp1 in settings.multipoles
            for  mp2 in settings.multipoles
                symi        = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
                symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                Klist       = oplus(line.finalLevel.J, line.initialLevel.J)
                for Jsym in symmetries
                    wk  = ComplexF64(0.);   wm = ComplexF64(0.)
                    for  K in Klist
                        wk =  wk + MultiPhotonDeExcitation.getReducedAmplitudeEmission(K, line.finalLevel, omega2, mp2, Jsym, omega1, mp1, 
                                                                                          line.initialLevel, gauge, line.channels) +
                               AngularMomentum.phaseFactor([K, +1, line.finalLevel.J, +1, line.initialLevel.J]) * 
                                   MultiPhotonDeExcitation.getReducedAmplitudeEmission(K, line.finalLevel, omega1, mp1, Jsym, omega2, mp2, 
                                                                                          line.initialLevel, gauge, line.channels)
                    end
                    dcs = dcs + abs( wk )^2
                end
            end
        end
        dcs = dcs * 2pi * Defaults.getDefaults("alpha")^2 / (AngularMomentum.twoJ(line.initialLevel.J) + 1) * omega1 * omega2
        
        return( dcs )
    end


    """
    `MultiPhotonDeExcitation.computeReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                                               omega1::Float64, multipole1::EmMultipole, initialLevel::Level,
                                                                               gauge::EmGauge, grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})`  
        ... to compute the reduced amplitude U^{K, 2gamma emission} (K, Jf, omega2, multipole2, Jsym, omega1, multipole1, Ji) by means of the
            given Green function channels. An amplitude::Complex{Float64} is returned.
    """
    function computeReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                               omega1::Float64, multipole1::EmMultipole, initialLevel::Level,
                                                               gauge::EmGauge, grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})
        U = Complex(0.);    found = false
        for channel in greenChannels
            if  Jsym == channel.symmetry
                found = true
                for (i, nuLevel) in enumerate(channel.gMultiplet.levels)
                    U = U + PhotoEmission.amplitude("emission", multipole2, gauge, omega2, finalLevel, nuLevel, grid) *
                            PhotoEmission.amplitude("emission", multipole1, gauge, omega1, nuLevel, initialLevel, grid) / 
                            (initialLevel.energy + omega1 - nuLevel.energy)
                end
            end
        end 
        
        if    found                                
              U = U * AngularMomentum.Wigner_6j(initialLevel.J, finalLevel.J, K, AngularJ64(multipole2.L), AngularJ64(multipole1.L), Jsym.J)
              ##x println(" ")
        else  println("No Green function cannel found for amplitude U^{K, 2gamma emission} (K, Jf, omega2, multipole2, Jsym = $Jsym, omega1, multipole1, Ji) ")
        end 
        
        return( U )
    end


    """
    `MultiPhotonDeExcitation.getReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                                           omega1::Float64, multipole1::EmMultipole, initialLevel::Level, 
                                                                           gauge::EmGauge, channels::Array{MultiPhotonDeExcitation.Channel_2pEmission,1})`  
        ... to get/return the reduced amplitude U^{K, 2gamma emission} (K, Jf, omega2, multipole2, Jsym, omega1, multipole1, Ji) from the calculated list
            of channels. An amplitude::Complex{Float64} is returned.
    """
    function getReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                           omega1::Float64, multipole1::EmMultipole, initialLevel::Level, 
                                                                           gauge::EmGauge, channels::Array{MultiPhotonDeExcitation.Channel_2pEmission,1})
        U = Complex(0.);    found = false
        for channel in channels
            if  K == channel.K  &&  omega1 == channel.omega1  &&  omega2 == channel.omega2  &&  Jsym == channel.Jsym  &&  
                multipole1 == channel.multipole1  &&  multipole2 == channel.multipole2      && 
                (gauge == channel.gauge  ||  EmGauge("Magnetic")  == channel.gauge)
                U = channel.amplitude;    found = true
            end
        end 
        
        if    found                                
              println("U^{K, 2gamma emission} (..) = $U found. ")
        else  println("No U^{K, 2gamma emission} (K, Jf, omega2, multipole2, Jsym = $Jsym, omega1, multipole1, Ji) amplitude found.")
        end 
        
        return( U )
    end


    """
    `MultiPhotonDeExcitation.determineChannels_2pEmission(finalLevel::Level, initialLevel::Level, energy::Float64, omegas::Array{Float64,1},
                                                          settings::MultiPhotonDeExcitation.Settings)`  
        ... to determine a list of MultiPhotonDeExcitation.Channel_2pEmission for a transitions from the initial to final level and by taking 
            into account the particular settings for this computation; an Array{MultiPhotonDeExcitation.Channel_2pEmission,1} is returned.
    """
    function determineChannels_2pEmission(finalLevel::Level, initialLevel::Level, energy::Float64, omegas::Array{Float64,1}, 
                                          settings::MultiPhotonDeExcitation.Settings)
        channels   = MultiPhotonDeExcitation.Channel_2pEmission[];   
        symi       = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  om1 in omegas
            om2 = energy - om1
            for  mp1 in settings.multipoles
                for  mp2 in settings.multipoles
                    symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                    Klist       = oplus(symf.J, symi.J)
                    for  symn in symmetries
                        hasMagnetic = false
                        for  gauge in settings.gauges
                            # Include further restrictions if appropriate
                            if     string(mp1)[1] == 'E' || string(mp2)[1] == 'E'  &&   gauge == Basics.UseCoulomb
                                for K in Klist  push!(channels, MultiPhotonDeExcitation.Channel_2pEmission(K, om1, om2, mp1, mp2, Basics.Coulomb, symn, 0.) )     end 
                            elseif string(mp1)[1] == 'E' || string(mp2)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                                for K in Klist  push!(channels, MultiPhotonDeExcitation.Channel_2pEmission(K, om1, om2, mp1, mp2, Basics.Babushkin, symn, 0.) )   end
                            elseif string(mp1)[1] == 'M' && string(mp2)[1] == 'M'
                                for K in Klist  push!(channels, MultiPhotonDeExcitation.Channel_2pEmission(K, om1, om2, mp1, mp2, Basics.Magnetic, symn, 0.) )    end
                            end
                        end
                    end 
                end
            end
        end

        return( channels )  
    end


    """
    `MultiPhotonDeExcitation.determineLines_2pEmission(finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                       settings::MultiPhotonDeExcitation.Settings)`
        ... to determine a list of MultiPhotonDeExcitation.Line_2pEmission's for transitions between the levels from the given 
            initial- and final-state multiplets and by taking into account the particular selections and settings for this computation; 
            an Array{MultiPhotonDeExcitation.Line_2pEmission,1} is returned. Apart from the level specification, all physical 
            properties are set to zero during this initialization process.  
    """
    function  determineLines_2pEmission(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::MultiPhotonDeExcitation.Settings)
        # The number of energy sharings are provided by the settings
        noSharings = settings.process.noSharings
        sharings   = QuadGK.gauss(noSharings)
    
        lines = MultiPhotonDeExcitation.Line_2pEmission[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    energy   = abs( iLevel.energy - fLevel.energy)
                    omegas   = Float64[];    weights = Float64[];    energyDiffCs = EmProperty[]
                    for  j = 1:noSharings
                        push!( omegas, sharings[1][j]*energy/2. + energy/2.);    push!( weights, sharings[1][j]*energy/2. )
                        push!( energyDiffCs, EmProperty(0., 0.))    
                    end
                    channels     = MultiPhotonDeExcitation.determineChannels_2pEmission(fLevel, iLevel, energy, omegas, settings) 
                    push!( lines, MultiPhotonDeExcitation.Line_2pEmission(iLevel, fLevel, energy, 
                                                               omegas, weights, energyDiffCs, EmProperty(0., 0.), true, channels) )
                end
            end
        end
        return( lines )
    end


    """
    `MultiPhotonDeExcitation.displayLines_2pEmission(lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines_2pEmission(lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1})
        nx = 170
        println(" ")
        println("  Selected two-photon emission lines with given photon splitting:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  ";   sc = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(12, "Energy"; na=2);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(12, "No channels"; na=2);                       sb = sb * TableStrings.hBlank(14)           
        sa = sa * TableStrings.flushleft(90, "Energy sharings/omega   " * TableStrings.inUnits("energy"); na=4)
        sc = sb[1:end-35] * TableStrings.hBlank(35)
        sb = sb * TableStrings.flushleft(90, "List of multipoles & intermediate level symmetries"; na=4)            
        sc = sc * TableStrings.flushleft(90, "(K-rank, multipole_1, Jsym, multipole_2, gauge), ..."; na=4)
        println(sa);    println(sb);    println(sc);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.energy)) * "   "
            sc = "         " * string(length(line.channels)) * "       "
            sa = sa * sc[end-12:end]
            sb = sa;  noEnergies = min(12, length(line.omegas))
            for  nn = 1:noEnergies    sb = sb * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", line.omegas[nn]))  * ",  "    end
            println( sb[1:end-3] )
            mpGaugeList = Tuple{AngularJ64, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge}[]
            for  channel in  line.channels
                push!( mpGaugeList, (channel.K, channel.multipole1, channel.Jsym, channel.multipole2, channel.gauge) )
            end
            wa = TableStrings.twoPhotonGaugeTupels(105, mpGaugeList)
            ##x if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
            for  i = 1:length(wa)
                sb = TableStrings.hBlank( length(sa) );    sb = sb * wa[i];    println( sb )
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `MultiPhotonDeExcitation.displayCrossSections_2pEmission(stream::IO, 
                                                  properties::Array{MultiPhotonDeExcitation.AbstractMultiPhotonProperty,1},
                                                  lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1})`  
        ... to display all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is 
            returned otherwise.
    """
    function  displayCrossSections_2pEmission(stream::IO, properties::Array{AbstractMultiPhotonProperty,1},
                                              lines::Array{Line_2pEmission,1})
        nx = 95
        println(stream, " ")
        println(stream, "  Two-photon emission with given splitting of photon energies:")
        println(stream, " ")
        println(stream, "  Energy-differential cross sections:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(10, "Energy"; na=4);              
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
        ##x sa = sa * TableStrings.center(10, "omega";  na=4);              
        ##x sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(28, "Cou -- cross section -- Bab"; na=4);              
        sb = sb * TableStrings.center(28, TableStrings.inUnits("cross section") * "           " * TableStrings.inUnits("cross section"); na=4)

        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=4)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.energy))        * "    "
            ##x sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega))         * "      "
            ##x sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic",   line.csLinear.Coulomb))          * "    "
            ##x sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic",   line.csLinear.Babushkin))        * "        "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end
