    #
    # Two-photon absorption by monochromatic and equally-polarized photons, usually from the same beam.
    #
    """
    `struct  MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic`  
        ... defines a type for a two-photon absorption channel for the absorption of monochromatic light with well-defined 
            multipolarities.

        + K              ::AngularJ64             ... Rank K of the channel.
        + omega          ::Float64                ... omega.
        + multipole1     ::EmMultipole            ... Multipole M1.
        + multipole2     ::EmMultipole            ... Multipole M2.
        + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
        + Jsym           ::LevelSymmetry          ... Symmetry of the Green function channel used in the summation.
        + amplitude      ::Complex{Float64}       ... reduced two-photon absorption amplitude U^(K, 2gamma, absorption) (..)
                                                      associated with the given channel.
    """
    struct  Channel_2pAbsorptionMonochromatic
        K                ::AngularJ64 
        omega            ::Float64
        multipole1       ::EmMultipole
        multipole2       ::EmMultipole
        gauge            ::EmGauge
        Jsym             ::LevelSymmetry
        amplitude        ::Complex{Float64}
    end


    """
    `struct  MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic`  
        ... defines a type for a two-photon absorption line by monochromatic light that may include the definition of channels.

        + initialLevel     ::Level          ... initial-(state) level
        + finalLevel       ::Level          ... final-(state) level
        + energy           ::Float64        ... Total transition energy.
        + omega            ::Float64        ... Energy of the incoming photons.
        + csLinear         ::EmProperty     ... Total cross section for linearly-polarized incident light.
        + csRightCircular  ::EmProperty     ... Total cross section for right-circularly polarized incident light.
        + csUnpolarized    ::EmProperty     ... Total cross section for unpolarized incident light.
        + hasChannels      ::Bool           ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                multipolarities, etc., or not.
        + channels         ::Array{MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic,1}  
                                            ... List of MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic's of this line.
    """
    struct  Line_2pAbsorptionMonochromatic
        initialLevel       ::Level
        finalLevel         ::Level
        energy             ::Float64
        omega              ::Float64
        csLinear           ::EmProperty
        csRightCircular    ::EmProperty
        csUnpolarized      ::EmProperty
        hasChannels        ::Bool
        channels           ::Array{MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic,1}
    end


    # `Base.show(io::IO, line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic)`  
    #   ... prepares a proper printout of the variable line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic.
    function Base.show(io::IO, line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "energy:            $(line.energy)  ")
        println(io, "omega:             $(line.omega)  ")
        println(io, "csLinear:          $(line.csLinear)  ")
        println(io, "csRightCircular:   $(line.csRightCircular)  ")
        println(io, "csUnpolarized:     $(line.csUnpolarized)  ")
        println(io, "hasChannels:       $(line.hasChannels)  ")
        println(io, "channels:          $(line.channels)  ")
    end


    """
    `MultiPhotonDeExcitation.computeAmplitudesProperties_2pAbsorptionMonochromatic(
                             line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic, grid::Radial.Grid, 
                             settings::MultiPhotonDeExcitation.Settings)` 
        ... to compute all amplitudes and properties of the given line; a line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic 
            is returned for which the amplitudes and properties are now evaluated.
    """
    function  computeAmplitudesProperties_2pAbsorptionMonochromatic(line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic, 
                                                                    grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings)
        newChannels = MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic[]
        for channel in line.channels
            amplitude = MultiPhotonDeExcitation.computeReducedAmplitudeAbsorption(channel.K, line.finalLevel, channel.multipole2, 
                                         channel.Jsym, channel.omega, channel.multipole1, line.initialLevel, channel.gauge, grid, settings.greenChannels)
            push!( newChannels, MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic(channel.K, channel.omega, channel.multipole1, channel.multipole2, 
                                                                                          channel.gauge, channel.Jsym, amplitude) )
        end
        # Calculate the requested cross sections, etc.
        csLinear        = EmProperty(0., 0.)
        csRightCircular = EmProperty(0., 0.)
        csUnpolarized   = EmProperty(0., 0.)
        for property in settings.process.properties
            if      typeof(property) == TotalCsLinear          
                totalCs_Cou = MultiPhotonDeExcitation.computeTotalCsLinear(line.omega, line, EmGauge("Coulomb"), settings)
                totalCs_Bab = MultiPhotonDeExcitation.computeTotalCsLinear(line.omega, line, EmGauge("Babushkin"), settings)
                csLinear    = EmProperty( totalCs_Cou,  totalCs_Bab)
            elseif  typeof(property) == TotalCsRightCircular   
                totalCs_Cou = MultiPhotonDeExcitation.computeTotalCsRightCircular(line.omega, line, EmGauge("Coulomb"), settings)
                totalCs_Bab = MultiPhotonDeExcitation.computeTotalCsRightCircular(line.omega, line, EmGauge("Babushkin"), settings)
                csRightCircular = EmProperty( totalCs_Cou,  totalCs_Bab)
            elseif  typeof(property) == TotalCsUnpolarized     
                totalCs_Cou = MultiPhotonDeExcitation.computeTotalCsUnpolarized(line.omega, line, EmGauge("Coulomb"), settings)
                totalCs_Bab = MultiPhotonDeExcitation.computeTotalCsUnpolarized(line.omega, line, EmGauge("Babushkin"), settings)
                csUnpolarized = EmProperty( totalCs_Cou,  totalCs_Bab)
            end
        end
        line = MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic( line.initialLevel, line.finalLevel, line.energy, line.omega, 
                                                                       csLinear, csRightCircular, csUnpolarized, true, newChannels)
        return( line )
    end



    """
    `MultiPhotonDeExcitation.computeLines(process::TwoPhotonAbsorptionMonochromatic, finalMultiplet::Multiplet, 
                                          initialMultiplet::Multiplet, grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings; 
                                          output=true)` 
        ... to compute the multiphoton transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic,1} is returned.
    """
    function  computeLines(process::TwoPhotonAbsorptionMonochromatic, finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                           grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings; output=true)
        println("")
        printstyled("MultiPhotonDeExcitation.computeLines(::TwoPhotonAbsorptionMonochromatic): The computation of amplitudes starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = MultiPhotonDeExcitation.determineLines_2pAbsorptionMonochromatic(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    MultiPhotonDeExcitation.displayLines_2pAbsorptionMonochromatic(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic[]
        for  line in lines
            newLine = MultiPhotonDeExcitation.computeAmplitudesProperties_2pAbsorptionMonochromatic(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        MultiPhotonDeExcitation.displayCrossSections_2pAbsorptionMonochromatic(stdout, settings.process.properties, lines)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    MultiPhotonDeExcitation.displayCrossSections_2pAbsorptionMonochromatic(iostream, settings.process.properties, lines)   end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `MultiPhotonDeExcitation.computeTotalCsLinear(omega::Float64, line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic, 
                                                  gauge::EmGauge, settings::MultiPhotonDeExcitation.Settings)`  
        ... to compute the total cross sections for linearly-polarized incident light. A tcs::Float64 is returned.
    """
    function computeTotalCsLinear(omega::Float64, line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic, 
                                  gauge::EmGauge, settings::MultiPhotonDeExcitation.Settings)
        tcs = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J)
        symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
        
        for  K in Klist
            qList = AngularMomentum.m_values(K)
            for  q in qList
                amp = ComplexF64(0.)
                for  lambda1  in [-1, 1]
                    for  lambda2  in [-1, 1]
                        for  mp1 in settings.multipoles
                            for  mp2 in settings.multipoles
                                if   mp1.electric   p1 = 1    else    p1 = 0    end
                                if   mp2.electric   p2 = 1    else    p2 = 0    end
                                symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                                for Jsym in symmetries
                                    wa = (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 
                                    wb = sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (2*AngularMomentum.twoJ(K) + 1)
                                    wc = AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)
                                    wd = MultiPhotonDeExcitation.getReducedAmplitudeAbsorption(K, line.finalLevel, mp2, Jsym, omega, mp1, 
                                                                                               line.initialLevel, gauge, line.channels) 
                                    println("computeTotalCsLinear: L1, L2, K, lambda1, lambda2, q = $(mp1.L), $(mp2.L), $K, $lambda1, $lambda2, $q")
                                    println("computeTotalCsLinear: wa, wb, wc, wd = $wa  $wb  $wc  $wd")
                                                                                                         
                                    amp = amp + (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2           *
                                                sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (2*AngularMomentum.twoJ(K) + 1)       *
                                                AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)             *
                                                MultiPhotonDeExcitation.getReducedAmplitudeAbsorption(K, line.finalLevel, mp2, Jsym, omega, mp1, 
                                                                                                         line.initialLevel, gauge, line.channels) 
                                end
                            end
                        end
                    end
                end
                tcs = tcs + abs( amp )^2
            end
        end
        
        println("computeTotalCsLinear: tcs = $tcs")
        tcs = tcs * 8*pi^5 * Defaults.getDefaults("alpha")^2 / (AngularMomentum.twoJ(line.initialLevel.J) + 1) / omega^2
        
        return( tcs )
    end


    """
    `MultiPhotonDeExcitation.computeTotalCsRightCircular(omega::Float64, line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic, 
                                                         gauge::EmGauge, settings::MultiPhotonDeExcitation.Settings)`  
        ... to compute the total cross sections for right-cicularly polarized incident light. A tcs::Float64 is returned.
    """
    function computeTotalCsRightCircular(omega::Float64, line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic, 
                                         gauge::EmGauge, settings::MultiPhotonDeExcitation.Settings)
        tcs = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J)
        symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
        
        # Need to be filled
        
        tcs = tcs * 8*pi^5 * Defaults.getDefaults("alpha")^2 / (AngularMomentum.twoJ(line.initialLevel.J) + 1) / omega^2
        
        return( tcs )
    end


    """
    `MultiPhotonDeExcitation.computeTotalCsUnpolarized(omega::Float64, line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic, 
                                                       gauge::EmGauge, settings::MultiPhotonDeExcitation.Settings)`  
        ... to compute the total cross sections for linearly-polarized incident light. A tcs::Float64 is returned.
    """
    function computeTotalCsUnpolarized(omega::Float64, line::MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic, 
                                       gauge::EmGauge, settings::MultiPhotonDeExcitation.Settings)
        tcs = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J)
        symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
        
        for  K in Klist
            qList = AngularMomentum.m_values(K)
            for  q in qList
                for  lambda1  in [-1, 1]
                    for  lambda2  in [-1, 1]
                        amp = ComplexF64(0.)
                        for  mp1 in settings.multipoles
                            for  mp2 in settings.multipoles
                                if   mp1.electric   p1 = 1    else    p1 = 0    end
                                if   mp2.electric   p2 = 1    else    p2 = 0    end
                                symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                                for Jsym in symmetries
                                    amp = amp + (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 *
                                                sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (2*AngularMomentum.twoJ(K) + 1)       *
                                                AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)             *
                                                MultiPhotonDeExcitation.getReducedAmplitudeAbsorption(K, line.finalLevel, mp2, Jsym, omega, mp1, 
                                                                                                         line.initialLevel, gauge, line.channels) 
                                end
                            end
                        end
                        tcs = tcs + abs( amp )^2
                    end
                end
            end
        end
        
        tcs = tcs * 8*pi^5 * Defaults.getDefaults("alpha")^2 / (AngularMomentum.twoJ(line.initialLevel.J) + 1) / omega^2
        
        return( tcs )
    end


    """
    `MultiPhotonDeExcitation.computeReducedAmplitudeAbsorption(K::AngularJ64, finalLevel::Level, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                                 omega::Float64, multipole1::EmMultipole, initialLevel::Level,
                                                                 gauge::EmGauge, grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})`  
        ... to compute the reduced amplitude U^{K, 2gamma emission} (K, Jf, multipole2, Jsym, omega, multipole1, Ji) by means of the
            given Green function channels. An amplitude::Complex{Float64} is returned.
    """
    function computeReducedAmplitudeAbsorption(K::AngularJ64, finalLevel::Level, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                omega::Float64, multipole1::EmMultipole, initialLevel::Level,
                                                                gauge::EmGauge, grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})
        U = Complex(0.);    found = false
        for channel in greenChannels
            if  Jsym == channel.symmetry
                found = true
                for (i, nuLevel) in enumerate(channel.gMultiplet.levels)
                    U = U + PhotoEmission.amplitude("absorption", multipole2, gauge, omega, finalLevel, nuLevel, grid) *
                            PhotoEmission.amplitude("absorption", multipole1, gauge, omega, nuLevel, initialLevel, grid) / 
                            (initialLevel.energy + omega - nuLevel.energy)
                end
            end
        end 
        
        if    found                                
              U = U * AngularMomentum.Wigner_6j(initialLevel.J, finalLevel.J, K, AngularJ64(multipole2.L), AngularJ64(multipole1.L), Jsym.J)
              ##x println(" ")
        else  println("No Green function cannel found for amplitude U^{K, 2gamma absorption} (K, Jf, multipole2, Jsym = $Jsym, omega, multipole1, Ji) ")
        end 
        
        return( U )
    end


    """
    `MultiPhotonDeExcitation.getReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, multipole2::EmMultipole, Jsym::LevelSymmetry, omega::Float64, 
                                                                                           multipole1::EmMultipole, initialLevel::Level, 
                                                         gauge::EmGauge, channels::Array{MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic,1})`  
        ... to get/return the reduced amplitude U^{K, 2gamma emission} (K, Jf, multipole2, Jsym, omega, multipole1, Ji) from the calculated list
            of channels. An amplitude::Complex{Float64} is returned.
    """
    function getReducedAmplitudeAbsorption(K::AngularJ64, finalLevel::Level, multipole2::EmMultipole, Jsym::LevelSymmetry, omega::Float64, 
                                                                             multipole1::EmMultipole, initialLevel::Level, 
                                           gauge::EmGauge, channels::Array{MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic,1})
        U = Complex(0.);    found = false
        for channel in channels
            if  K == channel.K  &&  omega == channel.omega  &&  Jsym == channel.Jsym       &&  
                multipole1 == channel.multipole1  &&  multipole2 == channel.multipole2     &&  
                (gauge == channel.gauge  ||  EmGauge("Magnetic")  == channel.gauge)
                U = channel.amplitude;    found = true
            end
        end 
        
        if    found                                
              println("U^{K, 2gamma absorption} (..) found. ")
        else  println("No U^{K, 2gamma absorption} (K, Jf, multipole2, Jsym = $Jsym, omega, multipole1, Ji) amplitude found.")
        end 
        
        return( U )
    end


    """
    `MultiPhotonDeExcitation.determineChannels_2pAbsorptionMonochromatic(omega::Float64, finalLevel::Level, initialLevel::Level, 
                                                                         settings::MultiPhotonDeExcitation.Settings)`  
        ... to determine a list of MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic for a transitions from the initial to 
            final level and by taking into account the particular settings of for this computation; 
            an Array{MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic,1} is returned.
    """
    function determineChannels_2pAbsorptionMonochromatic(omega::Float64, finalLevel::Level, initialLevel::Level, settings::MultiPhotonDeExcitation.Settings)
        channels   = MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic[];   
        symi       = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp1 in settings.multipoles
            for  mp2 in settings.multipoles
                symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                Klist       = oplus(symf.J, symi.J)
                for  symn in symmetries
                    hasMagnetic = false
                    for  gauge in settings.gauges
                        # Include further restrictions if appropriate
                        if     string(mp1)[1] == 'E' || string(mp2)[1] == 'E'  &&   gauge == Basics.UseCoulomb
                            for K in Klist  push!(channels, MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic(K, omega, mp1, mp2, Basics.Coulomb, symn, 0.) )     end 
                        elseif string(mp1)[1] == 'E' || string(mp2)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                            for K in Klist  push!(channels, MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic(K, omega, mp1, mp2, Basics.Babushkin, symn, 0.) )   end
                        elseif string(mp1)[1] == 'M' && string(mp2)[1] == 'M'
                            for K in Klist  push!(channels, MultiPhotonDeExcitation.Channel_2pAbsorptionMonochromatic(K, omega, mp1, mp2, Basics.Magnetic, symn, 0.) )    end
                        end
                    end 
                end
            end
        end

        return( channels )  
    end


    """
    `MultiPhotonDeExcitation.determineLines_2pAbsorptionMonochromatic(finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                                      settings::MultiPhotonDeExcitation.Settings)`
        ... to determine a list of MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic's for transitions between the levels from the given 
            initial- and final-state multiplets and by taking into account the particular selections and settings for this computation; 
            an Array{MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic,1} is returned. Apart from the level specification, all physical 
            properties are set to zero during this initialization process.  
    """
    function  determineLines_2pAbsorptionMonochromatic(finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                       settings::MultiPhotonDeExcitation.Settings)
        if    settings.selectLines    selectLines   = true;   selectedLines = Basics.determine("selected lines", settings.selectedLines)
        else                          selectLines   = false
        end
    
        lines = MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic[]
        for  i = 1:length(initialMultiplet.levels)
            for  f = 1:length(finalMultiplet.levels)
                if  selectLines  &&  !(haskey(selectedLines, (i,f)) )    continue   end
                energy   = abs( initialMultiplet.levels[i].energy - finalMultiplet.levels[f].energy)
                omega    = abs( initialMultiplet.levels[i].energy - finalMultiplet.levels[f].energy) / 2.
                channels = MultiPhotonDeExcitation.determineChannels_2pAbsorptionMonochromatic(omega, finalMultiplet.levels[f], initialMultiplet.levels[i], settings) 
                push!( lines, MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic(initialMultiplet.levels[i], finalMultiplet.levels[f], energy, omega,
                                                                                     EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), true, channels) )
            end
        end
        return( lines )
    end


    """
    `MultiPhotonDeExcitation.displayLines_2pAbsorptionMonochromatic(lines::Array{MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines_2pAbsorptionMonochromatic(lines::Array{MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic,1})
        println(" ")
        println("  Selected two-photon absorption lines by monochromatic and equally-polarized photons:")
        println(" ")
        println("  ", TableStrings.hLine(175))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(12, "Energy"; na=4);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"; na=4);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.flushleft(90, "List of multipoles & intermediate level symmetries"; na=4)            
        sb = sb * TableStrings.flushleft(90, "(K-rank, multipole_1, Jsym, multipole_2, gauge), ..."; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(175)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.energy)) * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.omega))  * "    "
            mpGaugeList = Tuple{AngularJ64, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge}[]
            for  channel in  line.channels
                push!( mpGaugeList, (channel.K, channel.multipole1, channel.Jsym, channel.multipole2, channel.gauge) )
            end
            ##x println("displayLines_2pAbsorptionMonochromatic: mpGaugeList = ", mpGaugeList)
            wa = TableStrings.twoPhotonGaugeTupels(105, mpGaugeList)
            if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
            for  i = 2:length(wa)
                sb = TableStrings.hBlank( length(sa) );    sb = sb * wa[i];    println( sb )
            end
        end
        println("  ", TableStrings.hLine(175))
        #
        return( nothing )
    end


    """
    `MultiPhotonDeExcitation.displayCrossSections_2pAbsorptionMonochromatic(stream::IO, 
                                                  properties::Array{MultiPhotonDeExcitation.AbstractMultiPhotonProperty,1},
                                                  lines::Array{MultiPhotonDeExcitation.Line_2pAbsorptionMonochromatic,1})`  
        ... to display all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is 
            returned otherwise.
    """
    function  displayCrossSections_2pAbsorptionMonochromatic(stream::IO, properties::Array{AbstractMultiPhotonProperty,1},
                                                             lines::Array{Line_2pAbsorptionMonochromatic,1})
        println(stream, " ")
        println(stream, "  Two-photon absorption by monochromatic and equally-polarized photons (usually from the same beam):")
        println(stream, " ")
        println(stream, "  Cross sections are given for:")
        noCs = 0  # Number of cross sections to be printed
        for property in properties
            if      typeof(property) == TotalCsLinear          
                noCs = noCs + 1;   println(stream, "    + total cross sections for linearly-polarized incident light")  
            elseif  typeof(property) == TotalCsRightCircular   
                noCs = noCs + 1;   println(stream, "    + total cross sections for right-circularly polarized incident light") 
            elseif  typeof(property) == TotalCsUnpolarized     
                noCs = noCs + 1;   println(stream, "    + total cross sections for unpolarized incident light") 
            end
        end
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(75 + 30noCs))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(10, "Energy"; na=4);              
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(10, "omega";  na=4);              
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
        for no = 1:noCs
            sa = sa * TableStrings.center(28, "Cou -- cross section -- Bab"; na=4);              
            sb = sb * TableStrings.center(28, TableStrings.inUnits("cross section") * "           " * TableStrings.inUnits("cross section"); na=4)
        end

        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(75 + 30noCs)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=4)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.energy))        * "    "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega))         * "      "
            for property in properties
                if      typeof(property) == TotalCsLinear          
                    sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic",   line.csLinear.Coulomb))          * "    "
                    sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic",   line.csLinear.Babushkin))        * "        "
                elseif  typeof(property) == TotalCsRightCircular
                    sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic",   line.csRightCircular.Coulomb))   * "    "
                    sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic",   line.csRightCircular.Babushkin)) * "        "
                elseif  typeof(property) == TotalCsUnpolarized 
                    sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic",   line.csUnpolarized.Coulomb))     * "    "
                    sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic",   line.csUnpolarized.Babushkin))   * "        "
                end
            end
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(75 + 30noCs))
        #
        return( nothing )
    end
