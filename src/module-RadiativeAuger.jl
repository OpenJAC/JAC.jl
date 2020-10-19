
"""
`module  JAC.RadiativeAuger`  
    ... a submodel of JAC that contains all methods for computing radiative Auger and autoionization amplitudes and rates.
"""
module RadiativeAuger

    using Printf, ..AngularMomentum, ..AutoIonization, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..PhotoEmission, ..TableStrings

    """
    `struct  RadiativeAuger.Settings`  ... defines a type for the settings in estimating radiative-Auger and autoionization rates.

        + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
        + green                   ::Array{GreenChannel,1}   ... Precalculated and user-specified Green function of the ion.
        + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
        + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
        + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
        + operator                ::AbstractEeInteraction   ... Auger operator that is to be used for evaluating the Auger amplitudes: 
                                                                allowed values are: CoulombInteraction(), BreitInteraction(), ...
        + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    """
    struct Settings
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        green                     ::Array{GreenChannel,1}
        NoEnergySharings          ::Int64 
        maxKappa                  ::Int64  
        printBefore               ::Bool   
        operator                  ::AbstractEeInteraction 
        lineSelection             ::LineSelection 
    end 


    """
    `RadiativeAuger.Settings()`  ... constructor for the default values of RadiativeAuger line computations
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], GreenChannel[], 0, 3, false, LineSelection())
    end


    # `Base.show(io::IO, settings::RadiativeAuger.Settings)`  ... prepares a proper printout of the variable settings::RadiativeAuger.Settings.
    function Base.show(io::IO, settings::RadiativeAuger.Settings) 
        println(io, "multipoles:                   $(settings.multipoles)  ")
        println(io, "gauges:                       $(settings.gauges)  ")
        println(io, "green:                         (settings.green)  ")
        println(io, "NoEnergySharings:             $(settings.NoEnergySharings)  ")
        println(io, "maxKappa:                     $(settings.maxKappa)  ")
        println(io, "printBefore:                  $(settings.printBefore)  ")
        println(io, "operator :                    $(settings.operator )  ")
        println(io, "lineSelection:                $(settings.lineSelection)  ")
    end


    """
    `struct  RadiativeAuger.ReducedChannel`  
        ... defines a type for a RadiativeAuger (reduced) channel to help characterize a scattering (continuum) state of many 
            electron-states with a single free electron.

        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + multipole      ::EmMultipole          ... Multipole of the photon absorption.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + omega          ::Float64              ... photon energy
        + epsilon        ::Float64              ... (free) electron energy
        + kappa          ::Int64                ... partial-wave of the free electron
        + phase          ::Float64              ... phase of the partial wave
        + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
    """
    struct  Channel
        symmetry         ::LevelSymmetry
        multipole        ::EmMultipole  
        gauge            ::EmGauge
        omega            ::Float64 
        epsilon          ::Float64   
        kappa            ::Int64 
        phase            ::Float64   
        amplitude        ::Complex{Float64}
    end


    """
    `struct  Sharing`  
        ... defines a type for a RadiativeAuger sharing to help characterize energy sharing between the emitted photon and
            the scattering (continuum) state with a single free electron.

        + photonEnergy   ::Float64         ... Energy of the emitted photon.
        + electronEnergy ::Float64         ... Energy of the (outgoing free) electron.
        + differentialCs ::EmProperty      ... differential cross section of this energy sharing.
        + channels       ::Array{RadiativeAuger.Channel,1}  ... List of RadiativeAuger channels of this line.
    """
    struct  Sharing
        photonEnergy     ::Float64
        electronEnergy   ::Float64
        differentialCs   ::EmProperty
        channels         ::Array{RadiativeAuger.Channel,1}
    end




    # RadiativeAuger line between an initial and final (bound-state) level
    """
    `struct  Line`  ... defines a type for a radiative Auger line that may include the definition of sublines and their corresponding amplitudes.

        + initialLevel   ::Level          ... initial-(state) level
        + finalLevel     ::Level          ... final-(state) level
        + totalRate      ::EmProperty     ... Total rate of this line.
        + sharings       ::Array{RadiativeAuger.Sharing,1}  ... List of RadiativeAuger sharings of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        totalRate        ::EmProperty
        sharings         ::Array{RadiativeAuger.Sharing,1}
    end 


    # `Base.show(io::IO, line::RadiativeAuger.Line)`  ... prepares a proper printout of the variable line::RadiativeAuger.Line.
     function Base.show(io::IO, line::RadiativeAuger.Line) 
        println(io, "initialLevel:           $(line.initialLevel)  ")
        println(io, "finalLevel:             $(line.finalLevel)  ")
        println(io, "totalRate:              $(line.totalRate)  ")
        println(io, "sharings:               $(line.sharings)  ")
    end


    """
    `RadiativeAuger.computeAmplitudesProperties(line::RadiativeAuger.Line, grid::Radial.Grid, settings::RadiativeAuger.Settings)`  
        ... to compute all amplitudes and properties of the given line; a line::RadiativeAuger.Line is returned for which the amplitudes 
            and properties have now been evaluated.
    """
    function  computeAmplitudesProperties(line::RadiativeAuger.Line, grid::Radial.Grid, settings::RadiativeAuger.Settings)
        newSharings = RadiativeAuger.Sharing[]
        for sharing in line.sharings
            newChannels = RadiativeAuger.Channel[]
            for channel in sharing.channels
                # Generate a continuum orbital
                amplitude = 1.0 
                push!( newChannels, RadiativeAuger.Channel( channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, amplitude) )
            end
            push!( newSharings, RadiativeAuger.Sharing( sharing.photonEnergy, sharing.electronEnergy, EmProperty(-1., -1.), true, newChannels) )
        end
        # Calculate the totalRate 
        totalRate = EmProperty(-1., -1.)
        line = RadiativeAuger.Line( line.initialLevel, line.finalLevel, totalRate, true, newSharings)
        return( line )
    end


    """
    `RadiativeAuger.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                 settings::RadiativeAuger.Settings; output=true)`  
        ... to compute the radiative Auger transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{RadiativeAuger.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                           settings::RadiativeAuger.Settings; output=true)
        println("")
        printstyled("RadiativeAuger.computeLines(): The computation of radiative Auger rates starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = RadiativeAuger.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    RadiativeAuger.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = RadiativeAuger.Line[]
        for  line in lines
            newLine = RadiativeAuger.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        RadiativeAuger.displayResults(lines)
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `RadiativeAuger.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RadiativeAuger.Settings)` 
        ... to determine a list of RadiativeAuger.Line's for transitions between levels from the initial- and final-state multiplets, 
            and  by taking into account the particular selections and settings for this computation; an Array{RadiativeAuger.Line,1} is 
            returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RadiativeAuger.Settings)
        lines = RadiativeAuger.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    energy    = iLevel.energy - fLevel.energy
                    if  energy < settings.minAugerEnergy  ||  energy > settings.maxAugerEnergy    continue   end  
                    channels = RadiativeAuger.determineSharingsAndChannels(fLevel, iLevel, energy, settings) 
                    push!( lines, RadiativeAuger.Line(iLevel, fLevel, EmProperty(0., 0.,), true, channels) )
                end
            end
        end
        return( lines )
    end


    """
    `RadiativeAuger.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::RadiativeAuger.Settings)`  
        ... to determine a list of RadiativeAuger Sharing's and Channel's for a transitions from the initial to final level and by taking into 
            account the particular settings of for this computation; an Array{RadiativeAuger.Sharing,1} is returned.
    """
    function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::RadiativeAuger.Settings)
        sharings  = RadiativeAuger.Sharing[];    eSharings = Basics.determineEnergySharings(energy, settings.NoEnergySharings) 
        for  es in eSharings
            pEnergy   = es[1];    eEnergy = es[2] 
            channels  = RadiativeAuger.Channel[];   
            symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
            for  mp in settings.multipoles
                for  gauge in settings.gauges
                    symList = AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                    for  symt in symList
                        kappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
                        for  kappa in kappaList
                            # Include further restrictions if appropriate
                            if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      
                                push!(channels, RadiativeAuger.Channel(mp, Basics.Coulomb,   kappa, symt, 0., Complex(0.)) )
                            elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                                push!(channels, RadiativeAuger.Channel(mp, Basics.Babushkin, kappa, symt, 0., Complex(0.)) )  
                            elseif string(mp)[1] == 'M'                                
                                push!(channels, RadiativeAuger.Channel(mp, Basics.Magnetic,  kappa, symt, 0., Complex(0.)) ) 
                            end 
                        end
                    end
                end
            end
            push!(sharings, RadiativeAuger.Sharing(pEnergy, eEnergy, EmProperty(0., 0.), true, channels) )
        end
        return( sharings )  
    end


    """
    `RadiativeAuger.displayLines(lines::Array{RadiativeAuger.Line,1})`  
        ... to display a list of lines, sharings and channels that have been selected due to the prior settings. A neat table 
            of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{RadiativeAuger.Line,1})
        nx = 170
        println(" ")
        println("  Selected radiative-Auger lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(34, "Energy  " * TableStrings.inUnits("energy"); na=4);              
        sb = sb * TableStrings.center(34, "i -- f        omega     e_Auger  "; na=4)
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
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "  "
            #
            for  sharing  in  line.sharings
                sb =      @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.photonEnergy))   * "  "
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.electronEnergy)) * "    "
                kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
                for  i in 1:length(sharing.channels)
                    push!( kappaMultipoleSymmetryList, (sharing.channels[i].kappa, sharing.channels[i].multipole, sharing.channels[i].gauge, 
                                                        sharing.channels[i].symmetry) )
                end
                wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
                sc = sa * sb * wa[1];    println( sc )  
                for  i = 2:length(wa)
                    sc = TableStrings.hBlank( length(sa*sb) ) * wa[i];    println( sc )
                end
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `RadiativeAuger.displayResults(lines::Array{RadiativeAuger.Line,1})`  
        ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(lines::Array{RadiativeAuger.Line,1})
        nx = 148
        println(" ")
        println("  Radiative-Auger rates:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(36, "Energy  " * TableStrings.inUnits("energy"); na=4);              
        sb = sb * TableStrings.center(36, "i -- f        omega     e_Auger  "; na=4)
        sa = sa * TableStrings.center(30, "Cou -- differ. rate -- Bab"; na=3)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("rate") * "          " * TableStrings.inUnits("rate"); na=3)
        sa = sa * TableStrings.center(30, "Cou -- total rate -- Bab"; na=3)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("rate") * "          " * TableStrings.inUnits("rate"); na=3)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #  
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            energy = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
            #
            first = true
            for  sharing  in  line.sharings
                sb =      @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.photonEnergy))   * "  "
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.electronEnergy)) * "    "
                sb = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", sharing.differentialCs.Coulomb))     * "    "
                sb = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", sharing.differentialCs.Babushkin))   * "    "
                if  first                  first = false
                sb = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.totalRate.Coulomb))     * "    "
                sb = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.totalRate.Babushkin))   * "    "
                end 
                println(sa*sb)
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

end # module
