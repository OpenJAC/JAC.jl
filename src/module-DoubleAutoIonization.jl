
"""
`module  JAC.DoubleAutoIonization`  
    ... a submodel of JAC that contains all methods for computing double Auger and autoionization amplitudes and rates.
"""
module DoubleAutoIonization

    using Printf, JAC, ..AngularMomentum, ..AtomicState, ..AutoIonization, ..Basics, ..ManyElectron, ..TableStrings


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

        + energy1        ::Float64              ... energy of the partial wave_1
        + kappa1         ::Int64                ... partial-wave of the free electron_1
        + phase1         ::Float64              ... phase of the partial wave_1
        + xSymmetry      ::LevelSymmetry        ... angular momentum and parity after coupling electron1.
        + energy2        ::Float64              ... energy of the partial wave_2
        + kappa2         ::Int64                ... partial-wave of the free electron_2
        + phase2         ::Float64              ... phase of the partial wave_2
        + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
    """
    struct  ReducedChannel
        energy1          ::Float64 
        kappa1           ::Int64 
        phase1           ::Float64
        xSymmetry        ::LevelSymmetry 
        energy2          ::Float64 
        kappa2           ::Int64   
        phase2           ::Float64 
        amplitude        ::Complex{Float64}
    end


    # `Base.show(io::IO, channel::DoubleAutoIonization.ReducedChannel)`  
    #       ... prepares a proper printout of the variable cannel::DoubleAutoIonization.ReducedChannel.
    function Base.show(io::IO, channel::DoubleAutoIonization.ReducedChannel) 
        ## sa = "reduced DA channel for J^P=(channel.tSymmetry)  with " * 
        ##      "kappa1=$(channel.kappa1), energy1=$(channel.energy1), phase1=$(channel.phase1), " *
        ##      "kappa2=$(channel.kappa2), energy2=$(channel.energy2), phase2=$(channel.phase2), amp=$(channel.amplitude) "
        ## println(io, sa)
        println(io, "multipole:                    $(channel.multipole)  ")
        println(io, "gauge:                        $(channel.gauge)  ")
        println(io, "omega:                        $(channel.omega)  ")
        println(io, "energy1:                      $(channel.energy1)  ")
        println(io, "kappa1:                       $(channel.kappa1)  ")
        println(io, "phase1:                       $(channel.phase1)  ")
        println(io, "xSymmetry:                    $(channel.xSymmetry)  ")
        println(io, "energy2:                      $(channel.energy2)  ")
        println(io, "kappa2:                       $(channel.kappa2)  ")
        println(io, "phase2:                       $(channel.phase2)  ")
        println(io, "amplitude:                    $(channel.amplitude)  ")
    end


    """
    `struct  DoubleAutoIonization.Sharing`  
        ... defines a type for a DoubleAutoIonization sharing to help characterize energy sharing between the two emitted electrons.

        + epsilon1         ::Float64         ... Energy of (free) electron 1.
        + epsilon2         ::Float64         ... Energy of (free) electron 2.
        + weight           ::Float64         ... Gauss-Lengendre weight of this sharing for energy-integrated quantities.
        + differentialRate ::EmProperty      ... differential cross section of this energy sharing.
        + channels         ::Array{DoubleAutoIonization.ReducedChannel,1}  ... List of DoubleAutoIonization (reduced) channels of this line.
    """
    struct  Sharing
        epsilon1           ::Float64
        epsilon2           ::Float64
        weight             ::Float64 
        differentialRate   ::EmProperty
        channels           ::Array{DoubleAutoIonization.ReducedChannel,1}
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
    `DoubleAutoIonization.computeAmplitudesProperties(line::DoubleAutoIonization.Line, grid::Radial.Grid, settings::DoubleAutoIonization.Settings)`  
        ... to compute all amplitudes and properties of the given line; a line::DoubleAutoIonization.Line is returned for which the amplitudes 
            and properties have now been evaluated.
    """
    function  computeAmplitudesProperties(line::DoubleAutoIonization.Line, grid::Radial.Grid, settings::DoubleAutoIonization.Settings)
        newSharings = DoubleAutoIonization.Sharing[]
        for sharing in line.sharings
            newChannels = DoubleAutoIonization.ReducedChannel[]
            for ch in sharing.channels
                # Generate a continuum orbital
                phase     = 1.3
                amplitude = 1.0im 
                push!( newChannels, ReducedChannel(ch.energy1, ch.kappa1, phase, ch.xSymmetry, ch.energy2, ch.kappa2, phase, amplitude) )
            end
            # Calculate the differential rate 
            diffRate = EmProperty(-1., -1.)
            push!( newSharings, DoubleAutoIonization.Sharing(sharing.epsilon1, sharing.epsilon2, sharing.weight, diffRate, newChannels) )
        end
        # Calculate the total rate
        totalRate = EmProperty(-1., -1.)
        line = DoubleAutoIonization.Line(line.initialLevel, line.finalLevel, totalRate, newSharings)
        return( line )
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
        printstyled("DoubleAutoIonization.computeLines(): The computation of double-Auger rates starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------- \n", color=:light_green)
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
        DoubleAutoIonization.displayTotalRates(stdout, lines, settings)
        DoubleAutoIonization.displayDifferentialRates(stdout, lines, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   DoubleAutoIonization.displayTotalRates(iostream, lines, settings)
                           DoubleAutoIonization.displayDifferentialRates(iostream, lines, settings)     end
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
                    sharings = DoubleAutoIonization.determineSharingsAndChannels(fLevel, iLevel, energy, settings) 
                    push!( lines, DoubleAutoIonization.Line(iLevel, fLevel, EmProperty(0., 0.,), sharings) )
                end
            end
        end
        return( lines )
    end


    """
    `DoubleAutoIonization.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::DoubleAutoIonization.Settings)`  
        ... to determine a list of DoubleAutoIonization Sharing's and Channel's for a transitions from the initial to final level and by taking into 
            account the particular settings of for this computation; an Array{DoubleAutoIonization.Sharing,1} is returned.
    """
    function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::DoubleAutoIonization.Settings)
        sharings = DoubleAutoIonization.Sharing[];    eSharings = Basics.determineEnergySharings(energy, settings.NoEnergySharings) 
        for  es in eSharings
            epsilon1  = es[1];    epsilon2 = es[2];    weight = es[3] 
            channels  = DoubleAutoIonization.ReducedChannel[];   
            symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
            couplings = AngularMomentum.allowedDoubleKappaCouplingSequence(symf, symi, settings.maxKappa)
            for  coupling  in  couplings
                push!(channels, ReducedChannel(epsilon1, coupling[1], 0., coupling[2], epsilon2, coupling[3], 0., Complex(0.)) ) 
            end
            push!(sharings, DoubleAutoIonization.Sharing(epsilon1, epsilon2, weight, EmProperty(0., 0.), channels) )
        end
        return( sharings )  
    end


    """
    `DoubleAutoIonization.displayDifferentialRates(stream::IO, lines::Array{DoubleAutoIonization.Line,1}, settings::DoubleAutoIonization.Settings)`  
        ... to display all differential rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayDifferentialRates(stream::IO, lines::Array{DoubleAutoIonization.Line,1}, settings::DoubleAutoIonization.Settings)
        #
        # First, print lines and sharings
        nx = 128
        println(stream, " ")
        println(stream, "  Energy-differential rates of selected double-Auger lines:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.flushleft(38, "Energies (all in " * TableStrings.inUnits("energy") * ")"; na=3);              
        sb = sb * TableStrings.flushleft(38, "  i -- f        epsilon_1     epsilon_2"; na=3)
        sa = sa * TableStrings.center(14, "Weight"; na=0);                            sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(34, "Cou --  diff. rate  -- Bab"; na=3)      
        sb = sb * TableStrings.center(34, TableStrings.inUnits("rate") * "        " * 
                                          TableStrings.inUnits("rate"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
            energy = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
            #
            for  (is, sharing)  in  enumerate(line.sharings)
                if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1))                 * "    "
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2))                 * "    "
                sb = sb * @sprintf("%.4e",                                              sharing.weight)                    * "       "
                sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", sharing.differentialRate.Coulomb))   * "   "
                sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", sharing.differentialRate.Babushkin)) * "   "
                println(stream,  sb )
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `DoubleAutoIonization.displayLines(lines::Array{DoubleAutoIonization.Line,1})`  
        ... to display a list of lines, sharings and channels that have been selected due to the prior settings. A neat table 
            of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{DoubleAutoIonization.Line,1})
        #
        # First, print lines and sharings
        nx = 82
        println(" ")
        println("  Selected double-Auger lines & sharings:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.flushleft(54, "Energies (all in )" * TableStrings.inUnits("energy") * ")"; na=5);              
        sb = sb * TableStrings.flushleft(54, "   i -- f       epsilon_1     epsilon_2"; na=5)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
            energy = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
            #
            for  (is, sharing)  in  enumerate(line.sharings)
                if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1)) * "    "
                sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2)) * "    "
                println( sb )
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        #
        # Second, print lines and channles
        nx = 123
        println(" ")
        println("  Selected photo-double ionization lines & channels:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.flushleft(77, "Channels (all energies in " * TableStrings.inUnits("energy") * ")" ; na=5);              
        sb = sb * TableStrings.flushleft(77, "(epsilon, kappa)_1   -->    J^P_x  -->   (epsilon, kappa)_2   -->  J^P_t = J^P_i"; na=5)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #
        for  line in lines
            sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            for (ic, ch) in enumerate(line.sharings[1].channels)
                if  ic == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
                skappa1 = "  " * string(ch.kappa1);    skappa2 = "  " * string(ch.kappa2);    sxsym = "      " * string(ch.xSymmetry)
                sb = sb * " (" * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", ch.energy1)) * ", " * skappa1[end-2:end] * ")    -->"
                sb = sb * sxsym[end-8:end] * "  -->   "
                sb = sb * " (" * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", ch.energy2)) * ", " * skappa2[end-2:end] * ")    -->     "
                sb = sb * string(isym) 
                println( sb )
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `DoubleAutoIonization.displayTotalRates(stream::IO, lines::Array{DoubleAutoIonization.Line,1}, settings::DoubleAutoIonization.Settings)`  
        ... to display all total rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayTotalRates(stream::IO, lines::Array{DoubleAutoIonization.Line,1}, settings::DoubleAutoIonization.Settings)
        #
        # First, print lines and sharings
        nx = 88
        println(stream, " ")
        println(stream, "  Total rates of selected double-Auger lines:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(21)
        sa = sa * TableStrings.center(12, "i--Energy--f"; na=3)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(30, "Cou --   total rate   -- Bab"; na=4)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("rate") * "          " * 
                                          TableStrings.inUnits("rate"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
            energy = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))                 * "      "
            #
            sb = sa * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", line.totalRate.Coulomb))   * "     "
            sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", line.totalRate.Babushkin)) * "   "
            println(stream,  sb )
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

end # module
