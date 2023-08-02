
"""
`module  JAC.InternalRecombination`  
    ... a submodel of JAC that contains all methods for computing internal-recombination rates between some initial 
        (Rydberg-type) and final-state multiplets.
"""
module InternalRecombination

    using  Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, 
                   ..Radial, ..SpinAngular, ..TableStrings
    
    """
    `struct  InternalRecombination.Settings  <:  AbstractProcessSettings`  
        ... defines a type for the details and parameters of computing internal-recombination lines.

        + rydbergShells         ::Array{Shell,1}     ... List of (Rydberg) shells which are to be coupled to the initial levels.
        + printBefore           ::Bool               ... True, if all energies and lines are printed before their evaluation.
        + resonanceEnergyShift  ::Float64            ... An overall energy shift for resonance energies [Hartree].
        + gamma                 ::Float64            ... An overall widths of resonances [Hartree].
        + lineSelection         ::LineSelection      ... Specifies the selected levels, if any.
        + operator              ::AbstractEeInteraction   
            ... Electron-interaction operator that is to be used for evaluating the internal-recombination amplitudes; 
                allowed values are: CoulombInteraction(), BreitInteraction(), ...
    """
    struct Settings  <:  AbstractProcessSettings
        rydbergShells           ::Array{Shell,1}
        printBefore             ::Bool 
        resonanceEnergyShift    ::Float64
        gamma                   ::Float64
        lineSelection           ::LineSelection 
        operator                ::AbstractEeInteraction 
    end 


    """
    `InternalRecombination.Settings()`  ... constructor for the default values of InternalRecombination line computations
    """
    function Settings()
        Settings(Shell[], false, 0., 0.1, LineSelection(), CoulombInteraction())
    end


    """
    `InternalRecombination.Settings(set::InternalRecombination.Settings;`
    
            rydbergShells=..,           printBefore=..,             resonanceEnergyShift=..,    
            gamma=..,                   lineSelection=..,           operator=..)
                        
        ... constructor for modifying the given InternalRecombination.Settings by 'overwriting' the previously selected parameters.
    """
    function Settings(set::InternalRecombination.Settings;    
        rydbergShells::Union{Nothing,Array{Shell,1}}=nothing,    printBefore::Union{Nothing,Bool}=nothing, 
        resonanceEnergyShift::Union{Nothing,Float64}=nothing,    gamma::Union{Nothing,Float64}=nothing,  
        lineSelection::Union{Nothing,LineSelection}=nothing,     operator::Union{Nothing,String}=nothing)  
        
        if  rydbergShells        == nothing  rydbergShellsx        = set.rydbergShells        else  rydbergShellsx        = rydbergShells        end 
        if  printBefore          == nothing  printBeforex          = set.printBefore          else  printBeforex          = printBefore          end 
        if  resonanceEnergyShift == nothing  resonanceEnergyShiftx = set.resonanceEnergyShift else  resonanceEnergyShiftx = resonanceEnergyShift end 
        if  gamma                == nothing  gammax                = set.gamma                else  gammax                = gamma                end 
        if  lineSelection        == nothing  lineSelectionx        = set.lineSelection        else  lineSelectionx        = lineSelection        end 
        if  operator             == nothing  operatorx             = set.operator             else  operatorx             = operator             end 

        Settings( rydbergShellsx, printBeforex, resonanceEnergyShiftx, gammax, lineSelectionx, operatorx)
    end


    # `Base.show(io::IO, settings::InternalRecombination.Settings)`  
    #   ... prepares a proper printout of the variable settings::InternalRecombination.Settings.
    function Base.show(io::IO, settings::InternalRecombination.Settings) 
        println(io, "rydbergShells:            $(settings.rydbergShells)  ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "resonanceEnergyShift:     $(settings.resonanceEnergyShift)  ")
        println(io, "gamma:                    $(settings.gamma)  ")
        println(io, "lineSelection:            $(settings.lineSelection)  ")
        println(io, "operator:                 $(settings.operator)  ")
    end


    """
    `struct  InternalRecombination.Channel`   
        ... defines a type for a InternalRecombination channel to help characterize am initial Rydberg state with 
            one electron at high n.

        + subshell       ::Subshell             ... subshell of the Rydberg electron.
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the overall Rydberg state.
        + amplitude      ::Complex{Float64}     ... Internal-recombination amplitude associated with the given channel.
    """
    struct  Channel
        subshell         ::Subshell 
        symmetry         ::LevelSymmetry
        amplitude        ::Complex{Float64}
    end


    """
    `struct  InternalRecombination.Line`  
        ... defines a type for a InternalRecombination line that may include the definition of sublines and their 
            corresponding amplitudes.

        + initialLevel   ::Level           ... initial-(state) level
        + finalLevel     ::Level           ... final-(state) level
        + deltaEnergy    ::Float64         ... Energy difference of the initial and final levels.
        + rateZ          ::Float64         ... Rate of this line as obtained from the level energies.
        + rate           ::Float64         ... Rate of this line as obtained with shifted resonance energy.
        + channels       ::Array{InternalRecombination.Channel,1}  ... List of internal-recombination channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        deltaEnergy      ::Float64
        rateZ            ::Float64 
        rate             ::Float64 
        channels         ::Array{InternalRecombination.Channel,1}
    end 


    """
    `InternalRecombination.Line()`  
        ... constructor for an empty InternalRecombination line.
    """
    function Line()
        Line(Level(), Level(), 0., 0., 0., InternalRecombination.Channel[])
    end


    # `Base.show(io::IO, line::InternalRecombination.Line)`  ... prepares a proper printout of the variable line::InternalRecombination.Line.
    function Base.show(io::IO, line::InternalRecombination.Line) 
        println(io, "initialLevel:           $(line.initialLevel)  ")
        println(io, "finalLevel:             $(line.finalLevel)  ")
        println(io, "deltaEnergy:            $(line.deltaEnergy)  ")
        println(io, "rateZ:                  $(line.rateZ)  ")
        println(io, "rate:                   $(line.rate)  ")
        println(io, "channels:               $(line.channels)  ")
    end


    """
    `InternalRecombination.computeAmplitudesProperties(line::InternalRecombination.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                                       settings::InternalRecombination.Settings; printout::Bool=true)` 
        ... to compute all amplitudes and properties of the given line; a line::InternalRecombination.Line is returned 
            for which the amplitudes and properties are now evaluated.
    """
    function computeAmplitudesProperties(line::InternalRecombination.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                         settings::InternalRecombination.Settings; printout::Bool=true) 
        newChannels = InternalRecombination.Channel[];   rateZ = 0.;   rate = 0.;   gHalf = settings.gamma / 2.
        # Define a common subshell list for both multiplets
        subshellList = Basics.generate("subshells: ordered list for two bases", line.finalLevel.basis, line.initialLevel.basis)
        Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
        
        for channel in line.channels
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, subshellList)
            newfLevel = Basics.generateLevelWithExtraSubshell(channel.subshell, newfLevel)
            # Generate Rydberg orbital in the field of the initial core
            ## newrLevel = Basics.generateLevelWithExtraElectron(rOrbital, channel.symmetry, newiLevel)
            ## amplitude = InternalRecombination.amplitude(settings.operator, newfLevel, newrLevel, grid, printout=printout)
            amplitude = Complex(1.)
            # Apply energy dominator + gamma
            factorZ   = gHalf / ( (line.initialLevel.energy - line.finalLevel.energy)^2 + gHalf^2 )
            factor    = gHalf / ( line.deltaEnergy^2 + gHalf^2 )
            rateZ     = rateZ + conj(amplitude) * amplitude * factorZ
            rate      = rate  + conj(amplitude) * amplitude * factor
            push!( newChannels, InternalRecombination.Channel(channel.subshell, channel.symmetry, amplitude) )
        end
        rateZ = 2pi* rateZ;   rate = 2pi* rate
        newLine   = InternalRecombination.Line(line.initialLevel, line.finalLevel, line.deltaEnergy, rateZ, rate, newChannels)
        
        return( newLine )
    end


    """
    `InternalRecombination.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                        settings::InternalRecombination.Settings; output=true, printout::Bool=true)`  
        ... to compute the internal-recombination transition amplitudes and rates as requested by the given settings. A list of 
            lines::Array{InternalRecombination.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::InternalRecombination.Settings; output=true, printout::Bool=true)
        println("")
        printstyled("InternalRecombination.computeLines(): The computation of internal-recombination rates starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        lines = InternalRecombination.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    InternalRecombination.displayLines(lines)    end  
        # Calculate all amplitudes and requested properties
        newLines = InternalRecombination.Line[]
        for  line in lines
            newLine = InternalRecombination.computeAmplitudesProperties(line, nm, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        InternalRecombination.displayRates(stdout, newLines, settings)
        InternalRecombination.displayTotalRates(stdout, newLines, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   InternalRecombination.displayRates(iostream, newLines, settings)    
                           InternalRecombination.displayTotalRates(iostream, newLines, settings)        end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `InternalRecombination.determineChannels(finalLevel::Level, initialLevel::Level, settings::InternalRecombination.Settings)`  
        ... to determine a list of internal-recombination Channel for a transitions from the initial to final level and by taking 
            into account the particular settings of for this computation; an Array{InternalRecombination.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::InternalRecombination.Settings)
        channels  = InternalRecombination.Channel[];   
        symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        subshells = Basics.generateSubshellList(settings.rydbergShells)
        for  subsh in subshells
            tSymmetries = AngularMomentum.allowedTotalSymmetries(symi, subsh.kappa)
            @show symi, subsh.kappa, tSymmetries
            for  symt in tSymmetries
                if  symt != symf      continue    end
                push!(channels, InternalRecombination.Channel(subsh, symt, Complex(0.)) )
            end
        end
        return( channels )  
    end


    """
    `InternalRecombination.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::InternalRecombination.Settings)`  
        ... to determine a list of InternalRecombination.Line's for transitions between levels from the initial- and final-state 
            multiplets, and by taking into account the particular selections and settings for this computation; 
            an Array{InternalRecombination.Line,1} is returned. Apart from the level specification, all physical properties are 
            set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::InternalRecombination.Settings)
        lines = InternalRecombination.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    dEnergy = iLevel.energy - fLevel.energy   + settings.resonanceEnergyShift
                    ##x if   energy < 0.01                                                             continue   end
                    channels = InternalRecombination.determineChannels(fLevel, iLevel, settings) 
                    push!( lines, InternalRecombination.Line(iLevel, fLevel, dEnergy, 0., 0., channels) )
                end
            end
        end
        return( lines )
    end


    """
    `InternalRecombination.displayLines(lines::Array{InternalRecombination.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{InternalRecombination.Line,1})
        nx = 150
        println(" ")
        println("  Selected internal-recombination lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(14, "Energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(14, "Delta Energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.flushleft(37, "List of subshells and total symmetries"; na=4)  
        sb = sb * TableStrings.flushleft(37, "subshell (total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.7e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "    "
            sa = sa * @sprintf("%.7e", Defaults.convertUnits("energy: from atomic", line.deltaEnergy))          * "    "
            subshellSymmetryList = Tuple{Subshell,LevelSymmetry}[]
            for  channel in line.channels
                push!( subshellSymmetryList, (channel.subshell, channel.symmetry) )
            end
            sa = sa * TableStrings.subshellSymmetryTupels(80, subshellSymmetryList)
            println( sa )
        end
        println("  ", TableStrings.hLine(nx), "\n")
        #
        return( nothing )
    end


    """
    `InternalRecombination.displayRates(stream::IO, lines::Array{InternalRecombination.Line,1}, settings::InternalRecombination.Settings)`  
        ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is returned 
            otherwise.
    """
    function  displayRates(stream::IO, lines::Array{InternalRecombination.Line,1}, settings::InternalRecombination.Settings)
        nx = 104
        println(stream, " ")
        println(stream, "  Internal recombination rates without and with the resonance shift: \n")
        println(stream, "  Resonance shift [Hartree]:    $(settings.resonanceEnergyShift) ")
        println(stream, "  Gamma [Hartree]:              $(settings.gamma)  \n")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "Energy"   ; na=2);               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(12, "Delta energy"   ; na=2);               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(14, "Rate (dE=0.)"   ; na=2);       
        sb = sb * TableStrings.center(14, TableStrings.inUnits("rate");  na=2)
        sa = sa * TableStrings.center(15, "Rate"           ; na=2);                           
        sb = sb * TableStrings.center(14, TableStrings.inUnits("rate");  na=2)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=2)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "  "
            sx = "     " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.deltaEnergy))
            sa = sa * sx[end-13:end]                                                                            * "  "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.rateZ))                  * "  "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.rate))                   * "  "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `InternalRecombination.displayTotalRates(stream::IO, lines::Array{InternalRecombination.Line,1}, 
                                             settings::InternalRecombination.Settings)`  
        ... to list all total rates for the initial levels. A neat table is printed but nothing is returned 
            otherwise.
    """
    function  displayTotalRates(stream::IO, lines::Array{InternalRecombination.Line,1}, settings::InternalRecombination.Settings)
        nx = 76
        println(stream, " ")
        println(stream, "  Total internal-recombination rates for initial levels without and with the resonance shift: \n")
        println(stream, "  Resonance shift [Hartree]:    $(settings.resonanceEnergyShift) ")
        println(stream, "  Gamma [Hartree]:              $(settings.gamma)  \n")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center( 8, "i-level"; na=2);                         sb = sb * TableStrings.hBlank(10)
        sa = sa * TableStrings.center( 6, "J^P";     na=2);                         sb = sb * TableStrings.hBlank( 8)
        sa = sa * TableStrings.center(12, "Energy";  na=2);               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(18, "Total rate (dE=0.)"; na=0);       
        sb = sb * TableStrings.center(18, TableStrings.inUnits("rate");  na=0)
        sa = sa * TableStrings.center(12, "Total rate"; na=2);                           
        sb = sb * TableStrings.center(12, TableStrings.inUnits("rate");  na=2)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #
        used = falses(length(lines))
        for  (i, line) in enumerate(lines)
            if  used[i]     continue    end
            totalZ = 0.;    total = 0.
            for  j = 1:length(lines)
                if  line.initialLevel.index == lines[j].initialLevel.index      used[j] = true
                    totalZ = totalZ + lines[j].rateZ
                    total  = total  + lines[j].rate 
                end
            end
            sa  = "       ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
            sa = sa * string(line.initialLevel.index) * "         " * string(isym)                              * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", totalZ))                      * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", total))                       * "    "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

end # module

