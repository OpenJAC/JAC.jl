
"""
`module  JAC.MultiPhotonDeExcitation`  ... a submodel of JAC that contains all methods for computing multi-photon excitation and decay rates;
                                           it is using JAC, JAC.Radial.
"""
module MultiPhotonDeExcitation

    using Printf, JAC.BasicTypes, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0


    """
    `struct  MultiPhotonDeExcitation.Settings`  ... defines a type for the settings in estimating multi-photon excitation and decay rates.

        + NoPhotons               ::Int64                        ... Number of photons in the multi-photon ionization
        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).

    """
    struct Settings 
        NoPhotons                 ::Int64
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        printBeforeComputation    ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1} 
    end 


    """
    `JAC.MultiPhotonDeExcitation.Settings()`  ... constructor for the default values of multi-photon excitation and decay rates.
    """
    function Settings()
        Settings(0, EmMultipole[], UseGauge[], false, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::MultiPhotonDeExcitation.Settings)`  
    #		... prepares a proper printout of the variable settings::MultiPhotonDeExcitation.Settings.
     function Base.show(io::IO, settings::MultiPhotonDeExcitation.Settings) 
        println(io, "NoPhotons:                $(settings.NoPhotons)  ")
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
    end


    """
    `struct  MultiPhotonDeExcitation.Channel`  ... defines a type for a multi-photon channel to help characterize multi-photon 
             (single-electron) ionization with well-defined multipolarities.

        + multipoles     ::Array{EmMultipole,1}   ... Multipoles of all N incoming/outgoing photons.
        + gauge          ::UseGauge               ... Gauge for dealing with the (coupled) radiation field.
        + amplitude      ::Complex{Float64}       ... multi-photon ionization amplitude associated with the given channel.
   """
    struct  Channel
        multipoles       ::Array{EmMultipole,1}
        gauge            ::UseGauge 
        amplitude        ::Complex{Float64}
    end


    """
    `struct  MultiPhotonDeExcitation.Line`  ... defines a type for a multi-photon ionization line that may include the definition of channels.

        + initialLevel     ::Level                  ... initial-(state) level
        + finalLevel       ::Level                  ... final-(state) level
        + photonEnergy     ::Float64                ... Energy of the incoming photons; all photons are assumed to have equal energy.
        + totalRate        ::EmProperty             ... Cross section for this multi-photon ionization.
        + hasChannels      ::Bool                   ... Determines whether the individual (sub-) channels are defined in terms of their 
                                                        multipolarities, etc., or not.
        + channels         ::Array{MultiPhotonDeExcitation.Channel,1}  ... List of MultiPhotonDeExcitation.Channels of this line.
    """
    struct  Line
        initialLevel       ::Level
        finalLevel         ::Level
        photonEnergy       ::Float64
        totalRate          ::EmProperty
        hasChannels        ::Bool
        channels           ::Array{MultiPhotonDeExcitation.Channel,1}
    end


    # `Base.show(io::IO, line::MultiPhotonDeExcitation.Line)`  ... prepares a proper printout of the variable line::MultiPhotonDeExcitation.Line.
    function Base.show(io::IO, line::MultiPhotonDeExcitation.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "photonEnergy:      $(line.photonEnergy)  ")
        println(io, "totalRate:         $(line.totalRate)  ")
        println(io, "hasChannels:       $(line.hasChannels)  ")
        println(io, "channels:          $(line.channels)  ")
    end


    """
    `JAC.MultiPhotonDeExcitation.computeAmplitudesProperties(line::MultiPhotonDeExcitation.Line, grid::Radial.Grid, 
                                                             settings::MultiPhotonDeExcitation.Settings)`  ... to compute all amplitudes and 
         properties of the given line; a line::Einstein.Line is returned for which the amplitudes and properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::MultiPhotonDeExcitation.Line, grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings)
        global JAC_counter
        newChannels = MultiPhotonDeExcitation.Channel[]
        for channel in line.channels
            # matrix    = JAC.MultiPhotonDeExcitation.computeMatrix(channel.multipole, channel.gauge, line.omega, line.finalLevel.basis, 
            #                                                                                         line.initialLevel.basis, grid, settings)
            # amplitude = line.finalLevel.mc * matrix * line.initialLevel.mc 
            JAC_counter = JAC_counter + 1
            if   JAC_counter < 20   println("MultiPhotonDeExcitation.computeAmplitudesProperties-aa: warning ... amplitude always -1.")    end
            amplitude = 1.0 
            push!( newChannels, MultiPhotonDeExcitation.Channel( channel.multipoles, channel.gauge, amplitude) )
        end
        # Calculate the photonrate and angular beta if requested
        JAC_counter = JAC_counter + 1
        if   JAC_counter < 20   println("MultiPhotonDeExcitation.computeAmplitudesProperties-ab: warning ... rate always -1.")    end
        totalRate = EmProperty(-1., -1.)
        line = MultiPhotonDeExcitation.Line( line.initialLevel, line.finalLevel, line.photonEnergy, totalRate, true, newChannels)
        return( line )
    end



    """
    `JAC.MultiPhotonDeExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                              settings::MultiPhotonDeExcitation.Settings; output=true)` ... to compute the multiphoton transition 
         amplitudes and all properties as requested by the given settings. A list of lines::Array{MultiPhotonDeExcitation.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                           settings::MultiPhotonDeExcitation.Settings; output=true)
        println("")
        printstyled("JAC.MultiPhotonDeExcitation.computeLines(): The computation of multiphoton transition amplitudes starts now ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = JAC.MultiPhotonDeExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.MultiPhotonDeExcitation.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = MultiPhotonDeExcitation.Line[]
        for  line in lines
            newLine = JAC.MultiPhotonDeExcitation.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        JAC.MultiPhotonDeExcitation.displayRates(stdout, lines)
        printSummary, iostream = Constants.give("summary flag/stream")
        if  printSummary    JAC.MultiPhotonDeExcitation.displayRates(iostream, lines)   end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `JAC.MultiPhotonDeExcitation.determineChannels(finalLevel::Level, initialLevel::Level, settings::MultiPhotonDeExcitation.Settings)`  
         ... to determine a list of MultiPhotonDeExcitation.Channel for a transitions from the initial to final level and by taking into account 
         the particular settings of for this computation; an Array{MultiPhotonDeExcitation.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::MultiPhotonDeExcitation.Settings)
        global JAC_counter
        channels = MultiPhotonDeExcitation.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        mpListList = Array{JAC.EmMultipole,1}[];    morePhotons = true;    NoPhotons = 1
        for  mp in settings.multipoles    push!( mpListList, [mp] )    end
        if  NoPhotons + 1 > settings.NoPhotons    error("stop a")      end
        while  morePhotons
            NoPhotons     = NoPhotons + 1
            newMpListList = Array{JAC.EmMultipole,1}[]
            for  mpList in mpListList
                for  mp in settings.multipoles
                    mppList = deepcopy(mpList)
                    push!( newMpListList, push!(mppList, mp) )
                end
            end
            mpListList = deepcopy(newMpListList)
            if  NoPhotons + 1 > settings.NoPhotons    morePhotons = false   end
        end
        #
        for  gauge in settings.gauges
            for  mpList in mpListList
                # Cycle if list of multipoles cannot connect the initial and final levels for parity reasons
                JAC_counter = JAC_counter + 1
                if   JAC_counter < 5   println("MultiPhotonDeExcitation.determineChannels-aa: warning ... no test of parities yet.")  end
                push!( channels, MultiPhotonDeExcitation.Channel(mpList, gauge, 0.) )
            end
        end
        ##x println("MultiPhotonDeExcitation.determineChannels-ad:  channels = $channels")
        return( channels )  
    end


    """
    `JAC.MultiPhotonDeExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::MultiPhotonDeExcitation.Settings)`
         ... to determine a list of MultiPhotonDeExcitation Line's for transitions between the levels from the given initial- and final-state 
         multiplets and by taking into account the particular selections and settings for this computation; an Array{MultiPhotonDeExcitation.Line,1}
         is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.  
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::MultiPhotonDeExcitation.Settings)
        if    settings.selectLines    selectLines   = true;   selectedLines = JAC.determine("selected lines", settings.selectedLines)
        else                          selectLines   = false
        end
    
        lines = MultiPhotonDeExcitation.Line[]
        for  i = 1:length(initialMultiplet.levels)
            for  f = 1:length(finalMultiplet.levels)
                if  selectLines  &&  !(haskey(selectedLines, (i,f)) )    continue   end
                ##x println("MultiPhotonDeExcitation.determineLines-aa: angular i = $i, f = $f")
                omega    = abs( initialMultiplet.levels[i].energy - finalMultiplet.levels[f].energy)

                channels = JAC.MultiPhotonDeExcitation.determineChannels(finalMultiplet.levels[f], initialMultiplet.levels[i], settings) 
                push!( lines, MultiPhotonDeExcitation.Line(initialMultiplet.levels[i], finalMultiplet.levels[f], omega, EmProperty(0., 0.), 
                                                           true, channels) )
            end
        end
        return( lines )
    end


    """
    `JAC.MultiPhotonDeExcitation.displayLines(lines::Array{Einstein.Line,1})`  ... to display a list of lines and channels that have been 
         selected due to the prior settings. A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{MultiPhotonDeExcitation.Line,1})
        println(" ")
        println("  Selected multi-photon excitation or decay lines:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(165))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.flushleft(70, "List of multipoles"; na=4)            
        sb = sb * JAC.TableStrings.flushleft(70, "Used gauge [multipolo_1, multipole_2, ...]"; na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(165)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", Constants.convert("energy: from atomic", line.photonEnergy)) * "    "
            mpGaugeList = Tuple{UseGauge, Array{JAC.EmMultipole,1}}[]
            for  i in 1:length(line.channels)
                push!( mpGaugeList, (line.channels[i].gauge, line.channels[i].multipoles) )
            end
            wa = JAC.TableStrings.gaugeMultipolesTupels(105, mpGaugeList)
            if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
            for  i = 2:length(wa)
                sb = JAC.TableStrings.hBlank( length(sa) );    sb = sb * wa[i];    println( sb )
            end
        end
        println("  ", JAC.TableStrings.hLine(165))
        #
        return( nothing )
    end


    """
    `JAC.MultiPhotonDeExcitation.displayRates(stream::IO, lines::Array{MultiPhotonDeExcitation.Line,1})`  ... to display all results, energies, rates, etc. 
         of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayRates(stream::IO, lines::Array{MultiPhotonDeExcitation.Line,1})
        println(stream, " ")
        println(stream, "  Multi-photon excitation or decay lines:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(95))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=3)
        sa = sa * JAC.TableStrings.center(32, "Cou -- Rate -- Bab"; na=4);              
        sb = sb * JAC.TableStrings.center(32, JAC.TableStrings.inUnits("rate") * "           " * JAC.TableStrings.inUnits("rate"); na=4)

        println(stream, sa);    println(stream, b);    println(stream, "  ", JAC.TableStrings.hLine(95)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", Constants.convert("energy: from atomic", line.photonEnergy))         * "    "
            sa = sa * @sprintf("%.8e", Constants.convert("rate: from atomic",   line.totalRate.Coulomb))    * "    "
            sa = sa * @sprintf("%.8e", Constants.convert("rate: from atomic",   line.totalRate.Babushkin))  * "    "
            println(stream, sa )
        end
        println(stream, "  ", JAC.TableStrings.hLine(95))
        #
        return( nothing )
    end

end # module
