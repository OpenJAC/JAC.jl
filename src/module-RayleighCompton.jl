
"""
`module  JAC.RayleighCompton`  
    ... a submodel of JAC that contains all methods for computing elastic Rayleigh and inelastic Compton photon scattering cross sections.
"""
module RayleighCompton

    using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..TableStrings


    """
    `struct  RayleighCompton.Settings`  ... defines a type for the settings in estimating lastic Rayleigh and inelastic Compton 
                                            photon scattering cross sections

        + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
        + photonEnergies          ::Array{Float64,1}        ... List of photon energies in units defined by the user.
        + green                   ::Float64                 ... Precalculated and user-specified Green function of the ion.
        + calcAngular             ::Bool                    ... Calculate the angular distribution of the fluorescence photon.
        + calcStokes              ::Bool                    ... Calculate the Stokes parameters of the Rayleigh-scattered photon.
        + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
        + incidentStokes          ::ExpStokes               ... Stokes parameters of the incident radiation.
        + solidAngles             ::Array{SolidAngle,1}     ... List of solid angles [(theta_1, pho_1), ...].  
        + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.

    """
    struct Settings 
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        photonEnergies            ::Array{Float64,1} 
        green                     ::Float64
        calcAngular               ::Bool
        calcStokes                ::Bool
        printBefore               ::Bool 
        incidentStokes            ::ExpStokes 
        solidAngles               ::Array{SolidAngle,1}
        lineSelection             ::LineSelection
    end 


    """
    `RayleighCompton.Settings()`  ... constructor for the default values of Rayleigh-Compton photon-scattering estimates.
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], Float64[], 0., false, false, false, Basics.ExpStokes(), SolidAngle[], LineSelection() )
    end


    # `Base.show(io::IO, settings::RayleighCompton.Settings)`  ... prepares a proper printout of the variable settings::RayleighCompton.Settings.
    function Base.show(io::IO, settings::RayleighCompton.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        println(io, "photonEnergies:           $(settings.photonEnergies)  ")
        println(io, "green:                    $(settings.green)  ")
        println(io, "calcAngular:              $(settings.calcAngular)  ")
        println(io, "calcStokes:               $(settings.calcStokes)  ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "incidentStokes:           $(settings.incidentStokes)  ")
        println(io, "solidAngles:              $(settings.solidAngles)  ")
        println(io, "lineSelection:            $(settings.lineSelection)  ")
    end


    """
    `struct  RayleighCompton.Channel`  
        ... defines a type for a Rayleigh-Compton channel to help characterize the scattering of a single photon with well-defined multipolarities;
            in JAC, a Rayleigh-Compton channel is characterized by a reduced amplitude 
            U^(K, Compton) (level_f, multipole_2, J_nu, omega, multipole_1, level_i).

        + K              ::AngularJ64           ... Rank K of the channel.
        + multipole1     ::EmMultipole          ... Multipole of the incoming photon.
        + multipole2     ::EmMultipole          ... Multipole of the outgoing, scattered photon.
        + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
        + Jsym           ::LevelSymmetry        ... Symmetry J_nu of the intermediate levels.
        + amplitude      ::Complex{Float64}     ... Two-photon reduced amplitude associated with the given channel.
   """
    struct  Channel
        K                ::AngularJ64
        multipole1       ::EmMultipole
        multipole2       ::EmMultipole
        gauge            ::EmGauge  
        Jsym             ::LevelSymmetry
        amplitude        ::Complex{Float64}
    end


    """
    `struct  RayleighCompton.Line`  ... defines a type for a Rayleigh-Compton scattering line that may include the definition of channels.

        + initialLevel     ::Level              ... initial-(state) level
        + finalLevel       ::Level              ... final-(state) level
        + inPhotonEnergy   ::Float64            ... Energy of the incoming photon.
        + outPhotonEnergy  ::Float64            ... Energy of the outgoing, scattered photon.
        + crossSection     ::EmProperty         ... Cross section for this photoionization.
        + channels         ::Array{RayleighCompton.Channel,1}  ... List of RayleighCompton.Channels of this line.
    """
    struct  Line
        initialLevel       ::Level
        finalLevel         ::Level
        inPhotonEnergy     ::Float64
        outPhotonEnergy    ::Float64
        crossSection       ::EmProperty
        channels           ::Array{RayleighCompton.Channel,1}
    end


    # `Base.show(io::IO, line::RayleighCompton.Line)`  ... prepares a proper printout of the variable line::RayleighCompton.Line.
    function Base.show(io::IO, line::RayleighCompton.Line) 
        println(io, "initialLevel:      $(line.initialLevel)  ")
        println(io, "finalLevel:        $(line.finalLevel)  ")
        println(io, "inPhotonEnergy:    $(line.inPhotonEnergy)  ")
        println(io, "outPhotonEnergy:   $(line.outPhotonEnergy)  ")
        println(io, "crossSection:      $(line.crossSection)  ")
        println(io, "channels:          $(line.channels)  ")
    end


    """
    `RayleighCompton.computeAmplitudesProperties(line::RayleighCompton.Line, grid::Radial.Grid, settings::RayleighCompton.Settings)` 
        ... to compute all amplitudes and properties of the given line; a line::RayleighCompton.Line is returned for which 
            the amplitudes and properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::RayleighCompton.Line, grid::Radial.Grid, settings::RayleighCompton.Settings)
        newChannels = RayleighCompton.Channel[]
        for channel in line.channels
            amplitude = RayleighCompton.computeReducedAmplitudeEmission(channel.K, line.finalLevel, line.outPhotonEnergy, channel.multipole2, 
                                channel.Jsym, line.inPhotonEnergy, channel.omega1, channel.multipole1, line.initialLevel, channel.gauge, 
                                grid, settings.greenChannels)
            push!( newChannels, RayleighCompton.Channel(channel.K, channel.multipole1, channel.multipole2, channel.gauge, channel.Jsym, amplitude) )
        end
        cs_Cou       = 0.;   cs_Bab = 0.
        crossSection = EmProperty(cs_Cou, cs_Bab)
        newLine      = RayleighCompton.Line(line.initialLevel, line.finalLevel, line.inPhotonEnergy, line.outPhotonEnergy, crossSection, newChannels)
        return( newLine )
    end


    """
    `RayleighCompton.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::RayleighCompton.Settings; output=true)`  
        ... to compute the Rayleigh-Compton transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{RayleighCompton.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::RayleighCompton.Settings; output=true)
        println("")
        printstyled("RayleighCompton.computeLines(): The computation of Rayleigh-Compton reduced amplitudes starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = RayleighCompton.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected pathways before the computations start
        if  settings.printBefore    RayleighCompton.displayLines(stdout, lines)    end
        # Calculate all reduced amplitudes and requested properties
        newLines = RayleighCompton.Line[]
        for  line in lines
            newLine = RayleighCompton.computeAmplitudesProperties(line, grid) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        RayleighCompton.displayCrossSections(stdout, newLines, settings)
        if  settings.calcAngular        RayleighCompton.displayAngular(stdout, newLines, settings)             end
        if  settings.calcStokes         RayleighCompton.displayStokesParameters(stdout, newLines, settings)    end
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   
            PhotoExcitationFluores.displayCrossSections(iostream, newLines, settings)     
            if  settings.calcAngular    RayleighCompton.displayAngular(stdout, newLines, settings)             end
            if  settings.calcStokes     RayleighCompton.displayStokesParameters(stdout, newLines, settings)    end
        end
        
        if    output    return( lines )
        else            return( nothing )
        end
    end
    

    """
    `RayleighCompton.determineChannels(finalLevel::Level, initialLevel::Level, settings::RayleighCompton.Settings)`  
        ... to determine a list of Rayleigh-Compton Channel's for a transitions from the initial to final level and by taking 
            into account the particular settings of for this computation; an Array{RayleighCompton.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::RayleighCompton.Settings)
        channels = RayleighCompton.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp1 in settings.multipoles
            for  mp2 in settings.multipoles
                symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                Klist       = oplus(symf.J, symi.J)
                for  symn in symmetries
                    for  gauge in settings.gauges
                        # Include further restrictions if appropriate
                        if     string(mp1)[1] == 'E' || string(mp2)[1] == 'E'  &&   gauge == Basics.UseCoulomb
                            for K in Klist  push!(channels, RayleighCompton.Channel(K, mp1, mp2, Basics.Coulomb, symn, 0.) )     end 
                        elseif string(mp1)[1] == 'E' || string(mp2)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                            for K in Klist  push!(channels, RayleighCompton.Channel(K, mp1, mp2, Basics.Babushkin, symn, 0.) )   end
                        elseif string(mp1)[1] == 'M' && string(mp2)[1] == 'M'
                            for K in Klist  push!(channels, RayleighCompton.Channel(K, mp1, mp2, Basics.Magnetic, symn, 0.) )    end
                        end
                    end
                end 
            end
        end
        return( channels )  
    end


    """
    `RayleighCompton.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RayleighCompton.Settings)`  
        ... to determine a list of RayleighCompton Line's for elastic scattering or Raman-type transitions between the levels from 
            the given initial- and final-state multiplets and by taking into account the particular selections and settings for 
            this computation; an Array{RayleighCompton.Line,1} is returned. Apart from the level specification, all physical 
            properties are set to zero during the initialization process.  
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RayleighCompton.Settings)
        lines = RayleighCompton.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    for  omega in settings.photonEnergies
                        inOmega  = Defaults.convertUnits("energy: to atomic", omega)    
                        outOmega = inOmega - (fLevel.energy - iLevel.energy);       if  outOmega <= 0.    continue   end  
                        channels = RayleighCompton.determineChannels(fLevel, iLevel, settings) 
                        if   length(channels) == 0   continue   end
                        push!( lines, RayleighCompton.Line(iLevel, fLevel, inOmega, outOmega, EmProperty(0., 0.), channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `RayleighCompton.displayLines(stream::IO, lines::Array{RayleighCompton.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(stream::IO, lines::Array{RayleighCompton.Line,1})
        nx = 157
        println(stream, " ")
        println(stream, "  Selected two-photon emission lines with given photon splitting:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "   ";   sb = "   ";   sc = "  "
        sa = sa * TableStrings.center(16, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(16, "i--J^P--f"; na=7);                         sb = sb * TableStrings.hBlank(23)
        sa = sa * TableStrings.center(20, "i -- Energy -- f"; na=3);              
        sb = sb * TableStrings.center(20, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(24, "omega-in   omega-out"; na=0);              
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=1)
        sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(11, "No channels"; na=2);                       sb = sb * TableStrings.hBlank(11)           
        sa = sa * TableStrings.flushleft(60, "List of multipoles & intermediate level symmetries";   na=1)            
        sb = sb * TableStrings.flushleft(60, "(K-rank, multipole_1, Jsym, multipole_2, gauge), ..."; na=1)
        println(stream, sa[1:end]);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * " "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.finalLevel.energy))   * "   "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy))      * " "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.outPhotonEnergy))     * "   "
            sc = "         " * string(length(line.channels)) * "      "
            sa = sa * sc[end-12:end]
            mpGaugeList = Tuple{AngularJ64, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge}[]
            for  channel in  line.channels
                push!( mpGaugeList, (channel.K, channel.multipole1, channel.Jsym, channel.multipole2, channel.gauge) )
            end
            wa = TableStrings.twoPhotonGaugeTupels(60, mpGaugeList)
            if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
            for  i = 1:length(wa)
                sb = TableStrings.hBlank( length(sa) );    sb = sb * wa[i];    println(stream, sb )
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

end # module
