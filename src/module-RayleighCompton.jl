
"""
`module  JAC.RayleighCompton`  
... a submodel of JAC that contains all methods for computing elastic Rayleigh and inelastic Compton photon scattering cross sections.
"""
module RayleighCompton

using Printf, ..AngularMomentum, ..AtomicState, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..TableStrings


"""
`struct  RayleighCompton.Settings  <:  AbstractProcessSettings`  ... defines a type for the settings in estimating lastic Rayleigh and inelastic Compton 
                                        photon scattering cross sections

    + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
    + photonEnergies          ::Array{Float64,1}        ... List of photon energies in units defined by the user.
    + green                   ::Array{GreenChannel,1}   ... Precalculated and user-specified Green function of the ion.
    + calcRayleighRaman       ::Bool                    ... Calculate Rayleigh-Raman cross sections and properties (if true).
    + calcAngular             ::Bool                    ... Calculate the angular distribution of the fluorescence photon (if true).
    + calcStokes              ::Bool                    ... Calculate the Stokes parameters of the Rayleigh-scattered photon (if true).
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
    + incidentStokes          ::ExpStokes               ... Stokes parameters of the incident radiation.
    + solidAngles             ::Array{SolidAngle,1}     ... List of solid angles [(theta_1, pho_1), ...].  
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.

"""
struct Settings  <:  AbstractProcessSettings 
    multipoles                ::Array{EmMultipole}
    gauges                    ::Array{UseGauge}
    photonEnergies            ::Array{Float64,1} 
    green                     ::Array{AtomicState.GreenChannel,1}
    calcRayleighRaman         ::Bool 
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
    Settings(EmMultipole[], UseGauge[], GreenChannel[], 0., false, false, false, false, ExpStokes(), SolidAngle[], LineSelection() )
end


# `Base.show(io::IO, settings::RayleighCompton.Settings)`  ... prepares a proper printout of the variable settings::RayleighCompton.Settings.
function Base.show(io::IO, settings::RayleighCompton.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "green:                     (settings.green)  ")
    println(io, "calcRayleighRaman:        $(settings.calcRayleighRaman)  ")
    println(io, "calcAngular:              $(settings.calcAngular)  ")
    println(io, "calcStokes:               $(settings.calcStokes)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "incidentStokes:           $(settings.incidentStokes)  ")
    println(io, "solidAngles:              $(settings.solidAngles)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
end


"""
`struct  RayleighCompton.ReducedChannel`  
    ... defines a type for a Rayleigh-Compton reduced channel/ U^(K, Compton) amplitude to help characterize the scattering of 
        a single photon with well-defined multipolarities; in JAC, a Rayleigh-Compton channel is characterized by a reduced amplitude 
        U^(K, Compton) (level_f, multipole_2, J_nu, omega, multipole_1, level_i).

    + K              ::AngularJ64           ... Rank K of the channel.
    + multipole1     ::EmMultipole          ... Multipole of the incoming photon.
    + multipole2     ::EmMultipole          ... Multipole of the outgoing, scattered photon.
    + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
    + omega1         ::Float64              ... Photon frequency of multipole1.
    + omega2         ::Float64              ... Photon frequency of multipole2.
    + omega          ::Float64              ... Photon frequency in denominator.
    + Jsym           ::LevelSymmetry        ... Symmetry J_nu of the intermediate levels.
    + amplitude      ::Complex{Float64}     ... Two-photon reduced amplitude associated with the given channel.
"""
struct  ReducedChannel
    K                ::AngularJ64
    multipole1       ::EmMultipole
    multipole2       ::EmMultipole
    gauge            ::EmGauge  
    omega1           ::Float64
    omega2           ::Float64
    omega            ::Float64
    Jsym             ::LevelSymmetry
    amplitude        ::Complex{Float64}
end


"""
`struct  RayleighCompton.Channel`  
    ... defines a type for a Rayleigh-Compton reduced channel/ U^(K, Compton) amplitude to help characterize the scattering of 
        a single photon with well-defined multipolarities; in JAC, a Rayleigh-Compton channel is characterized by a reduced amplitude 
        U^(K, Compton) (level_f, multipole_2, J_nu, omega, multipole_1, level_i).

    + isS12          ::Bool                 ... True, if S_12 is meant, and false for S_21
    + Jsym           ::LevelSymmetry        ... Symmetry J_nu of the intermediate levels.
    + multipole1     ::EmMultipole          ... Multipole of the incoming photon.
    + multipole2     ::EmMultipole          ... Multipole of the outgoing, scattered photon.
    + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
    + omega1         ::Float64              ... Photon frequency of multipole1.
    + omega2         ::Float64              ... Photon frequency of multipole2.
    + amplitude      ::Complex{Float64}     ... Two-photon reduced amplitude associated with the given channel.
"""
struct  Channel
    isS12            ::Bool
    Jsym             ::LevelSymmetry
    multipole1       ::EmMultipole
    multipole2       ::EmMultipole
    gauge            ::EmGauge  
    omega1           ::Float64
    omega2           ::Float64
    amplitude        ::Complex{Float64}
end


"""
`RayleighCompton.Channel()`  ... constructor for an `empty` instance of RayleighCompton.Channel.
"""
function Channel()
    Channel(true, LevelSymmetry(), E1, E1, Basics.Coulomb, Basics.EmProperty(0.), 0., 0., ComplexF64(0.) )
end


# `Base.show(io::IO, channel::RayleighCompton.Channel)`  ... prepares a proper printout of the variable channel::RayleighCompton.Channel.
function Base.show(io::IO, channel::RayleighCompton.Channel) 
    println(io, "isS12:        $(channel.isS12)  ")
    println(io, "Jsym:         $(channel.Jsym)  ")
    println(io, "multipole1:   $(channel.multipole1)  ")
    println(io, "multipole2:   $(channel.multipole2)  ")
    println(io, "gauge:        $(channel.gauge)  ")
    println(io, "omega1:       $(channel.omega1)  ")
    println(io, "omega2:       $(channel.omega2)  ")
    println(io, "amplitude:    $(channel.amplitude)  ")
end


"""
`struct  RayleighCompton.Line`  ... defines a type for a Rayleigh-Compton scattering line that may include the definition of channels.

    + initialLevel     ::Level              ... initial-(state) level
    + finalLevel       ::Level              ... final-(state) level
    + inOmega          ::Float64            ... Energy of the incoming photon.
    + outOmega         ::Float64            ... Energy of the outgoing, scattered photon.
    + crossSection     ::EmProperty         ... Cross section for this photoionization.
    + channels         ::Array{RayleighCompton.Channel,1}    ... List of (reduced) RayleighCompton.Channels of this line.
"""
struct  Line
    initialLevel       ::Level
    finalLevel         ::Level
    inOmega            ::Float64
    outOmega           ::Float64
    crossSection       ::EmProperty
    channels           ::Array{RayleighCompton.Channel,1}
end


# `Base.show(io::IO, line::RayleighCompton.Line)`  ... prepares a proper printout of the variable line::RayleighCompton.Line.
function Base.show(io::IO, line::RayleighCompton.Line) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "inOmega:           $(line.inOmega)  ")
    println(io, "outOmega:          $(line.outOmega)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "channels:          $(line.channels)  ")
end

        
#################################################################################################################################
#################################################################################################################################

"""
`RayleighCompton.computeAmplitudesProperties(line::RayleighCompton.Line, grid::Radial.Grid, settings::RayleighCompton.Settings)` 
    ... to compute all amplitudes and properties of the given line; a line::RayleighCompton.Line is returned for which 
        the amplitudes and properties are now evaluated.
"""
function  computeAmplitudesProperties(line::RayleighCompton.Line, grid::Radial.Grid, settings::RayleighCompton.Settings)
    newChannels = RayleighCompton.Channel[]
    for ch in line.channels
        amplitude = ComplexF64(1.0)
        ## amplitude = RayleighCompton.computeReducedAmplitude(ch.K, line.finalLevel, ch.multipole2, ch.omega2, ch.Jsym, ch.omega, 
        ##                                                     ch.multipole1, ch.omega1, line.initialLevel, ch.gauge, grid, settings.green)
        push!( newChannels, RayleighCompton.Channel(ch.isS12, ch.Jsym, ch.multipole1, ch.multipole2, ch.gauge, ch.omega1, ch.omega2, amplitude) )
    end
    cs_Cou       = 0.;   cs_Bab = 0.
    crossSection = EmProperty(cs_Cou, cs_Bab)
    newLine      = RayleighCompton.Line(line.initialLevel, line.finalLevel, line.inOmega, line.outOmega, crossSection, newChannels)
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
        newLine = RayleighCompton.computeAmplitudesProperties(line, grid, settings) 
        push!( newLines, newLine)
    end
    # Display all computed channel amplitudes in turn for all lines
    if  settings.printBefore    for  line  in  newLines  RayleighCompton.displayChannels(stdout, line, settings)    end   end
    # Print all results to screen
    RayleighCompton.displayCrossSections(stdout, newLines, settings)
    if  settings.calcAngular        RayleighCompton.displayAngular(stdout, newLines, settings)             end
    if  settings.calcStokes         RayleighCompton.displayStokesParameters(stdout, newLines, settings)    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   
        if  settings.calcAngular    RayleighCompton.displayAngular(stdout, newLines, settings)             end
        if  settings.calcStokes     RayleighCompton.displayStokesParameters(stdout, newLines, settings)    end
    end
    
    if    output    return( lines )
    else            return( nothing )
    end
end


"""
`RayleighCompton.computeChannelAmplitude(channel::RayleighCompton.Channel, finalLevel::Level, initialLevel::Level, 
                                         grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})` 
    ... to compute the channel amplitudes S_12^(RC) or S_21^(RC); a newChannel::RayleighCompton.Channel is returned. 
"""
function  computeChannelAmplitude(channel::RayleighCompton.Channel, finalLevel::Level, initialLevel::Level, 
                                  grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})
    function  computePoleInterval()
        # Compute from the given parameters the contribution of the pole-interval
        return( 0. )
    end
    
    # Extract the proper Green channel, if available
    hasGreenChannel = false;   gChannel = AtomicState.GreenChannel()
    for  grChannel in greenChannels
        if  grChannel.symmetry == channel.Jsym    hasGreenChannel = true;   gChannel = grChannel;   break   end
    end
    if      hasGreenChannel
    else    error("stop a: Green channel not found for symmetry $(channel.Jsym)")
    end
    
    # Extract whether the S_ik has a pole and, if yes, determine the index leftIdx of the resonance.
    # Analyse of whether there is a sign change in the denominator for the current and next function
    hasPole = false;    leftIdx = 0
    for  (ig, gLevel)  in  enumerate(gChannel.gMultiplet.levels)
        if   ig == 1   continue
        elseif   channel.isS12   &&    sign(initialLevel.energy - gChannel.gMultiplet.levels[ig-1].energy - channel.omega2)  !=  
                                       sign(initialLevel.energy - gLevel.energy - channel.omega2)   
                                       hasPole = true;      leftIdx = ig - 1;   break                              # for S_12
        elseif  !(channel.isS12) &&    sign(initialLevel.energy - gChannel.gMultiplet.levels[ig-1].energy + channel.omega1)  !=  
                                       sign(initialLevel.energy - gLevel.energy + channel.omega1)   
                                       hasPole = true;      leftIdx = ig - 1;   break                              # for S_21
        end
    end
    
    # Extract relevant parameters to (later) compute the contribution of the pole-interval
    if  hasPole      
    end
        
    # Sum over all terms of the Green function channel; distinghuis between S_12 and S_21 integrals
    amplitude  = 0.0im
    for  (ig, gLevel)  in  enumerate(gChannel.gMultiplet.levels)
        if   channel.isS12  # for S_12
            leftMe  = PhotonEmission.amplitude()
            rightMe = PhotonEmission.amplitude()
            me      = leftMe * rightMe / (initialLevel.energy - gLevel.energy - channel.omega2)
        else 
            leftMe  = PhotonEmission.amplitude()
            rightMe = PhotonEmission.amplitude()
            me      = leftMe * rightMe / (initialLevel.energy - gLevel.energy + channel.omega1)
        end
        # Add this term properly if a pole arise
        if      hasPole  &&  ig == leftIdx
            # add 1/2 of this term + contribution from pole interval
            amplitude = amplitude + 0.5 * me + computePoleInterval()
        elseif  hasPole  &&  ig == leftIdx + 1
            amplitude = amplitude + 0.5 * me   # add 1/2 of this term
        else
            amplitude = amplitude + 0.5 * me   # add full term
        end
    end
    
    newChannel = RayleighCompton.Channel(channel.isS12, channel.Jsym, channel.multipole1, channel.multipole2 , channel.gauge, 
                                         channel.omega1, channel.omega2, amplitude)
    
    return( newChannel )
end


"""
`RayleighCompton.determineChannels(finalLevel::Level, initialLevel::Level, inOmega::Float64, outOmega::Float64, settings::RayleighCompton.Settings)`  
    ... to determine a list of Rayleigh-Compton Channel's for a transitions from the initial to final level and by taking 
        into account the particular settings of for this computation; an Array{RayleighCompton.ChannelChannel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, inOmega::Float64, outOmega::Float64, settings::RayleighCompton.Settings)
    channels = RayleighCompton.Channel[];   
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp1 in settings.multipoles
        for  mp2 in settings.multipoles
            symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
            for  symn in symmetries
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp1)[1] == 'E' || string(mp2)[1] == 'E'  &&   gauge == Basics.UseCoulomb
                        push!(channels, Channel(true,  symn, mp1, mp2, Basics.Coulomb,   inOmega, outOmega, ComplexF64(0.)) ) 
                        push!(channels, Channel(false, symn, mp1, mp2, Basics.Coulomb,   inOmega, outOmega, ComplexF64(0.)) ) 
                    elseif string(mp1)[1] == 'E' || string(mp2)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                        push!(channels, Channel(true,  symn, mp1, mp2, Basics.Babushkin, inOmega, outOmega, ComplexF64(0.)) ) 
                        push!(channels, Channel(false, symn, mp1, mp2, Basics.Babushkin, inOmega, outOmega, ComplexF64(0.)) ) 
                    elseif string(mp1)[1] == 'M' && string(mp2)[1] == 'M'
                        push!(channels, Channel(true,  symn, mp1, mp2, Basics.Magnetic,  inOmega, outOmega, ComplexF64(0.)) ) 
                        push!(channels, Channel(false, symn, mp1, mp2, Basics.Magnetic,  inOmega, outOmega, ComplexF64(0.)) ) 
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
                    channels = RayleighCompton.determineChannels(fLevel, iLevel, inOmega, outOmega, settings) 
                    if   length(channels) == 0   continue   end
                    push!( lines, RayleighCompton.Line(iLevel, fLevel, inOmega, outOmega, Basics.EmProperty(0., 0.), channels) )
                end
            end
        end
    end
    return( lines )
end


"""
`RayleighCompton.displayAngular(stream::IO, lines::Array{RayleighCompton.Line,1}, settings::RayleighCompton.Settings)`  
    ... to display a list of lines and the angle-differential cross sections for the selected solid angles. 
        A neat table of all selected transitions, energies and solid angles is printed but nothing is returned otherwise.
"""
function  displayAngular(stream::IO, lines::Array{RayleighCompton.Line,1}, settings::RayleighCompton.Settings)
    nx = 127
    println(stream, " ")
    println(stream, "  Angular distribution of the selected Rayleigh-Compton scattering lines at seleced solid angles: ... to be adapted !!")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "   ";   sb = "   ";   sc = "  "
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=7);                         sb = sb * TableStrings.hBlank(23)
    sa = sa * TableStrings.center(20, "i -- Energy -- f"; na=3);              
    sb = sb * TableStrings.center(20, TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(24, "omega-in   omega-out"; na=2);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=1)
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=6)
    println(stream, sa[1:end]);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * " "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.finalLevel.energy))   * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.inOmega))      * " "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.outOmega))     * "       "
        println(stream, sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`RayleighCompton.displayCrossSections(stream::IO, lines::Array{RayleighCompton.Line,1}, settings::RayleighCompton.Settings)`  
    ... to display a list of lines and their cross sections as they were selected by the settings. A neat table of all selected 
        transitions, energies and cross sections is printed but nothing is returned otherwise.
"""
function  displayCrossSections(stream::IO, lines::Array{RayleighCompton.Line,1}, settings::RayleighCompton.Settings)
    nx = 127
    println(stream, " ")
    println(stream, "  Cross sections of the selected Rayleigh-Compton scattering lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "   ";   sb = "   ";   sc = "  "
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=7);                         sb = sb * TableStrings.hBlank(23)
    sa = sa * TableStrings.center(20, "i -- Energy -- f"; na=3);              
    sb = sb * TableStrings.center(20, TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(24, "omega-in   omega-out"; na=2);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=1)
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=6)
    sa = sa * TableStrings.center(28, "Cou -- cross section -- Bab"; na=0);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("cross section"); na=3)
    sb = sb * TableStrings.center(10, TableStrings.inUnits("cross section"); na=4)
    println(stream, sa[1:end]);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * " "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.finalLevel.energy))   * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.inOmega))      * " "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.outOmega))     * "       "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Coulomb))       * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Babushkin))     * "   "
        println(stream, sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`RayleighCompton.displayChannels(stream::IO, line::RayleighCompton.Line, settings::RayleighCompton.Settings)`  
    ... to display a list of all channels (amplitudes) of the given line as they were selected by the settings. 
        A neat table of all computed channels is printed but nothing is returned otherwise.
"""
function  displayChannels(stream::IO, line::RayleighCompton.Line, settings::RayleighCompton.Settings)
    nx = 160
    println(stream, " ")
    println(stream, "  Channels of the Rayleigh-Compton scattering line:")
    println(stream, " ")
    #
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "   ";   sb = "   ";   sc = "  "
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=7);                         sb = sb * TableStrings.hBlank(23)
    sa = sa * TableStrings.center(20, "i -- Energy -- f"; na=3);
    sc = repeat(" ", length(sa)-1)
    sb = sb * TableStrings.center(20, TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center( 6, "J^p"; na=2);                               sb = sb * TableStrings.hBlank(8)
    sa = sa * TableStrings.center( 6, "Mp1"; na=0);                               sb = sb * TableStrings.hBlank(4)
    sa = sa * TableStrings.center(12, "omega1"; na=2);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center( 6, "Mp2"; na=0);                               sb = sb * TableStrings.hBlank(6)
    sa = sa * TableStrings.center(12, "omega2"; na=2);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center( 6, "Sik"; na=1);                               sb = sb * TableStrings.hBlank(8)
    sa = sa * TableStrings.center(10, "Gauge"; na=4);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(20, "Sik amplitude"; na=0);              
    println(stream, sa[1:end]);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  (ic, ch) in enumerate(line.channels)
        if  ic == 1
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                         fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * " "
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.finalLevel.energy))   * " "
        else  sa = sc
        end 
        #
        sa = sa * " " * TableStrings.center(8, string(ch.Jsym) ) * "  " 
        sa = sa * string(ch.multipole1) * "    " 
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", ch.omega1))      * "    "
        sa = sa * string(ch.multipole2) * "    " 
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", ch.omega2))      * " "
        if  ch.isS12  sa = sa * "   S12   "   else  sa = sa * "   S21   "       end
        sa = sa * "  " * TableStrings.flushleft(10, string(ch.gauge); na=0)
        sa = sa * @sprintf("% .4e", ch.amplitude.re) * "  " 
        sa = sa * @sprintf("% .4e", ch.amplitude.im)
        println(stream, sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`RayleighCompton.displayLines(stream::IO, lines::Array{RayleighCompton.Line,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{RayleighCompton.Line,1})
    nx = 157
    println(stream, " ")
    println(stream, "  Selected Rayleigh-Compton scacttering lines and (reduced) channels:")
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
    sb = sb * TableStrings.flushleft(60, "(multipole_1, Jsym, multipole_2, gauge), ..."; na=1)
    println(stream, sa[1:end]);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * " "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.finalLevel.energy))   * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.inOmega))      * " "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.outOmega))     * "   "
        sc = "         " * string(length(line.channels)) * "      "
        sa = sa * sc[end-12:end]
        mpGaugeList = Tuple{String, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge}[]
        for  channel in  line.channels
            if channel.isS12   sd = "S12"   else   sd = "S21"  end
            push!( mpGaugeList, (sd, channel.multipole1, channel.Jsym, channel.multipole2, channel.gauge) )
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


"""
`RayleighCompton.displayStokesParameters(stream::IO, lines::Array{RayleighCompton.Line,1}, settings::RayleighCompton.Settings)`  
    ... to display a list of lines and the Stokes parameters at the selected solid angles. 
        A neat table of all selected transitions, energies and solid angles is printed but nothing is returned otherwise.
"""
function  displayStokesParameters(stream::IO, lines::Array{RayleighCompton.Line,1}, settings::RayleighCompton.Settings)
    nx = 127
    println(stream, " ")
    println(stream, "  StokesParameters of the selected Rayleigh-Compton scattering lines at seleced solid angles: ... to be adapted !!")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "   ";   sb = "   ";   sc = "  "
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=7);                         sb = sb * TableStrings.hBlank(23)
    sa = sa * TableStrings.center(20, "i -- Energy -- f"; na=3);              
    sb = sb * TableStrings.center(20, TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(24, "omega-in   omega-out"; na=2);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=1)
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=6)
    println(stream, sa[1:end]);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * " "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.finalLevel.energy))   * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.inOmega))      * " "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.outOmega))     * "       "
        println(stream, sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
