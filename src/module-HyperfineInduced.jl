
"""
`module  JAC.HyperfineInduced`  
... a submodel of JAC that contains all methods for computing HFS A and B coefficients, hyperfine sublevel representations, etc.
"""
module HyperfineInduced


using Printf, ..AngularMomentum, ..Basics,  ..Defaults, ..Hfs, ..ManyElectron, ..Radial, ..Nuclear, 
              ..TableStrings, ..PhotoEmission



"""
`struct  HyperfineInduced.IJF_Vector`  
    ... defines a type for a IJF-coupled basis vector with given nuclear spin and parity as well as
        an electronic ASF (level).

    + I             ::AngularJ64        ... Nuclear spin I.
    + F             ::AngularJ64        ... Total angular momentum F
    + nuclearParity ::Parity            ... Parity of the nuclear level.
    + levelJ        ::Level             ... Atomic level with well-defined total (electronic) angular momentum J.
    + isomer        ::Nuclear.Isomer    ... Isomeric level of the nucleus.
"""
struct IJF_Vector
    I               ::AngularJ64
    F               ::AngularJ64
    nuclearParity   ::Parity 
    levelJ          ::Level
    isomer          ::Nuclear.Isomer
end 


"""
`HyperfineInduced.IJF_Vector()`  ... constructor for an `empty` instance of the basis IJF_Vector.
"""
function IJF_Vector()
    IJF_Vector(AngularJ64(0), AngularJ64(0), Basics.plus, Level(), Nuclear.Isomer())
end


# `Base.show(io::IO, ijfVector::HyperfineInduced.IJF_Vector)`  ... prepares a proper printout of the variable ijfVector.
function Base.show(io::IO, ijfVector::HyperfineInduced.IJF_Vector) 
    println(io, "I:              $(ijfVector.I)  ")
    println(io, "F:              $(ijfVector.F)  ")
    println(io, "nuclearParity:  $(ijfVector.nuclearParity)  ")
    println(io, "levelJ:         $(ijfVector.levelJ)  ")
    println(io, "isomer:         $(ijfVector.isomer)  ")
end


"""
`struct  HyperfineInduced.IJF_Level`  
    ... defines a type for a IJF-coupled level that can be based on the coupling of different nuclear and electronic levels;
        such levels suitable to model hyperfine-induced transitions of either the nucleus and/or the electronic system.

    + I              ::AngularJ64       ... Total angular momentum I.
    + F              ::AngularJ64       ... Total angular momentum F.
    + M              ::AngularM64       ... Total projection M, only defined if a particular sublevel is referred to.
    + nuclearParity  ::Parity           ... nuclear parity of the level.
    + index          ::Int64            ... Index of this level in its original multiplet, or 0.
    + energy         ::Float64          ... energy of the level, including the nuclear and electronic contributions.
    + eLevel         ::Level            ... electronic level on which this hyperfine level is based; it includes the 
                                            electronic J and parity.
    + basis          ::Array{HyperfineInduced.IJF_Vector,1}      ... set of IJF-coupled basis states for this level.
    + mc             ::Vector{ComplexF64} ... Vector of mixing coefficients w.r.t basis.
"""
struct IJF_Level
    I                ::AngularJ64  
    F                ::AngularJ64
    M                ::AngularM64
    nuclearParity    ::Parity
    index            ::Int64 
    energy           ::Float64
    eLevel           ::Level  
    basis            ::Array{HyperfineInduced.IJF_Vector,1}
    mc               ::Vector{ComplexF64}
end 


"""
`HyperfineInduced.IJF_Level()`  ... constructor for an `empty` instance of IJF_Level.
"""
function IJF_Level()
    IJF_Level(AngularJ64(0), AngularJ64(0), AngularM64(0), Basics.plus, 0, 0., Level(), HyperfineInduced.IJF_Vector[], ComplexF64[])
end


# `Base.show(io::IO, ijfLevel::HyperfineInduced.IJF_Level)`  ... prepares a proper printout of the variable ijfLevel::HyperfineInduced.IJF_Level.
function Base.show(io::IO, ijfLevel::HyperfineInduced.IJF_Level) 
    println(io, "I:              $(ijfLevel.I)  ")
    println(io, "F:              $(ijfLevel.F)  ")
    println(io, "M:              $(ijfLevel.M)  ")
    println(io, "nuclearParity:  $(ijfLevel.nuclearParity)  ")
    println(io, "index:          $(ijfLevel.index)  ")
    println(io, "energy:         $(ijfLevel.energy)  ")
    println(io, "eLevel:         $(ijfLevel.eLevel)  ")
    println(io, "basis:          $(ijfLevel.basis)  ")
    println(io, "mc:             $(ijfLevel.mc)  ")
end


"""
`struct  HyperfineInduced.IJF_Multiplet`  ... defines a type for a multiplet of IJF-coupled levels.

    + name     ::String                                 ... A name associated to the multiplet.
    + levelFs  ::Array{HyperfineInduced.IJF_Level,1}    ... List of IJF-coupled levels.

"""
struct IJF_Multiplet
    name       ::String
    levelFs    ::Array{HyperfineInduced.IJF_Level,1}
end 


"""
`HyperfineInduced.IJF_Multiplet()`  ... constructor for an `empty` instance of HyperfineInduced.IJF_Multiplet.
"""
function IJF_Multiplet()
    IJF_Multiplet("", HyperfineInduced.IJF_Level[])
end


# `Base.show(io::IO, ijfMultiplet::HyperfineInduced.IJF_Multiplet)`  ... prepares a proper printout of this variable.
function Base.show(io::IO, ijfMultiplet::HyperfineInduced.IJF_Multiplet) 
    println(io, "name:           $(ijfMultiplet.name)  ")
    println(io, "levelFs:        $(ijfMultiplet.levelFs)  ")
end


"""
`struct  HyperfineInduced.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing radiative lines.

    + multipoles                ::Array{EmMultipoles,1}     ... Specifies the (radiat. field) multipoles to be included.
    + hfMultipoles              ::Array{EmMultipoles,1}     ... Specifies the multipoles of the hyperfine interaction [M1, E2, E3, ...]
    + gauges                    ::Array{UseGauge,1}         ... Gauges to be included into the computations.
    + iIsomer                   ::Nuclear.Isomer            ... Initial isomeric levels in the HFI transitions.
    + iLevelIndex               ::Int64                     
        ... level index of the initial electronic level for which the expansion in the initial hyperfine basis is done.
    + iAddIndices               ::Array{Int64,1}            
        ... level indices of the initial electronic levels which are included into the initial hyperfine basis.
    + iFvalues                  ::Array{AngularJ64,1}     
        ... List of allowed F-values for which hyperfine levels are included into the initial hyperfine basis.
    + fIsomer                   ::Nuclear.Isomer            ... Final isomeric levels in the HFI transitions.
    + fLevelIndex               ::Int64                     
        ... level index of the final electronic level for which the expansion in the final hyperfine basis is done.
    + fAddIndices               ::Array{Int64,1}            
        ... level indices of the final electronic levels which are included into the final hyperfine basis.
    + fFvalues                  ::Array{AngularJ64,1}     
        ... List of allowed F-values for which hyperfine levels are included into the final hyperfine basis.
    + printBefore               ::Bool                      ... True, if all energies and lines are printed before the computation.
    + photonEnergyShift         ::Float64                   ... An overall energy shift for all photon energies in the HFI transitions.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                  ::Array{EmMultipole,1}
    hfMultipoles                ::Array{EmMultipole,1}
    gauges                      ::Array{UseGauge}
    iIsomer                     ::Nuclear.Isomer   
    iLevelIndex                 ::Int64                     
    iAddIndices                 ::Array{Int64,1}            
    iFvalues                    ::Array{AngularJ64,1}     
    fIsomer                     ::Nuclear.Isomer       
    fLevelIndex                 ::Int64                     
    fAddIndices                 ::Array{Int64,1}            
    fFvalues                    ::Array{AngularJ64,1}     
    printBefore                 ::Bool        
    photonEnergyShift           ::Float64
end 


"""
`HyperfineInduced.Settings()`  ... constructor for the default values of radiative line computations
"""
function Settings()
    Settings(EmMultipole[E1], EmMultipole[M1], UseGauge[Basics.UseCoulomb], Nuclear.Isomer(), 0, Int64[], AngularJ64[], 
             Nuclear.Isomer(), 0, Int64[], AngularJ64[], false, 0.)
end


# `Base.show(io::IO, settings::HyperfineInduced.Settings)`  ... prepares a proper printout of the variable settings::HyperfineInduced.Settings.
function Base.show(io::IO, settings::HyperfineInduced.Settings) 
    println(io, "multipoles:                $(settings.multipoles)  ")
    println(io, "hfMultipoles:              $(settings.hfMultipoles)  ")
    println(io, "gauges:                    $(settings.gauges)  ")
    println(io, "iIsomer:                   $(settings.iIsomer)  ")
    println(io, "iLevelIndex:               $(settings.iLevelIndex)  ")
    println(io, "iAddIndices:               $(settings.iAddIndices)  ")
    println(io, "iFvalues:                  $(settings.iFvalues)  ")
    println(io, "fIsomer:                   $(settings.fIsomer)  ")
    println(io, "fLevelIndex:               $(settings.fLevelIndex)  ")
    println(io, "fAddIndices:               $(settings.fAddIndices)  ")
    println(io, "fFvalues:                  $(settings.fFvalues)  ")
    println(io, "printBefore:               $(settings.printBefore)  ")
    println(io, "photonEnergyShift:         $(settings.photonEnergyShift)  ")
end


"""
`struct  HyperfineInduced.Channel`  
    ... defines a type for a single radiative emission/absorption channel that specifies the multipole, gauge and amplitude.

    + multipole         ::EmMultipole        ... Multipole of the photon emission/absorption.
    + gauge             ::EmGauge            ... Gauge for dealing with the (coupled) radiation field.
    + amplitude         ::Complex{Float64}   ... Amplitude of this multiple channel.
"""
struct Channel 
    multipole           ::EmMultipole
    gauge               ::EmGauge
    amplitude           ::Complex{Float64}
end 


# `Base.show(io::IO, channel::HyperfineInduced.Channel)`  ... prepares a proper printout of the variable channel::HyperfineInduced.Channel.
function Base.show(io::IO, channel::HyperfineInduced.Channel) 
    print(io, "HyperfineInduced.Channel($(channel.multipole), $(channel.gauge), amp = $(channel.amplitude)) ") 
end


"""
`struct  HyperfineInduced.Line`  
    ... defines a type for a (hyperfine-induced) line with corresponding amplitudes. These lines must be treated differently
        for nuclear and electronic transitions. For nuclear transitions, we have I_i, I_f, J_i == J_f, F_i, F_f, while for 
        electronic transitions, we generally have I_i == I_f, J_i, J_f, F_i, F_f, respectively.

    + initialLevel   ::HyperfineInduced.IJF_Level  ... initial level (or state)
    + finalLevel     ::HyperfineInduced.IJF_Level  ... final level (or state)
    + omega          ::Float64             ... Transition frequency of emitted photons from this line.
    + photonRate     ::EmProperty          ... Total rate of this line due to both, nuclear and electronic admixtures of levels.
    + channels       ::Array{HyperfineInduced.Channel,1}  ... List of radiative (photon) channels with decay multipoles.
"""
struct  Line
    initialLevel     ::HyperfineInduced.IJF_Level
    finalLevel       ::HyperfineInduced.IJF_Level
    omega            ::Float64
    photonRate       ::EmProperty
    channels         ::Array{HyperfineInduced.Channel,1}
end 


# `Base.show(io::IO, line::HyperfineInduced.Line)`  ... prepares a proper printout of the variable line::HyperfineInduced.Line.
function Base.show(io::IO, line::HyperfineInduced.Line) 
    println(io, "initialLevel:         $(line.initialLevel)  ")
    println(io, "finalLevel:           $(line.finalLevel)  ")
    println(io, "omega:                $(line.omega)  ")
    println(io, "photonRate:           $(line.photonRate)  ")
    println(io, "channels:             $(line.channels)  ")
end

        
#################################################################################################################################
#################################################################################################################################


"""
`HyperfineInduced.computeAmplitudesProperties(line::HyperfineInduced.Line, grid::Radial.Grid, 
                                              settings::HyperfineInduced.Settings; printout::Bool=true)`  
    ... to compute all amplitudes and properties of the given line; a line::HyperfineInduced.Line is returned for which the 
        amplitudes and properties are now evaluated.
"""
function  computeAmplitudesProperties(line::HyperfineInduced.Line, grid::Radial.Grid, settings::HyperfineInduced.Settings; printout::Bool=true)
    newChannels = HyperfineInduced.Channel[];    rateC = rateB = 0.
    for channel in line.channels
        #
        ## amplitude = HyperfineInduced.amplitudeNuclear("emission", channel.multipole, channel.gauge, line.omega, 
        ##                                               line.finalLevel, line.initialLevel, grid, printout=printout)
        amplitude = ComplexF64(1.5)
        push!( newChannels, HyperfineInduced.Channel( channel.multipole, channel.gauge, amplitude) )
        #
        if       channel.gauge == Basics.Coulomb     rateC = rateC + abs(amplitude)^2
        elseif   channel.gauge == Basics.Babushkin   rateB = rateB + abs(amplitude)^2
        elseif   channel.gauge == Basics.Magnetic    rateB = rateB + abs(amplitude)^2;   rateC = rateC + abs(amplitude)^2
        end
    end
    #     
    # Calculate the photonrate and angular beta if requested 
    wa         = 8pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.F) + 1) 
    photonrate = EmProperty(wa * rateC, wa * rateB)    
    newLine    = HyperfineInduced.Line( line.initialLevel, line.finalLevel, line.omega, photonrate, newChannels)
    #
    ##x @show photonrate, newLine
    
    return( newLine )
end


"""
`HyperfineInduced.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                               settings::HyperfineInduced.Settings; output=true)`  
    ... to compute the hyperfine-induced (radiative) transition amplitudes and all properties as requested by the given settings. 
        The finalMultiplet and initialMultiplet always refer to the electronic systems, whereas all further nuclear information either 
        comes from the nuclear model nm or the settings.nuclearCompound.
        A list of lines::Array{HyperfineInduced.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                       settings::HyperfineInduced.Settings; output=true) 
    # Define a common subshell list for both multiplets
    subshellList = Basics.generate("subshells: ordered list for two bases", finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)
    println("")
    printstyled("HyperfineInduced.computeLines(): The computation of hyperfine-induced transition amplitudes starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    #
    # Generate the hyperfine basis for the initial and final hyperfine levels as well as the levels themselve
    initialHfBasis = HyperfineInduced.generateBasis(initialMultiplet, settings.iLevelIndex, settings.iAddIndices, settings.iFvalues, settings)
    finalHfBasis   = HyperfineInduced.generateBasis(finalMultiplet,   settings.fLevelIndex, settings.fAddIndices, settings.fFvalues, settings)
    
    initialLevels  = HyperfineInduced.determineIJFlevels(initialMultiplet, settings.iIsomer, settings.iLevelIndex, settings.iFvalues, 
                                                         initialHfBasis, grid,   settings)
    finalLevels    = HyperfineInduced.determineIJFlevels(finalMultiplet,   settings.fIsomer, settings.fLevelIndex, settings.fFvalues, 
                                                         finalHfBasis, grid,     settings)
    
    # Generate initial and final hyperfine multiplet
    initialHfMultiplet = HyperfineInduced.IJF_Multiplet("iMultiplet", initialLevels)
    finalHfMultiplet   = HyperfineInduced.IJF_Multiplet("fMultiplet", finalLevels)
    
    if  settings.printBefore
        HyperfineInduced.displayHfMultiplet(stdout, "Initial IJF-coupled multiplet:", initialHfMultiplet)
        HyperfineInduced.displayHfMultiplet(stdout, "Final IJF-coupled multiplet:",   finalHfMultiplet)
    end
    #
    # Determine the HFI transitions/lines and all associated channels
    lines = HyperfineInduced.determineLines(finalHfMultiplet, initialHfMultiplet, settings)
    if  settings.printBefore    HyperfineInduced.displayLines(stdout, lines)    end
    
    # Calculate all transition amplitudes and requested properties
    newLines = HyperfineInduced.Line[]
    for  line in lines
        newLine = HyperfineInduced.computeAmplitudesProperties(line, grid, settings) 
        push!( newLines, newLine)
    end
    
    # Print all results to screen
    HyperfineInduced.displayRates(stdout, newLines, settings)
 
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   HyperfineInduced.displayRates(iostream, newLines, settings)
    end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`HyperfineInduced.determineChannels(finalLevel::IJF_Level, initialLevel::IJF_Level, settings::HyperfineInduced.Settings)`  
    ... to determine a list of HyperfineInduced.Channel for a transitions from the initial to final level and by taking into 
        account the particular settings of for this computation; an Array{HyperfineInduced.Channel,1} is returned.
        All channels only refer to the (decay) multipoles of the radiation field.
"""
function determineChannels(finalLevel::IJF_Level, initialLevel::IJF_Level, settings::HyperfineInduced.Settings)
    channels = HyperfineInduced.Channel[];
    iParity  = initialLevel.nuclearParity * initialLevel.eLevel.parity
    fParity  = finalLevel.nuclearParity   * finalLevel.eLevel.parity
    symi = LevelSymmetry(initialLevel.F, iParity);    symf = LevelSymmetry(finalLevel.F, fParity) 
    
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symi, mp, symf)
            hasMagnetic = false
            for  gauge in settings.gauges
                # Include further restrictions if appropriate
                if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      push!(channels, HyperfineInduced.Channel(mp, Basics.Coulomb,   0.) )
                elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    push!(channels, HyperfineInduced.Channel(mp, Basics.Babushkin, 0.) )  
                elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)                  push!(channels, HyperfineInduced.Channel(mp, Basics.Magnetic,  0.) );
                                                    hasMagnetic = true; 
                end 
            end
        end
    end
    return( channels )  
end


"""
`HyperfineInduced.determineLines(finalMultiplet::IJF_Multiplet, initialMultiplet::IJF_Multiplet, settings::HyperfineInduced.Settings)`  
    ... to determine a list of hyperfine-induced nuclear Line's for transitions between the hyperfine levels from the given 
        IJF-coupled initial- and final-state multiplets and by taking into account the particular selections and 
        settings for this computation;  an Array{HyperfineInduced.Line,1} is returned. Apart from the level specification, 
        all physical properties are set to zero during the initialization process.  
"""
function  determineLines(finalMultiplet::IJF_Multiplet, initialMultiplet::IJF_Multiplet, settings::HyperfineInduced.Settings)
    lines = HyperfineInduced.Line[]
    for  iLevel  in  initialMultiplet.levelFs
        for  fLevel  in  finalMultiplet.levelFs
            omega = iLevel.energy - fLevel.energy
            if  omega > 0.   channels = HyperfineInduced.determineChannels(fLevel, iLevel, settings)   end 
            push!( lines, HyperfineInduced.Line(iLevel, fLevel, omega, EmProperty(0., 0.), channels) )
        end
    end
    
    return( lines )
end


"""
`HyperfineInduced.determineIJFlevels(multiplet::Multiplet, isomer::Nuclear.Isomer, index::Int64, Fvalues::Array{AngularJ64,1},
                                     basis::Array{HyperfineInduced.IJF_Vector,1}, grid::Radial.Grid, settings::HyperfineInduced.Settings)`  
    ... to determine all hyperfine levels that are associated with the given nuclear isomeric level, the electronic level as
        well as the given F-values. It also computes the expansion coefficients of the given hyperfine levels within the given
        hyperfine basis. Usually, a perturbative representation of this level in the given basis is also computed. 
        A list of hyperfine levels::Array{HyperfineInduced.IJF_Level,1} is returned.
"""
function  determineIJFlevels(multiplet::Multiplet, isomer::Nuclear.Isomer, index::Int64, Fvalues::Array{AngularJ64,1},
                             basis::Array{HyperfineInduced.IJF_Vector,1}, grid::Radial.Grid, settings::HyperfineInduced.Settings)
    levels = HyperfineInduced.IJF_Level[];   ndx = 0;    nuclearSymmetry = LevelSymmetry(isomer.spinI, isomer.parity)
    for F  in  Fvalues
        for eLevel  in  multiplet.levels
            I = isomer.spinI;   nuclearParity = isomer.parity;   nuclearEnergy = Defaults.convertUnits("energy: to atomic", isomer.energy)
            M = AngularM64(F.num, F.den);                        totalEnergy   = nuclearEnergy + eLevel.energy               # Set M_F = F
            if  eLevel.index  == index   &&   AngularMomentum.isTriangle(I, eLevel.J, F)
                ndx = ndx + 1;   mc = ComplexF64[];   mbr = 0;    aRenormalized = 0.
                # Compute the mixing coefficients for the given basis
                for  (mbi, bState)  in  enumerate(basis)
                    if  bState.I == I   &&   bState.nuclearParity == nuclearParity  &&   bState.F == F   &&  
                        bState.levelJ.J == eLevel.J                                                     mx = ComplexF64(1.);   mbr = mbi
                    elseif  AngularMomentum.isTriangle(bState.I, bState.levelJ.J, F)    &&
                            AngularMomentum.isTriangle(eLevel.J, I, F) 
                        mx          = ComplexF64(0.)
                        deltaEnergy = nuclearEnergy + eLevel.energy - bState.levelJ.energy
                        for  mp in settings.hfMultipoles
                            if  AngularMomentum.isTriangle(eLevel.J, bState.levelJ.J, AngularJ64(mp.L) )  &&
                                AngularMomentum.isTriangle(bState.I, I, AngularJ64(mp.L))                 &&
                                abs(deltaEnergy)  >  1.0e-10 
                                # Compute the relevant multipole M^M and T^M interaction amplitude
                                MM = Hfs.computeInteractionAmplitudeM(mp, bState.isomer, isomer)
                                TM = Hfs.computeInteractionAmplitudeT(mp, bState.levelJ, eLevel, grid)
                                mx = mx + AngularMomentum.phaseFactor([I, 1, bState.levelJ.J, 1, F]) *
                                          AngularMomentum.Wigner_6j(bState.I, bState.levelJ.J, F, eLevel.J, I, AngularJ64(mp.L)) *
                                          MM * TM  / deltaEnergy                                
                                ##x wa1 = AngularMomentum.phaseFactor([I, 1, bState.levelJ.J, 1, F])
                                ##x wa2 = AngularMomentum.Wigner_6j(bState.I, bState.levelJ.J, F, eLevel.J, I, AngularJ64(mp.L))
                                ##x @show MM, TM, deltaEnergy, wa1, wa2, mx
                            end 
                        end
                    else                                                                                mx = ComplexF64(0.)
                    end
                    push!(mc, mx)
                end
                # Renormalize the contribution of the leading basis state
                for  (mbi, bState)  in  enumerate(basis)
                    if   mbi == mbr   continue   
                    else  aRenormalized = aRenormalized + abs(mc[mbi])^2
                    end
                end
                mc[mbr] = 1.0 - sqrt(aRenormalized)
                #
                push!(levels, HyperfineInduced.IJF_Level(I, F, M, nuclearParity, index, totalEnergy, eLevel, basis, mc) )
            end 
        end 
    end 
    
    println(">>> Generate IJF levels for I^P = $nuclearSymmetry, index = $(index) and F-values = $(Fvalues)  " *
            "==>   No IJF levels = $(length(levels))")
    
    return( levels )
end
    

"""
`HyperfineInduced.displayLines(stream::IO, lines::Array{HyperfineInduced.Line,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all 
        selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{HyperfineInduced.Line,1})
    nx = 95
    println(stream, " ")
    println(stream, "  Selected hyperfine-induced radiative lines:")
    println(stream, " ")

    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--F^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(30, "List of multipoles"; na=4);             sb = sb * TableStrings.hBlank(34)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";       isym = LevelSymmetry( line.initialLevel.F, line.initialLevel.nuclearParity * line.initialLevel.eLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.F,   line.finalLevel.nuclearParity   * line.finalLevel.eLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
        mpGaugeList = Tuple{EmMultipole,EmGauge}[]
        for  i in 1:length(line.channels)
            push!( mpGaugeList, (line.channels[i].multipole, line.channels[i].gauge) )
        end
        sa = sa * TableStrings.multipoleGaugeTupels(50, mpGaugeList)
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, " ")
    #
    return( nothing )
end


"""
`HyperfineInduced.displayHfMultiplet(stream::IO, sa::String, multiplet::HyperfineInduced.IJF_Multiplet)`  
    ... to display the hyperfine levels from a IJF-coupled multiplet; sa is a string that is printed to characterize 
        the multiplet. The results are printed but nothing is returned otherwise.
"""
function  displayHfMultiplet(stream::IO, sa::String, multiplet::HyperfineInduced.IJF_Multiplet)
    nx = 149
    println(stream, " ")
    println(stream, "  $sa")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(8, "level"; na=2);                              sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(8, "F^P";   na=4);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(30, "Mixing coefficients"; na=4);            sb = sb * TableStrings.hBlank(34)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  level in multiplet.levelFs
        sa  = "  ";    sym = LevelSymmetry( level.F, level.nuclearParity * level.eLevel.parity)
        sa = sa * TableStrings.center(8, string(level.index); na=2)
        sa = sa * TableStrings.center(8, string(sym); na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", level.energy)) * "     "
        for  (ic, coeff)  in  enumerate(level.mc)
            if       ic < 12    sa = sa * @sprintf("% .3e", coeff) * " "   
            elseif   ic == 13   sa = sa * "..."
            else     break   
            end  
        end
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx)) 
    println(stream, " ")
    #
    return( nothing )
end


"""
`HyperfineInduced.displayRates(stream::IO, lines::Array{HyperfineInduced.Line,1}, settings::HyperfineInduced.Settings)`  
    ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayRates(stream::IO, lines::Array{HyperfineInduced.Line,1}, settings::HyperfineInduced.Settings)
    nx = 131
    println(stream, " ")
    println(stream, "  Einstein coefficients, transition rates and oscillator strengths for hyperfine-induced transitions:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--F^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "Energy"   ; na=4);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center( 9, "Multipole"; na=0);                         sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(11, "Gauge"    ; na=4);                         sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(26, "A--Einstein--B"; na=3);       
    sb = sb * TableStrings.center(26, TableStrings.inUnits("rate")*"          "*TableStrings.inUnits("rate"); na=2)
    sa = sa * TableStrings.center(12, "Decay widths"; na=3);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        for  ch in line.channels
            sa  = "";      isym = LevelSymmetry( line.initialLevel.F, line.initialLevel.nuclearParity * line.initialLevel.eLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.F,   line.finalLevel.nuclearParity   * line.finalLevel.eLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=4)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
            sa = sa * TableStrings.center(9,  string(ch.multipole); na=4)
            sa = sa * TableStrings.flushleft(11, string(ch.gauge);  na=2)
            chRate =  8pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.F) + 1) * (abs(ch.amplitude)^2) 
            sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to Einstein A",    line, chRate)) * "  "
            sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to Einstein B",    line, chRate)) * "    "
            sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to decay width",   line, chRate)) * "    "
            println(stream, sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`HyperfineInduced.generateBasis(multiplet::Multiplet, index::Int64, addIndices::Array{Int64,1}, Fvalues::Array{AngularJ64,1},
                                settings::HyperfineInduced.Settings)`
    ... generates a complete -- tensor product & IJF-coupled -- from the given nuclear isomeric levels, the selected electronic levels
        (selected by their indices in multiplet) as well as for the given Fvalues. A basis::Array{HyperfineInduced.IJF_Vector,1} 
        is returned.
"""
function  generateBasis(multiplet::Multiplet, index::Int64, addIndices::Array{Int64,1}, Fvalues::Array{AngularJ64,1},
                        settings::HyperfineInduced.Settings)
    # First determine whether one or two isomeric levels are involved into the HFI transitions
    isomers     = [settings.iIsomer];    if   settings.iIsomer != settings.fIsomer   push!(isomers, settings.fIsomer)    end 
    nSymmetries = LevelSymmetry[];   for  isomer in isomers   push!(nSymmetries, LevelSymmetry(isomer.spinI, isomer.parity) )   end
    
    basis    = HyperfineInduced.IJF_Vector[];   
    for F  in  Fvalues
        for  isomer in isomers
            for level  in  multiplet.levels
                if  level.index  in  index  ||   level.index  in  addIndices 
                    ##x @show level.J
                    if   AngularMomentum.isTriangle(isomer.spinI, level.J, F)
                        push!(basis, HyperfineInduced.IJF_Vector(isomer.spinI, F, isomer.parity, level, isomer) )
                    end 
                end
            end 
        end 
    end 
    
    println(">>> Generate IJF basis for nuclear symmetries I^P = $nSymmetries, No. electronic levels = $(length(addIndices)+1) " * 
            "and F-values = $(Fvalues)  ==>  No IJF-coupled basis states = $(length(basis))")
    
    return( basis )
end


end # module


