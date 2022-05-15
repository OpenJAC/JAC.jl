
"""
`module  JAC.ImpactExcitation`  
    ... a submodel of JAC that contains all methods for computing electron impact excitation cross sections and collision strengths.
"""
module ImpactExcitation 

    using Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron, 
                  ..Nuclear, ..Radial, ..SpinAngular, ..TableStrings

    """
    `struct  ImpactExcitationSettings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing electron-impact excitation lines.

        + electronEnergies        ::Array{Float64,1}             ... List of impact-energies of the incoming elecgtrons (in user-defined units).
        + includeBreit            ::Bool                         ... True if the Breit interaction is to be included, and false otherwise.
        + calcCollisionStrength   ::Bool                         ... True, if collision strength need to be calculated, and false otherwise.
        + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.
        + energyShift             ::Float64                      ... An overall energy shift for all transitions |i> --> |f>.
        + maxKappa                ::Int64                        ... Maximum kappa value of partial waves to be included.
        + operator                ::AbstractEeInteraction   
            ... Interaction operator that is to be used for evaluating the e-e interaction amplitudes; allowed values are: 
                CoulombInteraction(), BreitInteraction(), ...
    """
    struct Settings  <:  AbstractProcessSettings
        electronEnergies          ::Array{Float64,1}
        includeBreit              ::Bool 
        calcCollisionStrength     ::Bool
        printBefore               ::Bool 
        lineSelection             ::LineSelection  
        energyShift               ::Float64    
        maxKappa                  ::Int64
        operator                  ::AbstractEeInteraction 
    end 


    """
    `ImpactExcitation.Settings()`  ... constructor for the default values of electron-impact excitation line computations.
    """
    function Settings()
       Settings( Float64[], false, false, false, LineSelection(), 0., 0, CoulombInteraction())
    end


    # `Base.show(io::IO, settings::ImpactExcitation.Settings)`  ... prepares a proper printout of the variable settings::ImpactExcitation.Settings.
    function Base.show(io::IO, settings::ImpactExcitation.Settings) 
        println(io, "electronEnergies:           $(settings.electronEnergies)  ")
        println(io, "includeBreit:               $(settings.includeBreit)  ")
        println(io, "calcCollisionStrength:      $(settings.calcCollisionStrength)  ")
        println(io, "printBefore:                $(settings.printBefore)  ")
        println(io, "lineSelection:              $(settings.lineSelection)  ")
        println(io, "energyShift:                $(settings.energyShift)  ")
        println(io, "maxKappa:                   $(settings.maxKappa)  ")
        println(io, "operator:                   $(settings.operator)  ")
    end


    """
    `struct  ImpactExcitation.Channel`  
        ... defines a type for a electron-impact excitaiton channel to help characterize the incoming and outgoing (continuum) states of 
            many electron-states with a single free electron

        + initialKappa     ::Int64              ... partial-wave of the incoming free electron
        + finalKappa       ::Int64              ... partial-wave of the outgoing free electron
        + symmetry         ::LevelSymmetry      ... total angular momentum and parity of the scattering state
        + initialPhase     ::Float64            ... phase of the incoming partial wave
        + finalPhase       ::Float64            ... phase of the outgoing partial wave
        + amplitude        ::Complex{Float64}   ... Collision amplitude associated with the given channel.
    """
    struct  Channel
        initialKappa       ::Int64 
        finalKappa         ::Int64 
        symmetry           ::LevelSymmetry
        initialPhase       ::Float64
        finalPhase         ::Float64
        amplitude          ::Complex{Float64}
    end


    # `Base.show(io::IO, channel::ImpactExcitation.Channel)`  ... prepares a proper printout of the variable channel::ImpactExcitation.Channel.
    function Base.show(io::IO, channel::ImpactExcitation.Channel) 
        println(io, "initialKappa:       $(channel.initialKappa)  ")
        println(io, "finalKappa:         $(channel.finalKappa)  ")
        println(io, "symmetry:           $(channel.symmetry)  ")
        println(io, "initialPhase:       $(channel.initialPhase)  ")
        println(io, "finalPhase:         $(channel.finalPhase)  ")
        println(io, "amplitude:          $(channel.amplitude)  ")
    end


    """
    `struct  ImpactExcitation.Line`  
        ... defines a type for a electron-impact excitation line that may include the definition of channels and their corresponding amplitudes.

        + initialLevel           ::Level         ... initial- (bound-state) level
        + finalLevel             ::Level         ... final- (bound-state) level
        + initialElectronEnergy  ::Float64       ... energy of the incoming (initial-state) free-electron
        + finalElectronEnergy    ::Float64       ... energy of the outgoing (final-state) free-electron
        + crossSection           ::Float64       ... total cross section of this line
        + collisionStrength      ::Float64       ... total collision strength of this line
        + channels               ::Array{ImpactExcitation.Channel,1}  ... List of ImpactExcitation channels of this line.
    """
    struct  Line
        initialLevel             ::Level
        finalLevel               ::Level
        initialElectronEnergy    ::Float64
        finalElectronEnergy      ::Float64
        crossSection             ::Float64 
        collisionStrength        ::Float64 
        channels                 ::Array{ImpactExcitation.Channel,1}
    end 


    """
    `ImpactExcitation.Line()`  ... 'empty' constructor for an electron-impact excitation line between a specified initial and final level.
    """
    function Line()
        Line(Level(), Level(), 0., 0., 0., 0., ImpactExcitation.Channel[] )
    end


    """
    `ImpactExcitation.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`  
        ... constructor for an electron-impact excitation line between a specified initial and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)
        Line(initialLevel, finalLevel, 0., 0., crossSection, 0., ImpactExcitation.Channel[] )
    end


    # `Base.show(io::IO, line::ImpactExcitation.Line)`  ... prepares a proper printout of the variable line::ImpactExcitation.Line.
    function Base.show(io::IO, line::ImpactExcitation.Line) 
        println(io, "initialLevel:            $(line.initialLevel)  ")
        println(io, "finalLevel:              $(line.finalLevel)  ")
        println(io, "initialElectronEnergy:   $(line.initialElectronEnergy)  ")
        println(io, "finalElectronEnergy:     $(line.finalElectronEnergy)  ")
        println(io, "crossSection:            $(line.crossSection)  ")
        println(io, "collisionStrength:       $(line.collisionStrength)  ")
        println(io, "channels:                $(line.channels)  ")
    end


    """
    `ImpactExcitation.amplitude(kind::AbstractEeInteraction, channel::ImpactExcitation.Channel, cFinalLevel::Level, cInitialLevel::Level, 
                                grid::Radial.Grid; printout::Bool=true)`  
        ... to compute the kind in  CoulombInteraction(), BreitInteraction(), CoulombBreit() electron-impact interaction amplitude 
            <(alpha_f J_f, kappa_f) J_t || O^(e-e, kind) || (alpha_i J_i, kappa_i) J_t>  due to the interelectronic interaction for 
            the given final and initial (continuum) level. A value::ComplexF64 is returned.
    """
    function amplitude(kind::AbstractEeInteraction, channel::ImpactExcitation.Channel, cFinalLevel::Level, cInitialLevel::Level, 
                       grid::Radial.Grid; printout::Bool=true)
        nf = length(cFinalLevel.basis.csfs);      fPartial = Subshell(9,channel.finalKappa)         
        ni = length(cInitialLevel.basis.csfs);    iPartial = Subshell(9,channel.initialKappa)    
        
        if  printout  printstyled("Compute ($kind) e-e matrix of dimension $nf x $ni in the final- and initial-state (continuum) bases " *
                                  "for the transition [$(cInitialLevel.index)- ...] " * 
                                  "and for partial waves $(string(fPartial)[2:end]),  $(string(iPartial)[2:end])... ", color=:light_green)    end
        matrix = zeros(ComplexF64, nf, ni)
        #
        ##x @show cInitialLevel.basis.subshells
        ##x @show cFinalLevel.basis.subshells
        if  cInitialLevel.basis.subshells == cFinalLevel.basis.subshells
            iLevel = cInitialLevel;   fLevel = cFinalLevel
        else
            subshells = Basics.merge(cInitialLevel.basis.subshells, cFinalLevel.basis.subshells)
            iLevel    = Level(cInitialLevel, subshells)
            fLevel    = Level(cFinalLevel,   subshells)
        end
        #
        if      kind in [ CoulombInteraction(), BreitInteraction(), CoulombBreit()]        ## pure V^Coulomb interaction
        #--------------------------------------------------------------------------
            for  r = 1:nf
                for  s = 1:ni
                    if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity != iLevel.parity      continue    end 
                    subshellList = fLevel.basis.subshells
                    opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                    wa   = SpinAngular.computeCoefficients(opa, fLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
                    #
                    me = 0.
                    for  coeff in wa
                        if   kind in [ CoulombInteraction(), CoulombBreit()]    
                            me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, 
                                                    fLevel.basis.orbitals[coeff.a], fLevel.basis.orbitals[coeff.b],
                                                    iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid)   end
                        if   kind in [ BreitInteraction(), CoulombBreit()]    
                            me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, 
                                                    fLevel.basis.orbitals[coeff.a], fLevel.basis.orbitals[coeff.b],
                                                    iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid)   end
                    end
                    matrix[r,s] = me
                end
            end 
            if  printout  printstyled("done. \n", color=:light_green)    end
            amplitude = transpose(fLevel.mc) * matrix * iLevel.mc 
            amplitude = im^Basics.subshell_l(Subshell(102, channel.finalKappa))   * exp( -im*channel.finalPhase )   * 
                        im^Basics.subshell_l(Subshell(101, channel.initialKappa)) * exp(  im*channel.initialPhase ) * amplitude
            @show amplitude
            #
            #
         elseif  kind == "H-E"
        #--------------------
            amplitude = 0.;    error("stop a")
        else    error("stop b")
        end
        
        return( amplitude )
    end


    """
    `ImpactExcitation.computeAmplitudesProperties(line::ImpactExcitation.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                                  settings::ImpactExcitation.Settings; printout::Bool=true)`  
        ... to compute all amplitudes and properties of the given line; a line::ImpactExcitation.Line is returned for which the amplitudes and
            properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::ImpactExcitation.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                          settings::ImpactExcitation.Settings; printout::Bool=true)
        newChannels = ImpactExcitation.Channel[];   contSettings = Continuum.Settings(false, grid.NoPoints-50);   cross = 0.;   coll = 0.
        
        # Define a common subshell list for both multiplets
        subshellList = Basics.generate("subshells: ordered list for two bases", line.finalLevel.basis, line.initialLevel.basis)
        Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
        
        # First determine a common set of continuum orbitals for the incoming and outgoing electron
        ciOrbitals = Dict{Subshell, Orbital}();     ciPhases = Dict{Subshell, Float64}()
        cfOrbitals = Dict{Subshell, Orbital}();     cfPhases = Dict{Subshell, Float64}()
        for channel in line.channels
            # Generate the continuum orbitals if they were not yet generated before
            iSubshell  = Subshell(101, channel.initialKappa)
            if  !haskey(ciOrbitals, iSubshell)
                newiLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
                ciOrbital, iPhase     = Continuum.generateOrbitalForLevel(line.initialElectronEnergy, iSubshell, newiLevel, nm, grid, contSettings)
                ciOrbitals[iSubshell] = ciOrbital;      ciPhases[iSubshell] = iPhase
            end
            
            fSubshell  = Subshell(102, channel.finalKappa)
            if  !haskey(cfOrbitals, fSubshell)
                newfLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
                cfOrbital, fPhase     = Continuum.generateOrbitalForLevel(line.finalElectronEnergy,   fSubshell, newfLevel, nm, grid, contSettings)
                cfOrbitals[fSubshell] = cfOrbital;      cfPhases[fSubshell] = fPhase
            end
        end
        
        for channel in line.channels
            # Generate two continuum orbitals
            newiLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
            newfLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
            iSubshell  = Subshell(101, channel.initialKappa)
            fSubshell  = Subshell(102, channel.finalKappa)
            ciOrbital  = ciOrbitals[iSubshell];     iPhase = ciPhases[iSubshell]
            cfOrbital  = cfOrbitals[fSubshell];     fPhase = cfPhases[fSubshell]
            ##x ciOrbital, iPhase  = Continuum.generateOrbitalForLevel(line.initialElectronEnergy, iSubshell, newiLevel, nm, grid, contSettings)
            ##x cfOrbital, fPhase  = Continuum.generateOrbitalForLevel(line.finalElectronEnergy,   fSubshell, newfLevel, nm, grid, contSettings)
            newiLevel = Basics.generateLevelWithExtraElectron(ciOrbital, channel.symmetry, newiLevel)
            newiLevel = Basics.generateLevelWithExtraSubshell(fSubshell,   newiLevel)
            newfLevel = Basics.generateLevelWithExtraSubshell(iSubshell, newfLevel)
            newfLevel = Basics.generateLevelWithExtraElectron(cfOrbital, channel.symmetry, newfLevel)
            ##x @show newiLevel.basis.subshells
            ##x @show newfLevel.basis.subshells
            newChannel = ImpactExcitation.Channel(channel.initialKappa, channel.finalKappa, channel.symmetry, iPhase, fPhase, 0.)
            #
            amplitude  = ImpactExcitation.amplitude(settings.operator, newChannel, newfLevel, newiLevel, grid, printout=printout)
            coll       = coll + AngularMomentum.bracket([channel.symmetry.J]) * conj(amplitude) * amplitude
            push!( newChannels, ImpactExcitation.Channel(newChannel.initialKappa, newChannel.finalKappa, newChannel.symmetry, iPhase, fPhase, amplitude) )
        end
        # Calculate the electron-impact excitation strength and cross section
        cross   = coll
        newLine = ImpactExcitation.Line( line.initialLevel, line.finalLevel, line.initialElectronEnergy, line.finalElectronEnergy, coll.re, cross.re, newChannels)
        return( newLine )
    end


    """
    `ImpactExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                   settings::ImpactExcitation.Settings; output=true)`  
        ... to compute the electron-impact excitation transition amplitudes and all properties as requested by the given settings. 
            A list of lines::Array{ImpactExcitation.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactExcitation.Settings; 
                           output=true)
        println("")
        printstyled("ImpactExcitation.computePathways(): The computation of electron-impact excitation cross sections starts now ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = ImpactExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    ImpactExcitation.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = ImpactExcitation.Line[]
        for  line in lines
            newLine = ImpactExcitation.computeAmplitudesProperties(line, nm, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        ImpactExcitation.displayResults(newLines)
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end


    """
    `ImpactExcitation.computeMatrix(finalBasis::Basis, initialBasis::Basis, settings::ImpactExcitation.Settings)`  
        ... to compute the transition matrix  (<finalContinuumCSF_r|| V(e-e) ||initialContinuumCSF_s>)  between the CSF_r from the 
            finalContinuumBasis and the CSF_s from the initialContinuumBasis. A (non-quadratic) matrix::Array{Float64,2} with 
            dimensions [length(finalContinuumBasis.csfs) x length(initialContinuumBasis.csfs)] is returned. Note that this transition 
            matrix is typically specific to just one Eimex channel due to the different energies, partial waves and overall symmetry 
            of the scattering states. **Not yet implemented !**
    """
    function computeMatrix(finalBasis::Basis, initialBasis::Basis, settings::ImpactExcitation.Settings)   
        error("Not yet implemented.")
    end



    """
    `ImpactExcitation.determineChannels(finalLevel::Level, initialLevel::Level, settings::ImpactExcitation.Settings)`  
        ... to determine a list of electron-impact excitation Channels for a transitions from the initial to the final level and by 
            taking into account the particular settings of for this computation; an Array{ImpactExcitation.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::ImpactExcitation.Settings)
        channels = ImpactExcitation.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        #
        for  inKappa = -(settings.maxKappa+1):settings.maxKappa
            if  inKappa == 0    continue    end
            symtList = AngularMomentum.allowedTotalSymmetries(symi, inKappa)
            for  symt in symtList
                outKappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
                for  outKappa in outKappaList
                    push!(channels, ImpactExcitation.Channel( inKappa, outKappa, symt, 0., 0.,Complex(0.)) )
                end
            end
        end
        return( channels )  
    end


    """
    `ImpactExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ImpactExcitation.Settings)`  
        ... to determine a list of ImpactExcitation.Line's for transitions between levels from the initial- and final-state multiplets, 
            and by taking into account the particular selections and settings for this computation; an Array{ImpactExcitation.Line,1} is 
            returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ImpactExcitation.Settings)
        energyShift = Defaults.convertUnits("energy: to atomic", settings.energyShift)
        lines = ImpactExcitation.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    for  en in settings.electronEnergies
                        initialElectronEnergy  = Defaults.convertUnits("energy: to atomic", en)
                        finalElectronEnergy    = initialElectronEnergy - (fLevel.energy - iLevel.energy) + energyShift
                        if  finalElectronEnergy < 0    continue   end  
                        channels = ImpactExcitation.determineChannels(fLevel, iLevel, settings) 
                        push!( lines, ImpactExcitation.Line(iLevel, fLevel, initialElectronEnergy, finalElectronEnergy, 0., 0., channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `ImpactExcitation.displayLines(lines::Array{ImpactExcitation.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all 
            selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{ImpactExcitation.Line,1})
        nx = 180
        println(" ")
        println("  Selected electron-impact ionization lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(12, "Energy e_in"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(12, "Energy e_out"; na=4);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.flushleft(57, "List of partial waves and total symmetries"; na=4)  
        sb = sb * TableStrings.flushleft(57, "partial-in [total J^P] partial-out        "; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #
        NoLines = 0   
        for  line in lines
            NoLines = NoLines + 1
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "    "
            kappaInOutSymmetryList = Tuple{Int64,Int64,LevelSymmetry}[]
            for  i in 1:length(line.channels)
                push!( kappaInOutSymmetryList, (line.channels[i].initialKappa, line.channels[i].finalKappa, line.channels[i].symmetry) ) 
            end
            wa = TableStrings.kappaKappaSymmetryTupels(95, kappaInOutSymmetryList)
            sb = sa * wa[1];    println( sb )  
            for  i = 2:length(wa)
                NoLines = NoLines + 1
                sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
            # Avoid long table printouts
            if NoLines > 100  println("\n  ImpactExcitation.displayLines():  A maximum of 100 lines are printed in this table. \n")
               break   
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `ImpactExcitation.displayResults(lines::Array{ImpactExcitation.Line,1})`  
        ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is 
            returned otherwise.
    """
    function  displayResults(lines::Array{ImpactExcitation.Line,1})
        nx = 131
        println(" ")
        println("  Electron-impact excitation cross sections:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(12, "Energy e_in"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(12, "Energy e_out"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(15, "Cross section"; na=3)      
        sb = sb * TableStrings.center(15, TableStrings.inUnits("cross section"); na=3)
        sa = sa * TableStrings.center(15, "Collision strength"; na=3)      
        sb = sb * TableStrings.center(15, " ";                  na=3)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "     "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection))    * "        "
            sa = sa * @sprintf("%.6e", line.collisionStrength)                                                    * "    "
            println(sa)
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

end # module
