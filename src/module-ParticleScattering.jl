
"""
`module  JAC.ParticleScattering`  
    ... a submodel of JAC that contains all methods for computing scattering amplitudes and cross sections
        with regard to "transitions" between some initial and final-state multiplets. It covers different types
        of particles (electrons, protons, ...), different types of scattering processes (elastic, ...), 
        different types of incoming beams (plane-wave, Bessel, Laguerre-Gauss, ...) and different treatments
        (nonrelativistic, relativistic, ...).
"""
module ParticleScattering

    using  Printf, GSL,
           ..AngularMomentum, ..Basics, ..Beam, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, 
           ..Radial, ..SpinAngular, ..TableStrings

    
    """
    `abstract type ParticleScattering.AbstractProcessType` 
        ... defines an abstract type to distinguish different scattering processes of particles with atoms and ions; see also:
        
        + struct ParticleScattering.ElasticElectronNR    ... to model the elastic electron scattering.
        + struct ParticleScattering.InelasticElectronNR  ... to model the inelastic electron scattering (not yet).
    """
    abstract type  AbstractProcessType                                            end
    struct   ElasticElectronNR      <:  ParticleScattering.AbstractProcessType    end
    struct   InelasticElectronNR    <:  ParticleScattering.AbstractProcessType    end
        
    
    """
    `struct  ParticleScattering.Settings  <:  AbstractProcessSettings`  
        ... defines a type for the details and parameters in calculating the scattering amplitudes and cross sections
            for selected lines.

        + processType         ::ParticleScattering.AbstractProcessType
        + beamType            ::Beam.AbstractBeamType
        + polarization        ::Basics.AbstractPolarization
        + impactEnergies      ::Array{Float64,1}
        + polarThetas         ::Array{Float64,1}
        + polarPhis           ::Array{Float64,1}
        + printBefore         ::Bool               ... True, if all energies and lines are printed before their evaluation.
        + lineSelection       ::LineSelection      ... Specifies the selected levels, if any.
        + maxPartialWave      ::Int64              ... maximum partial wave (l_max or kappa_max)
    """
    struct Settings  <:  AbstractProcessSettings
        processType           ::ParticleScattering.AbstractProcessType
        beamType              ::Beam.AbstractBeamType
        polarization          ::Basics.AbstractPolarization
        impactEnergies        ::Array{Float64,1}
        polarThetas           ::Array{Float64,1}
        polarPhis             ::Array{Float64,1}
        printBefore           ::Bool 
        lineSelection         ::LineSelection 
        maxPartialWave        ::Int64
    end 


    """
    `ParticleScattering.Settings()`  ... constructor for the default ParticleScattering.Settings.
    """
    function Settings()
        Settings(ElasticElectron(), PlaneWave(), LinearX(), Float64[], Float64[], Float64[], false, LineSelection(), 2)
    end


    """
    `ParticleScattering.Settings(set::ParticleScattering.Settings;`
    
            processType=..,         beamType=..,                polarization=..,          
            impactEnergies=..,      polarThetas=..,             polarPhis=.., 
            printBefore=..,         lineSelection=..,           maxPartialWave=.. )
                        
        ... constructor for modifying the given ParticleScattering.Settings by 'overwriting' the previously selected parameters.
    """
    function Settings(set::ParticleScattering.Settings; 
        processType::Union{Nothing,ParticleScattering.AbstractProcessType}=nothing,
        beamType::Union{Nothing,Beam.AbstractBeamType}=nothing,
        polarization::Union{Nothing,Basics.AbstractPolarization}=nothing,
        impactEnergies::Union{Nothing,Array{Float64,1}}=nothing,    polarThetas::Union{Nothing,Array{Float64,1}}=nothing,
        polarPhis::Union{Nothing,Array{Float64,1}}=nothing,         printBefore::Union{Nothing,Bool}=nothing, 
        lineSelection::Union{Nothing,LineSelection}=nothing,        maxPartialWave::Union{Nothing,Int64}=nothing)  
        
        if  processTypey   == nothing   processTypex    = set.processType       else  processTypex    = processType       end
        if  beamType       == nothing   beamTypex       = set.beamType          else  beamTypex       = beamType          end
        if  polarization   == nothing   polarizationx   = set.polarization      else  polarizationx   = polarization      end
        if  impactEnergies == nothing   impactEnergiesx = set.impactEnergies    else  impactEnergiesx = impactEnergies    end
        if  polarThetas    == nothing   polarThetasx    = set.polarThetas       else  polarThetasx    = polarThetas       end
        if  polarPhis      == nothing   polarPhisx      = set.polarPhis         else  polarPhisx      = polarPhis         end
        if  printBefore    == nothing   printBeforex    = set.printBefore       else  printBeforex    = printBefore       end 
        if  lineSelection  == nothing   lineSelectionx  = set.lineSelection     else  lineSelectionx  = lineSelection     end 
        if  maxPartialWave == nothing   maxPartialWavex = set.maxPartialWave    else  maxPartialWavex = maxPartialWave    end

        Settings( processTypex, beamTypex, polarizationx, impactEnergiesx, polarThetasx, polarPhisx, printBeforex, 
                  lineSelectionx, maxPartialWavex )
    end


    # `Base.show(io::IO, settings::ParticleScattering.Settings)`  
    #       ... prepares a proper printout of the variable settings::ParticleScattering.Settings.
    function Base.show(io::IO, settings::ParticleScattering.Settings) 
        println(io, "processType:           $(settings.processType)  ")
        println(io, "beamType:              $(settings.beamType)  ")
        println(io, "polarization:          $(settings.polarization)  ")
        println(io, "impactEnergies :       $(settings.impactEnergies )  ")
        println(io, "polarThetas:           $(settings.polarThetas)  ")
        println(io, "polarPhis:             $(settings.polarPhis)  ")
        println(io, "printBefore:           $(settings.printBefore)  ")
        println(io, "lineSelection:         $(settings.lineSelection)  ")
        println(io, "maxPartialWave:        $(settings.maxPartialWave)  ")
    end
        
    
    """
    `struct  ParticleScattering.PartialWaveNR`  
        ... defines a type to represent a single partial-wave amplitude.

        + l                   ::Int64      ... OAM of the partial wave.
        + theta               ::Float64    ... polar theta of partial-wave amplitude.
        + phi                 ::Float64    ... polar phi of partial-wave amplitude.
        + phase               ::Float64    ... phase shift of partial-wave.
        + amplitude           ::ComplexF64 ... partial-wave amplitude
    """
    struct PartialWaveNR
        l                     ::Int64   
        theta                 ::Float64
        phi                   ::Float64
        phase                 ::Float64  
        amplitude             ::ComplexF64 
    end 


    """
    `ParticleScattering.PartialWaveNR()`  ... constructor for the default ParticleScattering.PartialWaveNR.
    """
    function PartialWaveNR()
        PartialWaveNR(0, 0., 0., 0., 0.)
    end


    # `Base.show(io::IO, pw::ParticleScattering.PartialWaveNR)`  ... printout of the variable pw::ParticleScattering.PartialWaveNR.
    function Base.show(io::IO, pw::ParticleScattering.PartialWaveNR) 
        println(io, "l:               $(pw.l)  ")
        println(io, "theta:           $(pw.theta)  ")
        println(io, "phi:             $(pw.phi)  ")
        println(io, "phase:           $(pw.phase)  ")
        println(io, "amplitude:       $(pw.amplitude)  ")
    end
        
    
    """
    `struct  ParticleScattering.LineNR`  
        ... defines a type to collect data & results for a (no-relativistic) scattering "line".

        + processType    ::ParticleScattering.AbstractProcessType
        + beamType       ::Beam.AbstractBeamType
        + initialLevel   ::Level           ... initial-(state) level
        + finalLevel     ::Level           ... final-(state) level
        + impactEnergy   ::Float64         ... Energy of the (incoming) particle.
        + partialWaves   ::Array{ParticleScattering.PartialWaveNR,1}
    """
    struct LineNR
        processType      ::ParticleScattering.AbstractProcessType
        beamType         ::Beam.AbstractBeamType
        initialLevel     ::Level  
        finalLevel       ::Level 
        impactEnergy     ::Float64  
        partialWaves     ::Array{ParticleScattering.PartialWaveNR,1}
    end 


    """
    `ParticleScattering.LineNR()`  ... constructor for the default ParticleScattering.LineNR.
    """
    function LineNR()
        LineNR(ElasticElectronNR(), Beam.PlaneWave(), Level(), Level(), 0., PartialWaveNR[])
    end


    # `Base.show(io::IO, pw::ParticleScattering.LineNR)`  ... printout of the variable pw::ParticleScattering.LineNR.
    function Base.show(io::IO, line::ParticleScattering.LineNR) 
        println(io, "processType:        $(line.processType)  ")
        println(io, "beamType:           $(line.beamType)  ")
        println(io, "initialLevel:       $(line.initialLevel)  ")
        println(io, "finalLevel:         $(line.finalLevel)  ")
        println(io, "impactEnergy:       $(line.impactEnergy)  ")
        println(io, "partialWaves:       $(line.partialWaves)  ")
    end

    
    """
    `ParticleScattering.amplitude(processType::ElasticElectronNR, beamType::Beam.PlaneWave, 
                                  l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)`
        ... to compute the partial-wave (l) amplitude for the nonrelativistic elastic electron scattering of plane-wave beam 
            at given theta and phi. An  amplitude::ComplexF64 is returned.
    """
    function amplitude(processType::ElasticElectronNR, beamType::Beam.PlaneWave, 
                       l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)
        amplitude = ComplexF64(0.)
        if  beamType.kx == beamType.ky == 0.
            amplitude = amplitude + (2*l + 1) * exp( im*lPhase ) * sin(lPhase) * GSL.sf_legendre_Pl_e(l, cos(theta)).val
            amplitude = amplitude / beamType.kz
        else  error("stop a")
        end
        @show processType, beamType, l, amplitude
        
        return( amplitude )
    end

    
    """
    `ParticleScattering.amplitude(processType::ElasticElectronNR, beamType::Beam.BesselBeam, 
                                  l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)`
        ... to compute the partial-wave (l) amplitude for the nonrelativistic elastic electron scattering of Bessel beam 
            at given theta and phi. An  amplitude::ComplexF64 is returned.
    """
    function amplitude(processType::ElasticElectronNR, beamType::Beam.BesselBeam, 
                       l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)
        amplitude = ComplexF64(0.);   
        #
        k         = beamType.kz / cos(beamType.kz)
        wa        = 2 * (-im)^beamType.mOAM / k * exp( im*lPhase ) * (-1)^beamType.mOAM             *
                    sqrt( (2*l + 1) * factorial(l - beamType.mOAM) / factorial(l + beamType.mOAM) ) * 
                    GSL.sf_legendre_Plm(l, beamType.mOAM, cos(theta))                               * 
                    AngularMomentum.sphericalYlm(l, beamType.mOAM, theta, phi)
        amplitude = amplitude + wa
        @show processType, beamType, l, amplitude
        
        return( amplitude )
    end


    """
    `ParticleScattering.computeAmplitudesProperties(processType::ElasticElectronNR, line::ParticleScattering.LineNR, 
                                                    nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                    settings::ParticleScattering.Settings; printout::Bool=true)` 
        ... to compute all amplitudes and properties of the given line; a line::ParticleScattering.Line is returned for which the amplitudes 
            and properties are now evaluated.
    """
    function computeAmplitudesProperties(processType::ElasticElectronNR, line::ParticleScattering.LineNR, nm::Nuclear.Model, 
                                         grid::Radial.Grid, nrContinuum::Int64, settings::ParticleScattering.Settings; printout::Bool=true) 
        newPws = ParticleScattering.PartialWaveNR[];   contSettings = Continuum.Settings(false, nrContinuum)
        
        for pw in line.partialWaves
            #==
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, subshellList)
            newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
            cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
            newChannel = AutoIonization.Channel(channel.kappa, channel.symmetry, phase, 0.) ==#
            
            lPhase     = 0.11
            amplitude  = ParticleScattering.amplitude(line.processType, line.beamType, pw.l, lPhase, pw.theta, pw.phi, grid)
            push!( newPws, ParticleScattering.PartialWaveNR(pw.l, pw.theta, pw.phi, lPhase, amplitude) )
        end
        newLine   = ParticleScattering.LineNR(line.processType, line.beamType, line.initialLevel, line.finalLevel, line.impactEnergy, newPws)

        return( newLine )
    end


    """
    `ParticleScattering.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                 settings::ParticleScattering.Settings; output=true, printout::Bool=true)`  
        ... to compute the particle scattering amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{ParticleScattering.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::ParticleScattering.Settings; output=true, printout::Bool=true)
        println("")
        printstyled("ParticleScattering.computeLines(): The computation of Auger rates and properties starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        if !( settings.processType in [ElasticElectronNR()] )   error("stop a")    end
        #
        # Distinguish different sub-procedures ... if relativistic scattering will be considered in the future
        lines = ParticleScattering.determineLinesNR(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    ParticleScattering.displayLines(settings.processType, lines)    end  
        # Determine maximum energy and check for consistency of the grid
        maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.impactEnergy)   end
        nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newLines = ParticleScattering.LineNR[]
        for  line in lines
            newLine = ParticleScattering.computeAmplitudesProperties(settings.processType, line, nm, grid, nrContinuum, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        ParticleScattering.displayAmplitudes(stdout, settings.processType, newLines, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   ParticleScattering.displayAmplitudes(iostream, settings.processType, newLines, settings)    end
        ==#
        if    output    return( newLines )
        else            return( nothing )
        end
    end


    """
    `ParticleScattering.determinePartialWavesNR(finalLevel::Level, initialLevel::Level, settings::ParticleScattering.Settings)`  
        ... to determine a list of (non-relativistic) partial waves for a transitions from the initial to final level and 
            by taking into account the particular settings of for this computation; 
            an Array{ParticleScattering.PartialWaveNR,1} is returned.
    """
    function determinePartialWavesNR(finalLevel::Level, initialLevel::Level, settings::ParticleScattering.Settings)
        pws = ParticleScattering.PartialWaveNR[]
        for  l in 0:settings.maxPartialWave
            for  theta in settings.polarThetas
                for  phi in settings.polarPhis
                    push!(pws, ParticleScattering.PartialWaveNR(l, theta, phi, 0., 0.) )
                end
            end
        end
        return( pws )  
    end


    """
    `ParticleScattering.determineLinesNR(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ParticleScattering.Settings)`  
        ... to determine a list of ParticleScattering.Line's for transitions between levels from the initial- and final-state multiplets, and  
            by taking into account the particular selections and settings for this computation; an Array{ParticleScattering.Line,1} is returned. 
            Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLinesNR(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ParticleScattering.Settings)
        lines = ParticleScattering.LineNR[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    pws = ParticleScattering.determinePartialWavesNR(fLevel, iLevel, settings) 
                    for en in settings.impactEnergies
                        push!( lines, ParticleScattering.LineNR(settings.processType, settings.beamType, iLevel, fLevel, en, pws) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `ParticleScattering.displayLines(processType::ElasticElectronNR, lines::Array{ParticleScattering.LineNR,1})`  
        ... to display a list of lines and partial waves that have been selected due to the prior settings. 
            A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(processType::ElasticElectronNR, lines::Array{ParticleScattering.LineNR,1})
        nx = 106
        println(" ")
        println("  Selected elastic electron scattering lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(20, "Beam type"; na=2);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(14, "Impact energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(18, "No partial waves"; na=2);                         sb = sb * TableStrings.hBlank(20)           
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(20, string(typeof(line.beamType)); na=2 ) 
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.impactEnergy))         * "    "
            sa = sa * TableStrings.center(18, string( length(line.partialWaves)); na=2 )
            println( sa )
        end
        println("  ", TableStrings.hLine(nx), "\n")
        #
        return( nothing )
    end


    """
    `ParticleScattering.displayAmplitudes(stream::IO, processType::ElasticElectronNR, lines::Array{ParticleScattering.LineNR,1},
                                          settings::ParticleScattering.Settings)`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayAmplitudes(stream::IO, processType::ElasticElectronNR, lines::Array{ParticleScattering.LineNR,1},
                                settings::ParticleScattering.Settings)
        nx = 158
        println(stream, " ")
        println(stream, "  Selected elastic electron scattering lines:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(20, "Beam type"; na=2);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(14, "Impact energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=6)
        sa = sa * TableStrings.center(10, "theta"; na=2);                                    sb = sb * TableStrings.hBlank(12)               
        sa = sa * TableStrings.center(10, "phi";   na=4);                                    sb = sb * TableStrings.hBlank(14)               
        sa = sa * TableStrings.center( 3, "l";     na=6);                                    sb = sb * TableStrings.hBlank( 9)               
        sa = sa * TableStrings.center(32, "Partial-wave amplitudes";  na=2);                  
        sb = sb * TableStrings.center(32, "Re-Amplitude-Im";          na=2)           
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #
        nsa = 0
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(20, string(typeof(line.beamType)); na=2 ) 
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.impactEnergy))  * "      ";     nsa = length(sa)
            sa = sa * @sprintf("%.3e", line.partialWaves[1].theta) * "  "
            sa = sa * @sprintf("%.3e", line.partialWaves[1].phi)   * "     "
            sa = sa * string(line.partialWaves[1].l)               * "   " 
            sa = sa * string(line.partialWaves[1].amplitude.re)    * "   " 
            sa = sa * string(line.partialWaves[1].amplitude.im)    * "   " 
            println(stream,  sa)
            for  j in 2:length(line.partialWaves)
                sa = " "^nsa
                sa = sa * @sprintf("%.3e", line.partialWaves[1].theta) * "  "
                sa = sa * @sprintf("%.3e", line.partialWaves[1].phi)   * "     "
                sa = sa * string(line.partialWaves[1].l)               * "   " 
                sa = sa * string(line.partialWaves[1].amplitude.re)    * "   " 
                sa = sa * string(line.partialWaves[1].amplitude.im)    * "   " 
                println(stream,  sa)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx), "\n")
        #
        return( nothing )
    end
    

end # module
