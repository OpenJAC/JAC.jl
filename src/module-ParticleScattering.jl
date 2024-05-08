
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
        for selected scattering event.

    + processType         ::ParticleScattering.AbstractProcessType
    + beamType            ::Beam.AbstractBeamType
    + polarization        ::Basics.AbstractPolarization
    + impactEnergies      ::Array{Float64,1}
    + polarThetas         ::Array{Float64,1}
    + polarPhis           ::Array{Float64,1}
    + printBefore         ::Bool               ... True, if all energies and events are printed before their evaluation.
    + lineSelection       ::LineSelection      ... Specifies the selected levels, if any.
    + epsPartialWave      ::Float64            
        ... accuracy criterion; all partial waves (l, kappa) with 
            abs(contribution_l) + abs(contribution_l+1) < epsPartialWave are neglected;
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
    epsPartialWave        ::Float64
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
        printBefore=..,         lineSelection=..,           epsPartialWave=.. )
                    
    ... constructor for modifying the given ParticleScattering.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::ParticleScattering.Settings; 
    processType::Union{Nothing,ParticleScattering.AbstractProcessType}=nothing,
    beamType::Union{Nothing,Beam.AbstractBeamType}=nothing,
    polarization::Union{Nothing,Basics.AbstractPolarization}=nothing,
    impactEnergies::Union{Nothing,Array{Float64,1}}=nothing,    polarThetas::Union{Nothing,Array{Float64,1}}=nothing,
    polarPhis::Union{Nothing,Array{Float64,1}}=nothing,         printBefore::Union{Nothing,Bool}=nothing, 
    lineSelection::Union{Nothing,LineSelection}=nothing,        epsPartialWave::Union{Nothing,Float64}=nothing)  
    
    if  processTypey   == nothing   processTypex    = set.processType       else  processTypex    = processType       end
    if  beamType       == nothing   beamTypex       = set.beamType          else  beamTypex       = beamType          end
    if  polarization   == nothing   polarizationx   = set.polarization      else  polarizationx   = polarization      end
    if  impactEnergies == nothing   impactEnergiesx = set.impactEnergies    else  impactEnergiesx = impactEnergies    end
    if  polarThetas    == nothing   polarThetasx    = set.polarThetas       else  polarThetasx    = polarThetas       end
    if  polarPhis      == nothing   polarPhisx      = set.polarPhis         else  polarPhisx      = polarPhis         end
    if  printBefore    == nothing   printBeforex    = set.printBefore       else  printBeforex    = printBefore       end 
    if  lineSelection  == nothing   lineSelectionx  = set.lineSelection     else  lineSelectionx  = lineSelection     end 
    if  epsPartialWave == nothing   epsPartialWavex = set.epsPartialWave    else  epsPartialWavex = epsPartialWave    end

    Settings( processTypex, beamTypex, polarizationx, impactEnergiesx, polarThetasx, polarPhisx, printBeforex, 
                lineSelectionx, epsPartialWavex )
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
    println(io, "epsPartialWave:        $(settings.epsPartialWave)  ")
end
    

"""
`struct  ParticleScattering.PartialWaveNR`  
    ... defines a type to represent a single partial-wave amplitude.

    + l                   ::Int64      ... OAM of the partial wave.
    + phase               ::Float64    ... phase shift of partial-wave.
    + amplitude           ::ComplexF64 ... partial-wave amplitude
"""
struct PartialWaveNR
    l                     ::Int64   
    phase                 ::Float64  
    amplitude             ::ComplexF64 
end 


"""
`ParticleScattering.PartialWaveNR()`  ... constructor for the default ParticleScattering.PartialWaveNR.
"""
function PartialWaveNR()
    PartialWaveNR(0, 0., 0.)
end


# `Base.show(io::IO, pw::ParticleScattering.PartialWaveNR)`  ... printout of the variable pw::ParticleScattering.PartialWaveNR.
function Base.show(io::IO, pw::ParticleScattering.PartialWaveNR) 
    println(io, "l:               $(pw.l)  ")
    println(io, "phase:           $(pw.phase)  ")
    println(io, "amplitude:       $(pw.amplitude)  ")
end
    

"""
`struct  ParticleScattering.EventNR`  
    ... defines a type to collect data & results for a (no-relativistic) scattering "event".

    + processType    ::ParticleScattering.AbstractProcessType
    + beamType       ::Beam.AbstractBeamType
    + initialLevel   ::Level           ... initial-(state) level
    + finalLevel     ::Level           ... final-(state) level
    + impactEnergy   ::Float64         ... Energy of the (incoming) particle.
    + theta          ::Float64         ... polar theta of partial-wave amplitude.
    + phi            ::Float64         ... polar phi of partial-wave amplitude.
    + d2Sigma        ::Float64         ... Triple-differential CS: d^2 sigma / d Omega
    + dSigmadE       ::Float64         ... Energy-differential CS: d sigma / d epsilon
    + partialWaves   ::Array{ParticleScattering.PartialWaveNR,1}
"""
struct EventNR
    processType      ::ParticleScattering.AbstractProcessType
    beamType         ::Beam.AbstractBeamType
    initialLevel     ::Level  
    finalLevel       ::Level 
    impactEnergy     ::Float64  
    theta            ::Float64
    phi              ::Float64
    d2Sigma          ::Float64 
    dSigmadE         ::Float64
    partialWaves     ::Array{ParticleScattering.PartialWaveNR,1}
end 


"""
`ParticleScattering.EventNR()`  ... constructor for the default ParticleScattering.EventNR.
"""
function EventNR()
    EventNR(ElasticElectronNR(), Beam.PlaneWave(), Level(), Level(), 0., 0., 0., 0., 0., PartialWaveNR[])
end


# `Base.show(io::IO, event::ParticleScattering.EventNR)`  ... printout of the variable event::ParticleScattering.EventNR.
function Base.show(io::IO, event::ParticleScattering.EventNR) 
    println(io, "processType:        $(event.processType)  ")
    println(io, "beamType:           $(event.beamType)  ")
    println(io, "initialLevel:       $(event.initialLevel)  ")
    println(io, "finalLevel:         $(event.finalLevel)  ")
    println(io, "impactEnergy:       $(event.impactEnergy)  ")
    println(io, "theta:              $(event.theta)  ")
    println(io, "phi:                $(event.phi)  ")
    println(io, "d2Sigma:            $(event.d2Sigma)  ")
    println(io, "dSigmadE:           $(event.dSigmadE)  ")
    println(io, "partialWaves:       $(event.partialWaves)  ")
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
        amplitude = amplitude + (2*l + 1) / (4*pi) * exp( im*lPhase ) * sin(lPhase) * GSL.sf_legendre_Pl_e(l, cos(theta)).val
        amplitude = 4pi * amplitude / beamType.kz
    else  error("stop a")
    end
    ##x @show processType, beamType, l, amplitude
    
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
    amplitude = ComplexF64(0.)
    # 
    if  l >= abs(beamType.mOAM)
        # Only partial amplitudes with l >= abs(mOAM) are nonzero
        k         = beamType.kz / cos(beamType.openingAngle)
        wa        = GSL.sf_legendre_Plm(l, abs(beamType.mOAM), cos(beamType.openingAngle));   wc = 1.0
        if   beamType.mOAM < 0
            m  = abs(beamType.mOAM)
            wc = (-1)^m * factorial(big(l-m)) / factorial(big(l+m))     
            wc = Float64(wc)
        end
        wb        = 4pi * (-1.0im)^beamType.mOAM / k * exp( im*lPhase ) * (-1.)^beamType.mOAM                  *
                    sqrt( (2*l + 1) * factorial( big(l - beamType.mOAM) ) / (4*pi * factorial( big(l + beamType.mOAM) ) ) ) * 
                    wa * wc * AngularMomentum.sphericalYlm(l, beamType.mOAM, theta, phi) * sin(lPhase)
        amplitude = amplitude + ComplexF64(wb)
    end
    @show processType, beamType, l, amplitude, typeof(amplitude)
    
    return( amplitude )
end


"""
`ParticleScattering.computeAmplitudesProperties(processType::ElasticElectronNR, event::ParticleScattering.EventNR, 
                                                nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                settings::ParticleScattering.Settings; printout::Bool=true)` 
    ... to compute all amplitudes and properties of the given event; a event::ParticleScattering.EventNR is returned 
        for which the amplitudes and properties are now evaluated.
"""
function computeAmplitudesProperties(processType::ElasticElectronNR, event::ParticleScattering.EventNR, nm::Nuclear.Model, 
                                        grid::Radial.Grid, nrContinuum::Int64, settings::ParticleScattering.Settings; printout::Bool=true) 
    newPws   = ParticleScattering.PartialWaveNR[];   contSettings = Continuum.Settings(false, nrContinuum)
    totalAmp = ComplexF64(0.);   amp1 = amp2 = amp3 = amp4 = 1.0e6;   maxamp = 0.
    
    # For elastic scattering, the initial and final levels are expected to be the same
    for  l = 0:1000
        @show  l, amp1, amp2, amp3, amp4, maxamp, settings.epsPartialWave
        kappa      = -l - 1;    
        ##x symf       = LevelSymmetry(event.finalLevel.J, event.finalLevel.parity)
        ##x symmetries = AngularMomentum.allowedTotalSymmetries(symf, kappa::Int64)
        ##x @show symmetries
        cOrbital, lPhase  = Continuum.generateOrbitalForLevel(event.impactEnergy, Subshell(101, kappa), event.finalLevel, 
                                                                nm, grid, contSettings)
        ##x newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, symmetries[1], event.finalLevel)
        amplitude  = ParticleScattering.amplitude(event.processType, event.beamType, l, lPhase, event.theta, event.phi, grid)
        totalAmp   = totalAmp + amplitude
        push!( newPws, ParticleScattering.PartialWaveNR(l, lPhase, amplitude) )
        amp1 = amp2;   amp2 = amp3;   amp3 = amp4;   amp4 = abs(amplitude)^2;   maxamp = max(maxamp, amp4)
        if  l > 10  &&  (amp1 + amp2 + amp3 + amp4) / maxamp < settings.epsPartialWave    break     end
    end
    @show typeof(totalAmp)
    d2Sigma  = conj(totalAmp) * totalAmp
    newEvent = ParticleScattering.EventNR(event.processType, event.beamType, event.initialLevel, event.finalLevel, event.impactEnergy, 
                                            event.theta, event.phi, d2Sigma, 0., newPws)

    return( newEvent )
end


"""
`ParticleScattering.computeEvents(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::ParticleScattering.Settings; output=true, printout::Bool=true)`  
    ... to compute the particle scattering amplitudes and all properties as requested by the given settings. A list of 
        events::Array{ParticleScattering.Events} is returned.
"""
function  computeEvents(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                        settings::ParticleScattering.Settings; output=true, printout::Bool=true)
    println("") 
    printstyled("ParticleScattering.computeEvents(): The computation of Auger rates and properties starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------ \n", color=:light_green)
    println("")
    #
    if !( settings.processType in [ElasticElectronNR()] )   error("stop a")    end
    #
    # Distinguish different sub-procedures ... if relativistic scattering will be considered in the future
    events = ParticleScattering.determineEventsNR(finalMultiplet, initialMultiplet, settings)
    # Display all selected events before the computations start
    if  settings.printBefore    ParticleScattering.displayEvents(settings.processType, events)    end  
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  event in events   maxEnergy = max(maxEnergy, event.impactEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newEvents = ParticleScattering.EventNR[]
    for  event in events
        newEvent = ParticleScattering.computeAmplitudesProperties(settings.processType, event, nm, grid, nrContinuum, settings) 
        push!( newEvents, newEvent)
    end
    # Print all results to screen
    ParticleScattering.displayAmplitudes(stdout, settings.processType, newEvents, settings)
    ParticleScattering.displayCrossSections(stdout, settings.processType, newEvents, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   ParticleScattering.displayAmplitudes(iostream, settings.processType, newEvents, settings)
                        ParticleScattering.displayCrossSections(iostream, settings.processType, newEvents, settings)    end

    if    output    return( newEvents )
    else            return( nothing )
    end
end


#==
"""
`ParticleScattering.determinePartialWavesNR(finalLevel::Level, initialLevel::Level, settings::ParticleScattering.Settings)`  
    ... to determine a list of (non-relativistic) partial waves for a transitions from the initial to final level and 
        by taking into account the particular settings of for this computation; 
        an Array{ParticleScattering.PartialWaveNR,1} is returned.
"""
function determinePartialWavesNR(finalLevel::Level, initialLevel::Level, settings::ParticleScattering.Settings)
    pws = ParticleScattering.PartialWaveNR[]
    println(">>> Partial waves are calculated from l = 0:settings.maxPartialWave independent of the incoming " * 
            "beamType = $(settings.beamType)")
    #
    for  l in 0:settings.maxPartialWave
        push!(pws, ParticleScattering.PartialWaveNR(l, 0., 0.) )
    end
    return( pws )  
end  ==#


"""
`ParticleScattering.determineEventsNR(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ParticleScattering.Settings)`  
    ... to determine a list of ParticleScattering.Event's for transitions between levels from the initial- and final-state multiplets, and  
        by taking into account the particular selections and settings for this computation; an Array{ParticleScattering.Event,1} is returned. 
        Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineEventsNR(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ParticleScattering.Settings)
    events = ParticleScattering.EventNR[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                ##x pws = ParticleScattering.determinePartialWavesNR(fLevel, iLevel, settings) 
                for en in settings.impactEnergies
                    enau    = Defaults.convertUnits("energy: to atomic", en)
                    newBeam = Beam.redefineEnergy(enau, settings.beamType)
                    for  theta in settings.polarThetas,   phi in settings.polarPhis
                        push!(events, ParticleScattering.EventNR(settings.processType, newBeam, iLevel, fLevel, enau, 
                                                            theta, phi, 0., 0., [ParticleScattering.PartialWaveNR(0, 0., 0.)]) )
                    end
                end
            end
        end
    end
    return( events )
end


"""
`ParticleScattering.displayEvents(processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1})`  
    ... to display a list of events and partial waves that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayEvents(processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1})
    nx = 106
    println(" ")
    println("  Selected elastic electron scattering events:")
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
    for  event in events
        sa  = "  ";    isym = LevelSymmetry( event.initialLevel.J, event.initialLevel.parity)
                        fsym = LevelSymmetry( event.finalLevel.J,   event.finalLevel.parity)
        sa = sa * TableStrings.center(20, string(typeof(event.beamType)); na=2 ) 
        sa = sa * TableStrings.center(18, TableStrings.levels_if(event.initialLevel.index, event.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", event.impactEnergy))         * "    "
        sa = sa * TableStrings.center(18, string( length(event.partialWaves)); na=2 )
        println( sa )
    end
    println("  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end


"""
`ParticleScattering.displayAmplitudes(stream::IO, processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1},
                                        settings::ParticleScattering.Settings)`  
    ... to display a list of events and partial waves that have been selected due to the prior settings. A neat table of all selected 
        transitions, energies, angles and partial-wave amplitudes is printed but nothing is returned otherwise.
"""
function  displayAmplitudes(stream::IO, processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1},
                            settings::ParticleScattering.Settings)
    nx = 158
    println(stream, " ")
    println(stream, "  Selected elastic electron scattering events:  Scattering amplitudes")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(20, "Beam type"; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Impact energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=6)
    sa = sa * TableStrings.center(10, "theta"; na=2);                                    sb = sb * TableStrings.hBlank(12)               
    sa = sa * TableStrings.center(10, "phi";   na=3);                                    sb = sb * TableStrings.hBlank(13)               
    sa = sa * TableStrings.center( 3, "l";     na=2);                                    sb = sb * TableStrings.hBlank( 3)               
    sa = sa * TableStrings.center(30, "Partial-wave amplitudes";  na=2);                  
    sb = sb * TableStrings.center(30, "Re-Amplitude-Im";          na=2)           
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    nsa = 0
    for  event in events
        sa  = "  ";    isym = LevelSymmetry( event.initialLevel.J, event.initialLevel.parity)
                        fsym = LevelSymmetry( event.finalLevel.J,   event.finalLevel.parity)
        sa = sa * TableStrings.center(18, string(typeof(event.beamType)); na=1 ) 
        sa = sa * TableStrings.center(18, TableStrings.levels_if(event.initialLevel.index, event.finalLevel.index); na=1)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", event.impactEnergy))  * "     ";     nsa = length(sa)
        sa = sa * @sprintf("%.3e", event.theta)                           * "  "
        sa = sa * @sprintf("%.3e", event.phi)                             * "     "
        sa = sa * string(event.partialWaves[1].l)                         * "   " 
        sa = sa * @sprintf("% 1.6e", event.partialWaves[1].amplitude.re)    * "   " 
        sa = sa * @sprintf("% 1.6e", event.partialWaves[1].amplitude.im)    * "   " 
        println(stream,  sa)
        for  j in 2:length(event.partialWaves)
            sa = " "^nsa
            sa = sa * @sprintf("%.3e", event.theta)                           * "  "
            sa = sa * @sprintf("%.3e", event.phi)                             * "     "
            sa = sa * string(event.partialWaves[j].l)                         * "   " 
            sa = sa * @sprintf("% 1.6e", event.partialWaves[j].amplitude.re)    * "   " 
            sa = sa * @sprintf("% 1.6e", event.partialWaves[j].amplitude.im)    * "   " 
            println(stream,  sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end


"""
`ParticleScattering.displayCrossSections(stream::IO, processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1},
                                            settings::ParticleScattering.Settings)`  
    ... to display a list of events and associated triple-differential and energy-differential cross sections. 
        A neat table of all selected transitions, energies and cross sections is printed but nothing is returned otherwise.
"""
function  displayCrossSections(stream::IO, processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1},
                                settings::ParticleScattering.Settings)
    nx = 158
    println(stream, " ")
    println(stream, "  Selected elastic electron scattering events:  Differential cross sections for $(settings.beamType)")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(20, "Beam type"; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Impact energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=6)
    sa = sa * TableStrings.center(10, "theta"; na=3);                                    sb = sb * TableStrings.hBlank(12)               
    sa = sa * TableStrings.center(10, "phi";   na=4);                                    sb = sb * TableStrings.hBlank(12)               
    sa = sa * TableStrings.center(22, "d^2 sigma /d angles";  na=2); 
    sc = TableStrings.inUnits("cross section") 
    sb = sb * TableStrings.center(22, sc * "/sr" ;  na=2)           
    sa = sa * TableStrings.center(14, "d sigma / dE";  na=2);                  
    sb = sb * TableStrings.center(14, sc * "/dE" ;  na=2)           
    sa = sa * TableStrings.center(10, "kz";   na=4);                                     sb = sb * TableStrings.hBlank(14)               
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    nsa = 0;    d2Convert = Defaults.convertUnits("cross section: from atomic to predefined unit", 1.0) 
                deConvert = d2Convert / Defaults.convertUnits("energy: from atomic to predefined unit", 1.0)
    for  event in events
        sa  = "  ";    isym = LevelSymmetry( event.initialLevel.J, event.initialLevel.parity)
                        fsym = LevelSymmetry( event.finalLevel.J,   event.finalLevel.parity)
        sa = sa * TableStrings.center(18, string(typeof(event.beamType)); na=1 ) 
        sa = sa * TableStrings.center(18, TableStrings.levels_if(event.initialLevel.index, event.finalLevel.index); na=1)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", event.impactEnergy))  * "     ";     nsa = length(sa)
        sa = sa * @sprintf("%.3e", event.theta)                   * "  "
        sa = sa * @sprintf("%.3e", event.phi)                     * "           "
        sa = sa * @sprintf("%.6e", event.d2Sigma * d2Convert)     * "        " 
        sa = sa * @sprintf("%.6e", event.dSigmadE * deConvert)    * "   " 
        sa = sa * @sprintf("%.3e", event.beamType.kz)             * "   " 
        println(stream,  sa)
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end

"""
`ParticleScattering.extractCrossSections(events::Array{ParticleScattering.EventNR,1}, energy, phi)`  
    ... to extract all triple and energy-differential scattering cross sections for (constant) energy and phi which are 
        provided by the events. A name tuple (thetas::Array{Float64,1}, d2Sigmas::Array{Float64,1}, , dSigmadEs::Array{Float64,1})
        is returned which provide the corresponding values.
"""
function  extractCrossSections(events::Array{ParticleScattering.EventNR,1}, energy, phi)
    thetas = Float64[];    d2Sigmas = Float64[];    dSigmadEs = Float64[];  
    enau   = Defaults.convertUnits("energy: to atomic", energy)
    for event in events
        if   enau == event.impactEnergy   &&   phi == event.phi 
            push!(thetas, event.theta);   push!(d2Sigmas, event.d2Sigma);   push!(dSigmadEs, event.dSigmadE)
        end
    end
    
    wa = (thetas=thetas, d2Sigmas=d2Sigmas, dSigmadEs=dSigmadEs)
    return(wa)
end    

end # module
