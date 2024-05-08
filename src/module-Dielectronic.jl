
"""
`module  JAC.Dielectronic`  
... a submodel of JAC that contains all methods for computing dielectronic recombination properties between 
    some initial, intermediate and final-state multiplets.
"""
module Dielectronic


using Printf, SpecialFunctions,
        ..AngularMomentum, ..AutoIonization, ..Basics, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear, 
        ..PhotoEmission, ..Radial, ..TableStrings


"""
`abstract type Dielectronic.AbstractCorrections` 
    ... defines an abstract type to distinguish different types of corrections to the decay rates and strength; see also:
    
    + struct Dielectronic.NoCorrections  
        ... don't apply any additional correction to the resonance strengths.
    + struct Dielectronic.HydrogenicCorrections  
        ... to add for missing final decay levels the photon decay rates for non-relativistic hydrogenic ions;
            this improves the total photon rate as well as the resonance strength.

"""
abstract type  AbstractCorrections       end


"""
`struct  Dielectronic.NoCorrections          <:  Dielectronic.AbstractCorrections`  
    ... don't apply any additional correction to the resonance strengths.
"""
struct   NoCorrections                       <:  Dielectronic.AbstractCorrections   end


"""
`struct  Dielectronic.HydrogenicCorrections  <:  Dielectronic.AbstractCorrections`  
    ... to add for missing final decay levels the photon decay rates for non-relativistic hydrogenic ions;
        this improves the total photon rate as well as the resonance strength.

    + nDetailed       ::Int64   
        ... principal quantum number nDetailed up to which the photon rates are calculated explicitly;
            for all n-shells with nDetailed < n <= n_ryd-1, hydrogenic photon rates are scaled to account
            for the radiative decay into these high-n shells.
    + effectiveZ      ::Float64      ... effective charge Z_eff for the hydrogenic correction.
    + energyScaling   ::Float64      ... scaling factor to modify the non-relativistic energies.
"""
struct   HydrogenicCorrections               <:  Dielectronic.AbstractCorrections
    nDetailed         ::Int64   
    effectiveZ        ::Float64
    energyScaling     ::Float64
end


# `Base.show(io::IO, corr::HydrogenicCorrections)`  ... prepares a proper printout of the corr::HydrogenicCorrections.
function Base.show(io::IO, corr::HydrogenicCorrections)
    println(io, "nDetailed:       $(corr.nDetailed)  ")
    println(io, "effectiveZ:      $(corr.effectiveZ)  ")
    println(io, "energyScaling:   $(corr.energyScaling)  ")
end     


"""
`struct  Dielectronic.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing dielectronic recombination pathways.

    + multipoles            ::Array{EmMultipoles}  ... Multipoles of the radiation field that are to be included.
    + gauges                ::Array{UseGauge}      ... Specifies the gauges to be included into the computations.
    + calcRateAlpha         ::Bool                 ... True, if the DR rate coefficients are to be calculated, 
                                                        and false o/w.
    + printBefore           ::Bool                 ... True, if all energies and pathways are printed before their 
                                                        evaluation.
    + pathwaySelection      ::PathwaySelection     ... Specifies the selected levels/pathways, if any.
    + electronEnergyShift   ::Float64              ... An overall energy shift for all electron energies (i.e. from 
                                                        the initial to the resonance levels [Hartree].
    + photonEnergyShift     ::Float64              ... An overall energy shift for all photon energies (i.e. from 
                                                        the resonance to the final levels.
    + mimimumPhotonEnergy   ::Float64              ... minimum transition energy for which photon transitions are 
                                                        included into the evaluation.
    + temperatures          ::Array{Float64,1}     
        ... list of temperatures for which plasma rate coefficients are displayed; however, these rate coefficients
            only include the contributions from those pathsways that are calculated here explicitly.
    + corrections           ::Dielectronic.AbstractCorrections
        ... Specify, if appropriate, the inclusion of additional corrections to the rates and DR strengths.
    + augerOperator         ::AbstractEeInteraction 
        ... Auger operator that is to be used for evaluating the Auger amplitude's; the allowed values are: 
            CoulombInteraction(), BreitInteration(), CoulombBreit(), CoulombGaunt().
"""
struct Settings  <:  AbstractProcessSettings 
    multipoles              ::Array{EmMultipole,1}
    gauges                  ::Array{UseGauge}
    calcRateAlpha           ::Bool
    printBefore             ::Bool 
    pathwaySelection        ::PathwaySelection
    electronEnergyShift     ::Float64
    photonEnergyShift       ::Float64
    mimimumPhotonEnergy     ::Float64
    temperatures            ::Array{Float64,1}
    corrections             ::Dielectronic.AbstractCorrections
    augerOperator           ::AbstractEeInteraction
end 


"""
`Dielectronic.Settings()`  
    ... constructor for the default values of dielectronic recombination pathway computations.
"""
function Settings()
    Settings([E1], UseGauge[], false, false, PathwaySelection(), 0., 0., 0., Float64[], Dielectronic.NoCorrections(),
                CoulombInteraction())
end


"""
` (set::Dielectronic.Settings;`

        multipoles=..,           gauges=..,                  calcRateAlpha=..,           printBefore=..,
        pathwaySelection=..,     electronEnergyShift=..,     photonEnergyShift=..,       mimimumPhotonEnergy=..,     
        temperatures=..,         corrections=..,             augerOperator=..)
                    
    ... constructor for modifying the given Dielectronic.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::Dielectronic.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                        gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
    calcRateAlpha::Union{Nothing,Bool}=nothing,                                     printBefore::Union{Nothing,Bool}=nothing, 
    pathwaySelection::Union{Nothing,PathwaySelection}=nothing,                      electronEnergyShift::Union{Nothing,Float64}=nothing,
    photonEnergyShift::Union{Nothing,Float64}=nothing,                              mimimumPhotonEnergy::Union{Nothing,Float64}=nothing,
    temperatures::Union{Nothing,Array{Float64,1}}=nothing,                          corrections::Union{Nothing,AbstractCorrections}=nothing,
    augerOperator::Union{Nothing,AbstractEeInteraction}=nothing)
    
    if  multipoles          == nothing   multipolesx          = set.multipoles            else  multipolesx          = multipoles           end 
    if  gauges              == nothing   gaugesx              = set.gauges                else  gaugesx              = gauges               end 
    if  calcRateAlpha       == nothing   calcRateAlphax       = set.calcRateAlpha         else  calcRateAlphax       = calcRateAlpha        end 
    if  printBefore         == nothing   printBeforex         = set.printBefore           else  printBeforex         = printBefore          end 
    if  pathwaySelection    == nothing   pathwaySelectionx    = set.pathwaySelection      else  pathwaySelectionx    = pathwaySelection     end 
    if  electronEnergyShift == nothing   electronEnergyShiftx = set.electronEnergyShift   else  electronEnergyShiftx = electronEnergyShift  end 
    if  photonEnergyShift   == nothing   photonEnergyShiftx   = set.photonEnergyShift     else  photonEnergyShiftx   = photonEnergyShift    end 
    if  mimimumPhotonEnergy == nothing   mimimumPhotonEnergyx = set.mimimumPhotonEnergy   else  mimimumPhotonEnergyx = mimimumPhotonEnergy  end 
    if  temperatures        == nothing   temperaturesx        = set.temperatures          else  temperaturesx        = temperatures         end 
    if  corrections         == nothing   correctionsx         = set.corrections           else  correctionsx         = corrections          end 
    if  augerOperator       == nothing   augerOperatorx       = set.augerOperator         else  augerOperatorx       = augerOperator        end 

    Settings( multipolesx, gaugesx, calcRateAlphax, printBeforex, pathwaySelectionx, electronEnergyShiftx, 
                photonEnergyShiftx, mimimumPhotonEnergyx, temperaturesx, correctionsx, augerOperatorx )
end


# `Base.show(io::IO, settings::Dielectronic.Settings)`  ... prepares a proper printout of the variable settings::Dielectronic.Settings.
function Base.show(io::IO, settings::Dielectronic.Settings) 
    println(io, "multipoles:                 $(settings.multipoles)  ")
    println(io, "use-gauges:                 $(settings.gauges)  ")
    println(io, "calcRateAlpha:              $(settings.calcRateAlpha)  ")
    println(io, "printBefore:                $(settings.printBefore)  ")
    println(io, "pathwaySelection:           $(settings.pathwaySelection)  ")
    println(io, "electronEnergyShift:        $(settings.electronEnergyShift)  ")
    println(io, "photonEnergyShift:          $(settings.photonEnergyShift)  ")
    println(io, "mimimumPhotonEnergy:        $(settings.mimimumPhotonEnergy)  ")
    println(io, "temperatures:               $(settings.temperatures)  ")
    println(io, "corrections:                $(settings.corrections)  ")
    println(io, "augerOperator:              $(settings.augerOperator)  ")
end


"""
`struct  Dielectronic.Pathway`  
    ... defines a type for a dielectronic recombination pathways that may include the definition of channels and 
        their corresponding amplitudes.

    + initialLevel      ::Level                   ... initial-(state) level
    + intermediateLevel ::Level                   ... intermediate-(state) level
    + finalLevel        ::Level                   ... final-(state) level
    + electronEnergy    ::Float64                 ... energy of the (incoming, captured) electron
    + photonEnergy      ::Float64                 ... energy of the (emitted) photon
    + captureRate       ::Float64                 ... rate for the electron capture (Auger rate)
    + photonRate        ::EmProperty              ... rate for the photon emission
    + angularBeta       ::EmProperty              ... beta parameter of the photon emission
    + reducedStrength   ::EmProperty              ... reduced resonance strength S(i -> d -> f) * Gamma_d of this pathway;
                                                        this reduced strength does not require the knowledge of Gamma_d for each pathway.
    + captureChannels   ::Array{AutoIonization.Channel,1}   ... List of |i> -->  |n>   dielectronic (Auger) capture channels.
    + photonChannels    ::Array{PhotoEmission.Channel,1}    ... List of |n> -->  |f>   radiative stabilization channels.
"""
struct  Pathway
    initialLevel        ::Level
    intermediateLevel   ::Level
    finalLevel          ::Level
    electronEnergy      ::Float64
    photonEnergy        ::Float64 
    captureRate         ::Float64
    photonRate          ::EmProperty
    angularBeta         ::EmProperty
    reducedStrength     ::EmProperty
    captureChannels     ::Array{AutoIonization.Channel,1} 
    photonChannels      ::Array{PhotoEmission.Channel,1} 
end 


"""
`Dielectronic.Pathway()`  
    ... constructor for an 'empty' instance of a dielectronic recombination pathway between a specified 
        initial, intermediate and final level.
"""
function Pathway()
    em = EmProperty(0., 0.)
    Pathway(initialLevel, intermediateLevel, finalLevel, 0., 0., 0., em, em, em, AutoIonization.Channel[], PhotoEmission.Channel[])
end


# `Base.show(io::IO, pathway::Dielectronic.Pathway)`  ... prepares a proper printout of the variable pathway::Dielectronic.Pathway.
function Base.show(io::IO, pathway::Dielectronic.Pathway) 
    println(io, "initialLevel:               $(pathway.initialLevel)  ")
    println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
    println(io, "finalLevel:                 $(pathway.finalLevel)  ")
    println(io, "electronEnergy:             $(pathway.electronEnergy)  ")
    println(io, "photonEnergy:               $(pathway.photonEnergy)  ")
    println(io, "captureRate:                $(pathway.captureRate)  ")
    println(io, "photonRate:                 $(pathway.photonRate)  ")
    println(io, "angularBeta:                $(pathway.angularBeta)  ")
    println(io, "reducedStrength:            $(pathway.reducedStrength)  ")
    println(io, "captureChannels:            $(pathway.captureChannels)  ")
    println(io, "photonChannels:             $(pathway.photonChannels)  ")
end


"""
`struct  Dielectronic.Resonance`  
    ... defines a type for a dielectronic resonance as defined by a given initial and resonance level but by summing over all final levels

    + initialLevel      ::Level             ... initial-(state) level
    + intermediateLevel ::Level             ... intermediate-(state) level
    + resonanceEnergy   ::Float64           ... energy of the resonance w.r.t. the inital-state
    + resonanceStrength ::EmProperty        ... strength of this resonance due to the stabilization into any of the allowed final levels.
    + captureRate       ::Float64           ... capture (Auger) rate to form the intermediate resonance, starting from the initial level.
    + augerRate         ::Float64           ... total (Auger) rate for an electron emission of the intermediate resonance
    + photonRate        ::EmProperty        ... total photon rate for a photon emission, i.e. for stabilization.
"""
struct  Resonance
    initialLevel        ::Level
    intermediateLevel   ::Level
    resonanceEnergy     ::Float64 
    resonanceStrength   ::EmProperty
    captureRate         ::Float64
    augerRate           ::Float64
    photonRate          ::EmProperty
end 


"""
`Dielectronic.Resonance()`  
    ... constructor for an 'empty' instance of a dielectronic resonance as defined by a given initial and resonance 
        level but by summing over all final levels.
"""
function Resonance()
    em = EmProperty(0., 0.)
    Resonance(initialLevel, intermediateLevel, 0., em, 0., 0., em)
end


# `Base.show(io::IO, resonance::Dielectronic.Resonance)`  ... prepares a proper printout of the variable resonance::Dielectronic.Resonance.
function Base.show(io::IO, resonance::Dielectronic.Resonance) 
    println(io, "initialLevel:               $(resonance.initialLevel)  ")
    println(io, "intermediateLevel:          $(resonance.intermediateLevel)  ")
    println(io, "resonanceEnergy:            $(resonance.resonanceEnergy)  ")
    println(io, "resonanceStrength:          $(resonance.resonanceStrength)  ")
    println(io, "captureRate:                $(resonance.captureRate)  ")
    println(io, "augerRate:                  $(resonance.augerRate)  ")
    println(io, "photonRate:                 $(resonance.photonRate)  ")
end


"""
`struct  Dielectronic.ResonanceSelection`  
    ... defines a type for selecting classes of resonances in terms of leading configurations.

    + active          ::Bool              ... initial-(state) level
    + fromShells      ::Array{Shell,1}    ... List of shells from which excitations are to be considered.
    + toShells        ::Array{Shell,1}    ... List of shells to which (core-shell) excitations are to be considered.
    + intoShells      ::Array{Shell,1}    ... List of shells into which electrons are initially placed (captured).
"""
struct  ResonanceSelection
    active            ::Bool  
    fromShells        ::Array{Shell,1} 
    toShells          ::Array{Shell,1} 
    intoShells        ::Array{Shell,1} 
end 


"""
`Dielectronic.ResonanceSelection()`  
    ... constructor for an 'empty' instance of a ResonanceSelection()
"""
function ResonanceSelection()
    ResonanceSelection(false, Shell[], Shell[], Shell[] )
end


# `Base.show(io::IO, resonance::Dielectronic.ResonanceSelection)`  ... prepares a proper printout of resonance::Dielectronic.ResonanceSelection.
function Base.show(io::IO, rSelection::Dielectronic.ResonanceSelection) 
    println(io, "active:           $(rSelection.active)  ")
    println(io, "fromShells:       $(rSelection.fromShells)  ")
    println(io, "toShells:         $(rSelection.toShells)  ")
    println(io, "intoShells:       $(rSelection.intoShells)  ")
end





"""
`Dielectronic.checkConsistentMultiplets(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet)`  
    ... to check that the given initial-, intermediate- and final-state levels and multiplets are consistent to each other and
        to avoid later problems with the computations. An error message is issued if an inconsistency occurs,
        and nothing is returned otherwise.
"""
function  checkConsistentMultiplets(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet)
    initialSubshells      = initialMultiplet.levels[1].basis.subshells;             ni = length(initialSubshells)
    intermediateSubshells = intermediateMultiplet.levels[1].basis.subshells
    finalSubshells        = finalMultiplet.levels[1].basis.subshells
    
    if initialSubshells[1:end] == intermediateSubshells[1:ni]   &&
        intermediateSubshells   == finalSubshells
    else
        @show initialSubshells
        @show intermediateSubshells
        @show finalSubshells
        error("\nThe order of subshells must be equal for the initial-, intermediate and final states, \n"     *
                "and the same subshells must occur in the definition of the intermediate and final states. \n" *
                "Only the initial states can have less subshells; this limitation arises from the angular coefficients.")
    end
        
    return( nothing )
end


"""
`Dielectronic.checkOrbitalRepresentation(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet)`  
    ... to check (and analyze) that all high nl orbitals in these multiplets are properly represented on the given grid.
        The function prints for each symmetry block kappa the high-nl orbitals and checks that they are all bound.
        An error message is issued if this is not the case, and nothing is returned otherwise.
"""
function  checkOrbitalRepresentation(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet)
    
    println("\n  Check the energies and orbital representation:")
    println("\n  ----------------------------------------------")
    subshells = initialMultiplet.levels[1].basis.subshells;     basis = initialMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.    error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Initial occupied subshells:      $(subshells)")
    
    subshells = intermediateMultiplet.levels[1].basis.subshells;     basis = intermediateMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.    error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Intermediate occupied subshells: $(subshells)")
    
    subshells = finalMultiplet.levels[1].basis.subshells;     basis = finalMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.    error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Final occupied subshells:        $(subshells)")
        
    return( nothing )
end



"""
`Dielectronic.computeAmplitudesProperties(pathway::Dielectronic.Pathway, nm::Nuclear.Model, grid::Radial.Grid, 
                nrContinuum::Int64, settings::Dielectronic.Settings, hasCaptureChannels::Bool, lastCaptureChannels::Array{AutoIonization.Channel,1})` 
    ... to compute all amplitudes and properties of the given line; a line::Dielectronic.Pathway is returned for which the amplitudes and 
        properties have now been evaluated.
"""
function  computeAmplitudesProperties(pathway::Dielectronic.Pathway, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                        settings::Dielectronic.Settings, hasCaptureChannels::Bool, 
                                        lastCaptureChannels::Array{AutoIonization.Channel,1})
    rateA = 0.
    # println(">> pathway: $(pathway.initialLevel.index)--$(pathway.intermediateLevel.index)--$(pathway.finalLevel.index) ...")
    
    if hasCaptureChannels
        # Simply copy the results from previous computation of the same channels
        newcChannels = deepcopy(lastCaptureChannels)
        for cChannel in newcChannels
            rateA     = rateA + conj(cChannel.amplitude) * cChannel.amplitude
        end
    else
        newcChannels      = AutoIonization.Channel[];   contSettings = Continuum.Settings(false, nrContinuum)  
        initialLevel      = deepcopy(pathway.initialLevel)
        intermediateLevel = deepcopy(pathway.intermediateLevel)
        Defaults.setDefaults("relativistic subshell list", intermediateLevel.basis.subshells; printout=false)
        for cChannel in pathway.captureChannels
            newnLevel = Basics.generateLevelWithSymmetryReducedBasis(intermediateLevel, intermediateLevel.basis.subshells)
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(initialLevel, newnLevel.basis.subshells)
            newnLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, cChannel.kappa), newnLevel)
            cOrbital, phase  = Continuum.generateOrbitalForLevel(pathway.electronEnergy, Subshell(101, cChannel.kappa), 
                                                                    newiLevel, nm, grid, contSettings)
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, cChannel.symmetry, newiLevel)
            newcChannel = AutoIonization.Channel( cChannel.kappa, cChannel.symmetry, phase, Complex(0.))
            amplitude   = AutoIonization.amplitude(settings.augerOperator, cChannel, newcLevel, newnLevel, grid)
            rateA       = rateA + conj(amplitude) * amplitude
            newcChannel = AutoIonization.Channel( cChannel.kappa, cChannel.symmetry, phase, amplitude)
            push!( newcChannels, newcChannel)
        end
        println(">> pathway: $(pathway.initialLevel.index)--$(pathway.intermediateLevel.index)--$(pathway.finalLevel.index) ...")
    end
    #
    finalLevel        = deepcopy(pathway.finalLevel)
    intermediateLevel = deepcopy(pathway.intermediateLevel)
    Defaults.setDefaults("relativistic subshell list", intermediateLevel.basis.subshells; printout=false)
    
    newpChannels = PhotoEmission.Channel[];    rateC = 0.;    rateB = 0.
    for pChannel in pathway.photonChannels
        amplitude   = PhotoEmission.amplitude("emission", pChannel.multipole, pChannel.gauge, pathway.photonEnergy, 
                                                finalLevel, intermediateLevel, grid, display=false, printout=false)
        newpChannel = PhotoEmission.Channel( pChannel.multipole, pChannel.gauge, amplitude)
        push!( newpChannels, newpChannel)
        if       newpChannel.gauge == Basics.Coulomb     rateC = rateC + abs(amplitude)^2
        elseif   newpChannel.gauge == Basics.Babushkin   rateB = rateB + abs(amplitude)^2
        elseif   newpChannel.gauge == Basics.Magnetic    rateB = rateB + abs(amplitude)^2;   rateC = rateC + abs(amplitude)^2
        end
    end
    captureRate     = 2pi * rateA
    wa              = 8.0pi * Defaults.getDefaults("alpha") * pathway.photonEnergy / (Basics.twice(pathway.intermediateLevel.J) + 1) * 
                                                                                        (Basics.twice(pathway.finalLevel.J) + 1)
    wa              = wa / pi  ## modified for test with Xe^53+ (March/2024)
    photonRate      = EmProperty(wa * rateC, wa * rateB)  
    angularBeta     = EmProperty(-9., -9.)
    #  Factor due to UserGuide
    wa              = Defaults.convertUnits("kinetic energy to wave number: atomic units", pathway.electronEnergy)
    wa              = pi*pi / (wa*wa) * captureRate  * ##  2 * # factor 2 is not really clear.
                        ((Basics.twice(pathway.intermediateLevel.J) + 1) / (Basics.twice(pathway.initialLevel.J) + 1))
    #  Factor due to Tu et al. (Plasma Phys., 2016)
    ## wa              = pi*pi / 2. / pathway.electronEnergy * captureRate * 
    ##                   (Basics.twice(pathway.intermediateLevel.J) + 1) / (Basics.twice(pathway.initialLevel.J) + 1)
    reducedStrength = EmProperty(wa * photonRate.Coulomb, wa * photonRate.Babushkin)
    pathway = Dielectronic.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, pathway.electronEnergy, 
                                    pathway.photonEnergy, captureRate, photonRate, angularBeta, reducedStrength, newcChannels, newpChannels)
    return( pathway )
end


"""
`Dielectronic.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, 
                                    grid::Radial.Grid, settings::Dielectronic.Settings; output=true)`  
    ... to compute the dielectronic recombination amplitudes and all properties as requested by the given settings. 
        A list of pathways::Array{Dielectronic.Pathway,1} is returned.
"""
function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, 
                            grid::Radial.Grid, settings::Dielectronic.Settings; output=true)
    println("")
    printstyled("Dielectronic.computePathways(): The computation of dielectronic resonance strength, etc. starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    # First check that the initial, intermediate and final-state multiplets are consistent with each other to allow all requested computations
    Dielectronic.checkConsistentMultiplets(finalMultiplet, intermediateMultiplet, initialMultiplet)
    # Second, analyze that the high nl orbitals are properly represented with the given grid
    Dielectronic.checkOrbitalRepresentation(finalMultiplet, intermediateMultiplet, initialMultiplet)
    #
    pathways = Dielectronic.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    # Display all selected pathways before the computations start
    if  settings.printBefore    Dielectronic.displayPathways(pathways)    end
    # Determine maximum (electron) energy and check for consistency of the grid
    maxEnergy = 0.;   for  pathway in pathways   maxEnergy = max(maxEnergy, pathway.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    # Calculate all amplitudes and requested properties; simply copy if the captureChannels have been computed before
    newPathways = Dielectronic.Pathway[];    
    hasCaptureChannels = false;    lastCaptureIndices = (0,0);   lastCaptureChannels = AutoIonization.Channel[]
    for  pathway in pathways
        if  (pathway.initialLevel.index, pathway.intermediateLevel.index) == lastCaptureIndices   hasCaptureChannels = true
        else                                                                                      hasCaptureChannels = false
        end
        newPathway = Dielectronic.computeAmplitudesProperties(pathway, nm, grid, nrContinuum, settings, hasCaptureChannels, lastCaptureChannels) 
        push!( newPathways, newPathway)
        #
        lastCaptureIndices  = (newPathway.initialLevel.index, newPathway.intermediateLevel.index)
        lastCaptureChannels = newPathway.captureChannels
    end
    # 
    # Calculate all corresponding resonance
    resonances = Dielectronic.computeResonances(newPathways, settings)
    # Print all results to screen
    Dielectronic.displayResults(stdout, newPathways, settings)
    Dielectronic.displayResults(stdout, resonances,  settings)
    Dielectronic.displayRateCoefficients(stdout, resonances,  settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Dielectronic.displayResults(iostream, newPathways, settings)
                        Dielectronic.displayResults(iostream, resonances,  settings)
                        Dielectronic.displayRateCoefficients(iostream, resonances,  settings)    end
    #
    if    output    return( newPathways )
    else            return( nothing )
    end
end


"""
`Dielectronic.computeRateCoefficient(resonance::Dielectronic.Resonance, temp::Float64)`  
    ... computes for a delta-like resonance the DR rate coefficient alpha_d (i, Te) from the given resonance strength
        and temperature [K], and for both, Coulomb and Babushkin gauge. All values are directly returned in [cm^3/s].
        An alphaDR::EmProperty is returned.
"""
function computeRateCoefficient(resonance::Dielectronic.Resonance, temp::Float64)
    temp_au = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
    factor  = 4 / sqrt(2pi) * temp_au^(-3/2) * resonance.resonanceEnergy * exp(- resonance.resonanceEnergy/temp_au)
    alphaDR = factor * resonance.resonanceStrength 
    # Convert units into cm^3 / s
    factor  = Defaults.convertUnits("length: from atomic to fm", 1.0)^3 * 1.0e-39 * 
                Defaults.convertUnits("rate: from atomic to 1/s", 1.0) 
    alphaDR = factor * alphaDR
                
    return( alphaDR )
end


"""
`Dielectronic.computeResonances(pathways::Array{Dielectronic.Pathway,1}, settings::Dielectronic.Settings)`  
    ... to compute the data for all resonances (resonance lines) as defined by the given pathways and and settings. 
        A list of resonances::Array{Dielectronic.Resonance,1} is returned.
"""
function  computeResonances(pathways::Array{Dielectronic.Pathway,1}, settings::Dielectronic.Settings)
    # Determine all pre-defined resonances in pathways
    resonances = Dielectronic.Resonance[];    idxTuples = Tuple{Int64,Int64}[]
    for  pathway in pathways    idxTuples = union(idxTuples, [(pathway.initialLevel.index, pathway.intermediateLevel.index)] )    end
    # Determine the capture rate for each idxTuple
    captureRates = zeros( length(idxTuples) )
    for  (idx, idxTuple)  in  enumerate(idxTuples)
        for  pathway in pathways    
            if  idxTuple == (pathway.initialLevel.index, pathway.intermediateLevel.index)
                captureRates[idx] = pathway.captureRate
            end
        end
    end
    #
    # Determine the resonances
    electronEnergyShift = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    for  idxTuple in idxTuples
        iLevel = Level();    nLevel = Level();    resonanceEnergy = 0.;    reducedStrength = EmProperty(0., 0.)
        captureRate = 0.;    augerRate = 0.;      photonRate =  EmProperty(0., 0.)
        for  pathway in pathways    
            if  idxTuple == (pathway.initialLevel.index, pathway.intermediateLevel.index)
                iLevel            = pathway.initialLevel
                nLevel            = pathway.intermediateLevel
                resonanceEnergy   = pathway.intermediateLevel.energy - pathway.initialLevel.energy + electronEnergyShift
                captureRate       = pathway.captureRate 
                photonRate        = photonRate + pathway.photonRate
                reducedStrength   = reducedStrength + pathway.reducedStrength
            end
        end
        # Determine total Auger rates by summation over all idxTupels with the same m
        for  (idx, idxxTuple)  in  enumerate(idxTuples)
            if idxTuple[2] == idxxTuple[2]   augerRate = captureRates[idx]   end
        end
        # Correct the photon rate if requested
        if  typeof(settings.corrections) == Dielectronic.HydrogenicCorrections
            println(">>> Add hydrogenic corrections from n_low = $(settings.corrections.nDetailed+1) upwards")
            nDetailed     = settings.corrections.nDetailed
            Zeff          = settings.corrections.effectiveZ
            # Determine ni, li for the given resonance
            rydbergSubshs = Basics.extractRydbergSubshellList(nLevel, nDetailed+1, 1.0e-1)
            rydbergShells = Basics.extractNonrelativisticShellList(rydbergSubshs)
            if  length(rydbergShells) == 1  rShell = rydbergShells[1];   ni = rShell.n;   li = rShell.l
            else   error("Inappropriate number of Rydberg shells = $rydbergShells ")
            end
            # Compute and add hydrogenic rates A(ni,li --> nDetailed < n <= ni-1, li +- 1)
            hydrogenicRate = 0.
            for  nf = nDetailed+1:ni-1
                hydrogenicRate = hydrogenicRate + computeHydrogenicRate(ni, li, nf, li-1, Zeff) + 
                                                    computeHydrogenicRate(ni, li, nf, li+1, Zeff)
            end
            photonRate = photonRate + hydrogenicRate * settings.corrections.energyScaling
            println(">>> Total photon rate & total hydrogenic rate for A($ni, $li;  -> ...) = $photonRate, " * 
                    "$(hydrogenicRate * settings.corrections.energyScaling)")

        end
        
        resonanceStrength = reducedStrength / (augerRate + photonRate)
        push!( resonances, Dielectronic.Resonance( iLevel, nLevel, resonanceEnergy, resonanceStrength, captureRate, augerRate, photonRate) )
    end
    
    return( resonances )
end


"""
`Dielectronic.computeHydrogenicRate(ni::Int64, li::Int64, nf::Int64, lf::Int64,  Zeff::Float64)`
    ... to compute the nonrelativistic electric-dipole rate for the transition from shell ni,li --> nf,lf of a
        hydrogenic ion with effective charge Zeff. The recursion formulas by Infeld an Hull (1951) are used
        together with the absorption oscillator strength. This makes the overall formulation/computation rather
        obscure, unfortunately. Uses SpecialFunctions.logfactioial. A rate::Float64 [a.u.] is returned.
        This procedure has been worked out by Stefan Schippers (2023).
"""
function  computeHydrogenicRate(ni::Int64, li::Int64, nf::Int64, lf::Int64,  Zeff::Float64)
    # Compute A(n,l) coeffient in the recursion formulas
    function computeA(n::Int64, l::Int64)
        A = 0.0
        if n>l  &&  n*l > 0     A = sqrt( (n+l) * (n-l)) / (n*l)     end
        return ( A )
    end
    # Compute I(n,l; n',l') integral in the recursion formulas
    function computeI(n::Int64, l::Int64, np::Int64, lp::Int64)
        wi = 0.0
        #@show wi, n,l,np,lp
        if       l>=n  ||  l<0  || lp>n  ||  lp>=np  ||  lp<0  || abs(l-lp)!=1    wi = 0.0
        elseif   l == lp-1
            if  lp == n
                wi = (n+2)*log(4*n*np)+(np-n-2)*log(np-n)-(np+n+2)*log(np+n) + 
                        0.5 * ( SpecialFunctions.logfactorial(np+n) - SpecialFunctions.logfactorial(np-n-1) -
                                SpecialFunctions.logfactorial(2*n-1) )
                wi = 0.25* exp(wi)
            else
                wi = (2*lp+1)   * computeA(np, lp+1) * computeI(n,lp, np,lp+1) + computeA(n, lp+1) * computeI(n,lp+1, np,lp)
                wi = wi / (2*lp * computeA(n, lp))
            end
        elseif   l == lp+1
            wi = computeA(np, l+1) * computeI(n,l, np,l+1) + (2*l+1) * computeA(n, l+1) * computeI(n,l+1, np,l)
            wi = wi / (2*l * computeA(np, l))
        else
            error("Unexpected set of quantum number n=$n l=$l  np=$np  lp=$lp ")
        end
        return ( wi )
    end
    # Compute absorption oscillator strength
    function computeOsc(n::Int64, l::Int64, np::Int64, lp::Int64)
        # This oscillator strength is used in absorption
        wx = 0.0
        #@show wx, n, l, np, lp
        if       l == lp+1   wx = (1/n^2 - 1/np^2) * (lp+1) / (2*(lp+1) +1) * computeI(n,lp+1, np,lp)^2
        elseif   l == lp-1   wx = (1/n^2 - 1/np^2) *  lp    / (2*(lp-1) +1) * computeI(n,lp-1, np,lp)^2
        else     error("stop a")
        end
        return( wx / 3.0 )
    end
    #
    rate = 0.0;
    if  abs(li-lf)!=1  ||  ni <= nf  ||  li >= ni  ||  li<0  ||  lf >= nf  ||  lf<0   return( rate )   
    elseif  nf > 40    error("Don't use a recursive scheme ... but make a new implementation for ni = $ni ")
    end
    #
    if      lf == (li + 1)     rate = (2*(li+1) + 1) / (2*li+1) * computeOsc(nf, li+1, ni, li)
    elseif  lf == (li - 1)     rate = (2*(li-1) + 1) / (2*li+1) * computeOsc(nf, li-1, ni, li)
    else    error("stop b")
    end
    ##x alpha = 137.036
    ##x rate = rate * Zeff^4 / (2.0*alpha^3) * (1/nf^2 - 1/ni^2)^2
    rate = rate * Zeff^4 / 2.0 * Defaults.getDefaults("alpha")^3 * (1/nf^2 - 1/ni^2)^2
    ## println(">>>> Hydrogenic rate for A($ni, $li;  -> $nf, $lf) = $rate")

    return( rate )
end


"""
`Dielectronic.determineCaptureChannels(intermediateLevel::Level, initialLevel::Level, settings::Dielectronic.Settings)` 
    ... to determine a list of AutoIonization.Channel for a (Auger) capture transitions from the initial to an intermediate level, and by 
        taking into account the particular settings of for this computation;  an Array{AutoIonization.Channel,1} is returned.
"""
function determineCaptureChannels(intermediateLevel::Level, initialLevel::Level, settings::Dielectronic.Settings)
    channels = AutoIonization.Channel[];   
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)
    kappaList = AngularMomentum.allowedKappaSymmetries(symi, symn)
    for  kappa in kappaList
        push!( channels, AutoIonization.Channel(kappa, symn, 0., Complex(0.)) )
    end

    return( channels )  
end


"""
`Dielectronic.determinePhotonChannels(finalLevel::Level, intermediateLevel::Level, settings::Dielectronic.Settings)` 
    ... to determine a list of PhotoEmission.Channel for the photon transitions from the intermediate and to a final level, and by 
        taking into account the particular settings of for this computation;  an Array{PhotoEmission.Channel,1} is returned.
"""
function determinePhotonChannels(finalLevel::Level, intermediateLevel::Level, settings::Dielectronic.Settings)
    channels = PhotoEmission.Channel[];   
    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symn, mp, symf)
            hasMagnetic = false
            for  gauge in settings.gauges
                # Include further restrictions if appropriate
                if     string(mp)[1] == 'E'  &&   gauge == UseCoulomb      push!(channels, PhotoEmission.Channel(mp, Basics.Coulomb,   0.) )
                elseif string(mp)[1] == 'E'  &&   gauge == UseBabushkin    push!(channels, PhotoEmission.Channel(mp, Basics.Babushkin, 0.) )  
                elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)           push!(channels, PhotoEmission.Channel(mp, Basics.Magnetic,  0.) );
                                                    hasMagnetic = true; 
                end 
            end
        end
    end

    return( channels )  
end


"""
`Dielectronic.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                settings::Dielectronic.Settings)`  
    ... to determine a list of dielectronic-recombination pathways between the levels from the given initial-, intermediate- and 
        final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{Dielectronic.Pathway,1} is returned. Apart from the level specification, all physical properties are set to zero 
        during the initialization process.  
"""
function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, settings::Dielectronic.Settings)
    pathways = Dielectronic.Pathway[]
    electronEnergyShift = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    photonEnergyShift   = Defaults.convertUnits("energy: to atomic", settings.photonEnergyShift)
    #
    for  iLevel  in  initialMultiplet.levels
        for  nLevel  in  intermediateMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)
                    eEnergy = nLevel.energy - iLevel.energy + electronEnergyShift
                    pEnergy = nLevel.energy - fLevel.energy + photonEnergyShift
                    if  pEnergy < 0.   ||   eEnergy < 0.    continue    end
                    cChannels = Dielectronic.determineCaptureChannels(nLevel, iLevel, settings) 
                    pChannels = Dielectronic.determinePhotonChannels( fLevel, nLevel, settings) 
                    push!( pathways, Dielectronic.Pathway(iLevel, nLevel, fLevel, eEnergy, pEnergy, 0., EmProperty(0., 0.), 
                                                            EmProperty(0., 0.), EmProperty(0., 0.), cChannels, pChannels) )
                end
            end
        end
    end
    return( pathways )
end


"""
`Dielectronic.displayPathways(pathways::Array{Dielectronic.Pathway,1})`  
    ... to display a list of pathways and channels that have been selected due to the prior settings. A neat table of all selected 
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayPathways(pathways::Array{Dielectronic.Pathway,1})
    nx = 180
    println(" ")
    println("  Selected dielectronic-recombination pathways:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "     ";   sb = "     "
    sa = sa * TableStrings.center(23, "Levels"; na=4);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=4);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=3);
    sa = sa * TableStrings.center(26, "Energies  " * TableStrings.inUnits("energy"); na=5);              
    sb = sb * TableStrings.center(26, "electron        photon "; na=5)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=5)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "   "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy))   * "    "
        kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
        for  cChannel in pathway.captureChannels
            for  pChannel in pathway.photonChannels
                push!( kappaMultipoleSymmetryList, (cChannel.kappa, pChannel.multipole, pChannel.gauge, cChannel.symmetry) )
            end
        end
        wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
        if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
        end
    end
    println("  ", TableStrings.hLine(nx))
    println("\n>> A total of $(length(pathways)) dielectronic-recombination pathways will be calculated. \n")
    #
    return( nothing )
end


"""
`Dielectronic.displayResults(stream::IO, pathways::Array{Dielectronic.Pathway,1}, settings::Dielectronic.Settings)`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, pathways::Array{Dielectronic.Pathway,1}, settings::Dielectronic.Settings)
    nx = 150
    println(stream, " ")
    println(stream, "  Partial (Auger) capture and radiative decay rates:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);
    sa = sa * TableStrings.center(38, "Energies  " * TableStrings.inUnits("energy"); na=4);              
    sb = sb * TableStrings.center(38, "electron        m--i       photon  "; na=1)
    sa = sa * TableStrings.center(10, "Multipoles"; na=6);        sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(36, "Rates  " * TableStrings.inUnits("rate"); na=2);   
    sb = sb * TableStrings.center(36, "(Auger) capture    Cou--photon--Bab";        na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=3)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        en_mi = pathway.intermediateLevel.energy - pathway.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", en_mi))                  * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy))   * "    "
        multipoles = EmMultipole[]
        for  pch in pathway.photonChannels
            multipoles = push!( multipoles, pch.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles)     * "                   "
        sa = sa * TableStrings.flushleft(16, mpString[1:15];  na=0)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", pathway.captureRate))            * "     "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", pathway.photonRate.Coulomb))     * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", pathway.photonRate.Babushkin))   * "  "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    #
    nx = 137
    println(stream, " ")
    println(stream, "  Partial resonance strength:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);
    sa = sa * TableStrings.center(38, "Energies  " * TableStrings.inUnits("energy"); na=4);              
    sb = sb * TableStrings.center(38, "electron        m--i       photon  "; na=1)
    sa = sa * TableStrings.center(10, "Multipoles"; na=5);        sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(26, "S * Gamma_m  " * TableStrings.inUnits("reduced strength"); na=2);   
    sb = sb * TableStrings.center(26, " Cou -- photon -- Bab";        na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=3)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        en_mi = pathway.intermediateLevel.energy - pathway.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", en_mi))                  * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy))   * "    "
        multipoles = EmMultipole[]
        for  pch in pathway.photonChannels
            multipoles = push!( multipoles, pch.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles)                          * "                   "
        sa = sa * TableStrings.flushleft(16, mpString[1:15];  na=0)
        wa = Defaults.convertUnits("strength: from atomic", 1.0) * Defaults.convertUnits("energy: from atomic", 1.0)
        sa = sa * @sprintf("%.4e", wa * pathway.reducedStrength.Coulomb)     * "     "
        sa = sa * @sprintf("%.4e", wa * pathway.reducedStrength.Babushkin)   * "     "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
    + (stream::IO, resonances::Array{Dielectronic.Resonance,1}, settings::Dielectronic.Settings)`  
    ... to list all results for the resonances. A neat table is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, resonances::Array{Dielectronic.Resonance,1}, settings::Dielectronic.Settings)
    nx = 160
    println(stream, " ")
    println(stream, "  Total Auger rates, radiative rates and resonance strengths:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(42, "Auger rate     Cou -- rad. rates -- Bab"; na=1);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("rate"); na=1)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("rate"); na=0)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("rate"); na=6)
    sa = sa * TableStrings.center(30, "Cou -- res. strength -- Bab"; na=3);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("strength");  na=0)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("strength");  na=2)
    sa = sa * TableStrings.center(18, "Widths Gamma_m"; na=2);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("energy"); na=6)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  resonance in resonances
        sa  = "";      isym = LevelSymmetry( resonance.initialLevel.J,      resonance.initialLevel.parity)
                        msym = LevelSymmetry( resonance.intermediateLevel.J, resonance.intermediateLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(resonance.initialLevel.index, resonance.intermediateLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", resonance.resonanceEnergy))          * "      "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.augerRate))                  * "      "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.photonRate.Coulomb))         * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.photonRate.Babushkin))       * "        "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", resonance.resonanceStrength.Coulomb))    * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", resonance.resonanceStrength.Babushkin))  * "     "
        wa = resonance.augerRate + resonance.photonRate.Coulomb 
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wa))                                 * "   "
        wa = resonance.augerRate + resonance.photonRate.Babushkin 
        sa = sa * "(" * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wa)) * ")"                     * "   "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`Dielectronic.displayRateCoefficients(stream::IO, resonances::Array{Dielectronic.Resonance,1}, settings::Dielectronic.Settings)`  
    ... to list, if settings.calcRateAlpha, all rate coefficients for the selected temperatures. Both, the individual as well as
        the total DR plasma rate coefficients are printed in neat tables, though nothing is returned otherwise.
"""
function  displayRateCoefficients(stream::IO, resonances::Array{Dielectronic.Resonance,1}, settings::Dielectronic.Settings)
    ntemps = length(settings.temperatures)
    if  !settings.calcRateAlpha  ||  ntemps == 0     return(nothing)     end
    #
    nx = 54 + 17 * min(ntemps, 7)
    println(stream, " ")
    println(stream, "  Rate coefficients for delta-like resonances [cm^3/s]:        ... all results in Babushkin gauge")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    for  nt = 1:min(ntemps, 7)
        sa = sa * TableStrings.center(14, "T = " * @sprintf("%.2e", settings.temperatures[nt]); na=3);       
        sb = sb * TableStrings.center(14, "[K]"; na=3)
    end
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  resonance in resonances
        sa  = "";      isym = LevelSymmetry( resonance.initialLevel.J,      resonance.initialLevel.parity)
                        msym = LevelSymmetry( resonance.intermediateLevel.J, resonance.intermediateLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(resonance.initialLevel.index, resonance.intermediateLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", resonance.resonanceEnergy))    * "      "
        for  nt = 1:min(ntemps, 7)
            alphaDR      = Dielectronic.computeRateCoefficient(resonance, settings.temperatures[nt])
            sa = sa * @sprintf("%.4e", alphaDR.Babushkin)  * "       "
        end
        println(stream, sa)
    end
    #
    println(stream, "  ")
    sa = "       alpha^DR (T, i; Coulomb gauge):                      " 
    sb = "       alpha^DR (T, i; Babushkin gauge):                    " 
    for  nt = 1:min(ntemps, 7) 
        alphaDRtotal = EmProperty(0.)
        for  resonance in resonances
            alphaDR      = Dielectronic.computeRateCoefficient(resonance, settings.temperatures[nt])
            alphaDRtotal = alphaDRtotal + alphaDR
        end
        sa = sa * @sprintf("%.4e", alphaDRtotal.Coulomb)    * "       "
        sb = sb * @sprintf("%.4e", alphaDRtotal.Babushkin)  * "       "
    end
    println(stream, sa);    println(stream, sb)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`Dielectronic.extractRateCoefficients(resonances::Array{Dielectronic.Resonance,1}, settings::Dielectronic.Settings)`  
    ... to extract, if settings.calcRateAlpha, the total DR rate coefficients for all temperatures. 
        A list of total rate coefficients [cm^3/s] alphaDR::Array{EmProperty,1} is returned.
"""
function  extractRateCoefficients(resonances::Array{Dielectronic.Resonance,1}, settings::Dielectronic.Settings)
    ntemps = length(settings.temperatures);     alphaDRs = EmProperty[]
    if  !settings.calcRateAlpha  ||  ntemps == 0     return( EmProperty[] )     end
    #
    for  nt = 1:ntemps
        alphaDRtotal = EmProperty(0.)
        for  resonance in resonances
            alphaDR      = Dielectronic.computeRateCoefficient(resonance, settings.temperatures[nt])
            alphaDRtotal = alphaDRtotal + alphaDR
        end
        push!(alphaDRs, alphaDRtotal)
    end
    #
    return( alphaDRs )
end


"""
`Dielectronic.isResonanceToBeExcluded(level::Level, refLevel, rSelection::ResonanceSelection)`  
    returns true, if level is to be excluded from the valid resonances, and false otherwise.
    It returns false if the ResonanceSelection() is inactive or if level belongs to the selected resoances.
    It is true only of ResonanceSelection() is active but the level does not belong to the selected resonances.
"""
function  isResonanceToBeExcluded(level::Level, refLevel, rSelection::ResonanceSelection)
    wa = false
    if !rSelection.active    return( wa )
        # This is the standard case if no additional limitations are specified by the user
    else
        # Analyze of whether level belongs to the selected resonances; it determines of whether level has 
        # electrons in either toShells or intoShells, and if there is one electron less in the fromShells
        fromwb    = false;    towb    = false;    intowb    = true
        confList  = Basics.extractNonrelativisticConfigurations(level.basis)
        occShells = Basics.extractShellList(confList) 
        for  shell in rSelection.toShells   
            if  shell in occShells   towb   = true;   break    end
        end
        for  shell in occShells   
            if  !(shell in rSelection.intoShells)   intowb = false;   break    end
        end
        # Compare the occupation of the given resonance level with those of the refLevel; there should be (at least)
        # one electron more in the shells of leadingConfig than in the same shells of level
        refConfig = Basics.extractLeadingConfiguration(refLevel)
        levConfig = confList[1];   NoElectrons = 0
        for  (k,v) in refConfig.shells
            NoElectrons = NoElectrons + levConfig.shells[k]
        end
        @show NoElectrons, refConfig.NoElectrons
        if  NoElectrons < refConfig.NoElectrons     fromwb = true   end
        
        wa = !(fromwb &&  towb  && intowb)
    end
        
    
    return( wa )
end

end # module
