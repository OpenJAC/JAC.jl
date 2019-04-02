
"""
`module  JAC.Dielectronic`  ... a submodel of JAC that contains all methods for computing dielectronic recombination properties between 
                                some initial, intermediate and final-state multiplets; it is using JAC, JAC.ManyElectron, JAC.Radial, 
                                JAC.Radiative, JAC.Auge.
"""
module Dielectronic

    using Printf, JAC, JAC.ManyElectron, JAC.Radial, JAC.Radiative, JAC.Auger
    global JAC_counter = 0


    """
    `struct  Dielectronic.Settings`  ... defines a type for the details and parameters of computing dielectronic recombination pathways.

        + multipoles              ::Array{EmMultipoles}   ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}       ... Specifies the gauges to be included into the computations.
        + printBeforeComputation  ::Bool                  ... True, if all energies and pathways are printed before their evaluation.
        + selectPathways          ::Bool                               ... True if particular pathways are selected for the computations.
        + selectedPathways        ::Array{Tuple{Int64,Int64,Int64},1}  ... List of list of pathways, given by tupels (inital, inmediate, final).
        + electronEnergyShift     ::Float64               ... An overall energy shift for all electron energies (i.e. from the initial to
                                                              the resonance levels.
        + photonEnergyShift       ::Float64               ... An overall energy shift for all photon energies (i.e. from the resonance to
                                                              the final levels.
        + mimimumPhotonEnergy     ::Float64               ... minimum transition energy for which photon transitions are included into the
                                                              evaluation.
        + augerOperator           ::String                ... Auger operator that is to be used for evaluating the Auger amplitudes: allowed values
                                                              are: "Coulomb", "Breit", "Coulomb+Breit"
    """
    struct Settings 
        multipoles                ::Array{EmMultipole,1}
        gauges                    ::Array{UseGauge}
        printBeforeComputation    ::Bool
        selectPathways            ::Bool
        selectedPathways          ::Array{Tuple{Int64,Int64,Int64},1}
        electronEnergyShift       ::Float64
        photonEnergyShift         ::Float64
        mimimumPhotonEnergy       ::Float64
        augerOperator             ::String
        
    end 


    """
    `JAC.Dielectronic.Settings()`  ... constructor for the default values of dielectronic recombination pathway computations.
    """
    function Settings()
        Settings([E1], UseGauge[], false, false, Tuple{Int64,Int64,Int64}[], 0., 0., 0., "Coulomb")
    end


    """
    `Base.show(io::IO, settings::Dielectronic.Settings)`  ... prepares a proper printout of the variable settings::Dielectronic.Settings.
    """
    function Base.show(io::IO, settings::Dielectronic.Settings) 
        println(io, "multipoles:                 $(settings.multipoles)  ")
        println(io, "use-gauges:                 $(settings.gauges)  ")
        println(io, "printBeforeComputation:     $(settings.printBeforeComputation)  ")
        println(io, "selectPathways:             $(settings.selectPathways)  ")
        println(io, "selectedPathways:           $(settings.selectedPathways)  ")
        println(io, "electronEnergyShift:        $(settings.electronEnergyShift)  ")
        println(io, "photonEnergyShift:          $(settings.photonEnergyShift)  ")
        println(io, "mimimumPhotonEnergy:        $(settings.mimimumPhotonEnergy)  ")
        println(io, "augerOperator:              $(settings.augerOperator)  ")
    end


    #===
   """
    `struct  Dielectronic.Channel`  ... defines a type for a single dielectronic recombination channel that specifies all quantum numbers, 
                                        phases and amplitudes.

         + augerChannel      ::JAC.Auger.Channel        ... Channel that describes the capture, i.e. inverse Auger process.
         + radiativeChannel  ::JAC.Radiative.Channel    ... Channel that describes the stabilization, i.e. photon emission.
    """
    struct  Channel
        augerChannel         ::JAC.Auger.Channel
        radiativeChannel     ::JAC.Radiative.Channel
    end   ==#


    """
    `struct  Dielectronic.Pathway`  ... defines a type for a dielectronic recombination pathways that may include the definition of channels and 
                                        their corresponding amplitudes.

        + initialLevel      ::Level                   ... initial-(state) level
        + intermediateLevel ::Level                   ... intermediate-(state) level
        + finalLevel        ::Level                   ... final-(state) level
        + electronEnergy    ::Float64                 ... energy of the (incoming, captured) electron
        + photonEnergy      ::Float64                 ... energy of the (emitted) photon
        + captureRate       ::Float64                 ... rate for the electron capture (Auger rate)
        + photonRate        ::EmProperty              ... rate for the photon emission
        + angularBeta       ::EmProperty              ... beta parameter of the photon emission
        + resonanceStrength ::EmProperty              ... partial resonance strength of this pathway
        + hasChannels       ::Bool                    ... Determines whether the individual channels are defined in terms of their possible
                                                          Auger and radiative channels, or not.
        + captureChannels   ::Array{JAC.Auger.Channel,1}      ... List of |i> -->  |n>   dielectronic (Auger) capture channels.
        + photonChannels    ::Array{JAC.Radiative.Channel,1}  ... List of |n> -->  |f>   radiative stabilization channels.
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
        resonanceStrength   ::EmProperty
        hasChannels         ::Bool
        captureChannels     ::Array{JAC.Auger.Channel,1} 
        photonChannels      ::Array{JAC.Radiative.Channel,1} 
    end 


    """
    `JAC.Dielectronic.Pathway()`  ... constructor for an 'empty' instance of a dielectronic recombination pathway between a specified 
                                      initial, intermediate and final level.
    """
    function Pathway()
        em = EmProperty(0., 0.)
        Pathway(initialLevel, intermediateLevel, finalLevel, 0., 0., 0., em, em, em, false, Auger.Channel[], Radiative.Channel[])
    end


    """
    `Base.show(io::IO, pathway::Dielectronic.Pathway)`  ... prepares a proper printout of the variable pathway::Dielectronic.Pathway.
    """
    function Base.show(io::IO, pathway::Dielectronic.Pathway) 
        println(io, "initialLevel:               $(pathway.initialLevel)  ")
        println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
        println(io, "finalLevel:                 $(pathway.finalLevel)  ")
        println(io, "electronEnergy:             $(pathway.electronEnergy)  ")
        println(io, "photonEnergy:               $(pathway.photonEnergy)  ")
        println(io, "captureRate:                $(pathway.captureRate)  ")
        println(io, "photonRate:                 $(pathway.photonRate)  ")
        println(io, "angularBeta:                $(pathway.angularBeta)  ")
        println(io, "resonanceStrength:          $(pathway.resonanceStrength)  ")
        println(io, "hasChannels:                $(pathway.hasChannels)  ")
        println(io, "captureChannels:            $(pathway.captureChannels)  ")
        println(io, "photonChannels:             $(pathway.photonChannels)  ")
    end


    """
    `struct  Dielectronic.Resonance`  ... defines a type for a dielectronic resonance as defined by a given initial and resonance level but 
                                          by summing over all final levels

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
    `JAC.Dielectronic.Resonance()`  ... constructor for an 'empty' instance of a dielectronic resonance as defined by a given initial 
                                        and resonance level but by summing over all final levels.
    """
    function Resonance()
        em = EmProperty(0., 0.)
        Resonance(initialLevel, intermediateLevel, 0., em, 0., 0., em)
    end


    """
    `Base.show(io::IO, resonance::Dielectronic.Resonance)`  ... prepares a proper printout of the variable resonance::Dielectronic.Resonance.
    """
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
    `JAC.Dielectronic.computeAmplitudesProperties(pathway::Dielectronic.Pathway, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                                                  nrContinuum::Int64, settings::Dielectronic.Settings)` 
        ... to compute all amplitudes and properties of the given line; a line::Dielectronic.Pathway is returned for which the amplitudes and 
            properties have now been evaluated.
    """
    function  computeAmplitudesProperties(pathway::Dielectronic.Pathway, nm::JAC.Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                          settings::Dielectronic.Settings)
        newcChannels = Auger.Channel[];   contSettings = JAC.Continuum.Settings(false, nrContinuum)
        for cChannel in pathway.captureChannels
            newnLevel   = JAC.generateLevelWithSymmetryReducedBasis(pathway.intermediateLevel)
            newnLevel   = JAC.generateLevelWithExtraSubshell(Subshell(101, cChannel.kappa), newnLevel)
            newiLevel   = JAC.generateLevelWithSymmetryReducedBasis(pathway.initialLevel)
            cOrbital, phase  = JAC.Continuum.generateOrbitalForLevel(pathway.electronEnergy, Subshell(101, cChannel.kappa), newiLevel, nm, grid, contSettings)
            newcLevel   = JAC.generateLevelWithExtraElectron(cOrbital, cChannel.symmetry, newiLevel)
            newcChannel = Auger.Channel( cChannel.kappa, cChannel.symmetry, phase, Complex(0.))
            amplitude   = JAC.Auger.amplitude(settings.augerOperator, cChannel, newnLevel, newcLevel, grid)
            newcChannel = Auger.Channel( cChannel.kappa, cChannel.symmetry, phase, amplitude)
            push!( newcChannels, newcChannel)
        end
        #
        newpChannels = Radiative.Channel[]
        for pChannel in pathway.photonChannels
            amplitude   = JAC.Radiative.amplitude("emission", pChannel.multipole, pChannel.gauge, pathway.photonEnergy, 
                                                  pathway.finalLevel, pathway.intermediateLevel, grid)
            newpChannel = Radiative.Channel( pChannel.multipole, pChannel.gauge, amplitude)
            push!( newpChannels, newpChannel)
        end
        captureRate = 1.0
        photonRate     = EmProperty(-1., -1.);   angularBeta    = EmProperty(-1., -1.)
        polarizationP0 = EmProperty(-1., -1.);   polarizationP2 = EmProperty(-1., -1.);   resonanceStrength = EmProperty(-1., -1.)
        pathway = Dielectronic.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, pathway.electronEnergy, 
                                        pathway.photonEnergy, captureRate, photonRate, angularBeta, resonanceStrength, true, newcChannels, newpChannels)
        return( pathway )
    end


    """
    `JAC.Dielectronic.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, nm::JAC.Nuclear.Model, 
                                      grid::Radial.Grid, settings::Dielectronic.Settings; output=true)`  ... to compute the dielectronic 
         recombination amplitudes and all properties as requested by the given settings. A list of pathways::Array{Dielectronic.Pathway,1} is 
         returned.
    """
    function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, 
                              settings::Dielectronic.Settings; output=true)
        println("")
        printstyled("JAC.Dielectronic.computePathways(): The computation of dielectronic resonance strength, etc. starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        pathways = JAC.Dielectronic.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
        # Display all selected pathways before the computations start
        if  settings.printBeforeComputation    JAC.Dielectronic.displayPathways(pathways)    end
        # Determine maximum (electron) energy and check for consistency of the grid
        maxEnergy = 0.;   for  pathway in pathways   maxEnergy = max(maxEnergy, pathway.electronEnergy)   end
        nrContinuum = JAC.Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newPathways = Dielectronic.Pathway[]
        for  pathway in pathways
            newPathway = JAC.Dielectronic.computeAmplitudesProperties(pathway, nm, grid, nrContinuum, settings) 
            push!( newPathways, newPathway)
        end
        # 
        # Calculate all corresponding resonance
        resonances = JAC.Dielectronic.computeResonances(newPathways, settings)
        # Print all results to screen
        JAC.Dielectronic.displayResults(stdout, newPathways, settings)
        JAC.Dielectronic.displayResults(stdout, resonances,  settings)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary   JAC.Dielectronic.displayResults(iostream, newPathways, settings)
                           JAC.Dielectronic.displayResults(iostream, resonances,  settings)    end
        #
        if    output    return( newPathways )
        else            return( nothing )
        end
    end


    """
    `JAC.Dielectronic.computeResonances(pathways::Array{Dielectronic.Pathway,1}, settings::Dielectronic.Settings)`  
        ... to compute the data for all resonances (resonance lines) as defined by the given pathways and and settings. 
            A list of resonances::Array{Dielectronic.Resonance,1} is returned.
    """
    function  computeResonances(pathways::Array{Dielectronic.Pathway,1}, settings::Dielectronic.Settings)
        # Determine all pre-defined resonances in pathways
        resonances = Dielectronic.Resonance[];    idxTuples = Tuple{Int64,Int64}[]
        for  pathway in pathways    idxTuples = union(idxTuples, [(pathway.initialLevel.index, pathway.intermediateLevel.index)] )    end
        ##x println("idxTuples = $idxTuples")
        for  idxTuple in idxTuples
            iLevel = Level();    nLevel = Level();    resonanceEnergy = 0.;    resonanceStrength = EmProperty(0., 0.)
            captureRate = 0.;    augerRate = 0.;    photonRate =  EmProperty(0., 0.)
            for  pathway in pathways    
                ##x println("idxTuple = $idxTuple")
                if  idxTuple == (pathway.initialLevel.index, pathway.intermediateLevel.index)
                    iLevel            = pathway.initialLevel
                    nLevel            = pathway.intermediateLevel
                    resonanceEnergy   = pathway.intermediateLevel.energy - pathway.initialLevel.energy
                    resonanceStrength = resonanceStrength ## + pathway.resonanceStrength 
                    captureRate       = pathway.captureRate 
                    augerRate         = augerRate  ## + pathway.captureRate 
                    photonRate        = photonRate ## + pathway.photonRate
                end
            end
            push!( resonances, Dielectronic.Resonance( iLevel, nLevel, resonanceEnergy, resonanceStrength, captureRate, augerRate, photonRate) )
        end
        
        return( resonances )
    end


    """
    `JAC.Dielectronic.determineCaptureChannels(intermediateLevel::Level, initialLevel::Level, settings::Dielectronic.Settings)` 
        ... to determine a list of Auger.Channel for a (Auger) capture transitions from the initial to an intermediate level, and by 
            taking into account the particular settings of for this computation;  an Array{Auger.Channel,1} is returned.
    """
    function determineCaptureChannels(intermediateLevel::Level, initialLevel::Level, settings::Dielectronic.Settings)
        channels = Auger.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)
        kappaList = JAC.AngularMomentum.allowedKappaSymmetries(symi, symn)
        for  kappa in kappaList
            push!( channels, Auger.Channel(kappa, symn, 0., Complex(0.)) )
        end

        return( channels )  
    end


    """
    `JAC.Dielectronic.determinePhotonChannels(finalLevel::Level, intermediateLevel::Level, settings::Dielectronic.Settings)` 
        ... to determine a list of Radiative.Channel for the photon transitions from the intermediate and to a final level, and by 
            taking into account the particular settings of for this computation;  an Array{Radiative.Channel,1} is returned.
    """
    function determinePhotonChannels(finalLevel::Level, intermediateLevel::Level, settings::Dielectronic.Settings)
        channels = Radiative.Channel[];   
        symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp in settings.multipoles
            if   JAC.AngularMomentum.isAllowedMultipole(symn, mp, symf)
                hasMagnetic = false
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp)[1] == 'E'  &&   gauge == JAC.UseCoulomb      push!(channels, Radiative.Channel(mp, JAC.Coulomb,   0.) )
                    elseif string(mp)[1] == 'E'  &&   gauge == JAC.UseBabushkin    push!(channels, Radiative.Channel(mp, JAC.Babushkin, 0.) )  
                    elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)               push!(channels, Radiative.Channel(mp, JAC.Magnetic,  0.) );
                                                        hasMagnetic = true; 
                    end 
                end
            end
        end

        return( channels )  
    end


    """
    `JAC.Dielectronic.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                        settings::Dielectronic.Settings)`  ... to determine a list of dielectronic-recombination pathways 
         between the levels from the given initial-, intermediate- and final-state multiplets and by taking into account the particular 
         selections and settings for this computation; an Array{Dielectronic.Pathway,1} is returned. Apart from the level specification, 
         all physical properties are set to zero during the initialization process.  
    """
    function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                settings::Dielectronic.Settings)
        if    settings.selectPathways    selectPathways = true;    
              selectedPathways = JAC.determineSelectedPathways(settings.selectedPathways, initialMultiplet, intermediateMultiplet, finalMultiplet)
        else                             selectPathways = false
        end
    
        pathways = Dielectronic.Pathway[]
        for  i = 1:length(initialMultiplet.levels)
            for  n = 1:length(intermediateMultiplet.levels)
                for  f = 1:length(finalMultiplet.levels)
                    if  selectPathways  &&  !((i,n,f) in selectedPathways )    continue   end
                    eEnergy = intermediateMultiplet.levels[n].energy - initialMultiplet.levels[i].energy
                    pEnergy = intermediateMultiplet.levels[n].energy - finalMultiplet.levels[f].energy
                    if  pEnergy < 0.   ||   eEnergy < 0.    continue    end

                    cChannels = JAC.Dielectronic.determineCaptureChannels(intermediateMultiplet.levels[n], initialMultiplet.levels[i], settings) 
                    pChannels = JAC.Dielectronic.determinePhotonChannels(finalMultiplet.levels[f], intermediateMultiplet.levels[n], settings) 
                    push!( pathways, Dielectronic.Pathway(initialMultiplet.levels[i], intermediateMultiplet.levels[n], finalMultiplet.levels[f], 
                                                          eEnergy, pEnergy, 0., EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), 
                                                          true, cChannels, pChannels) )
                end
            end
        end
        return( pathways )
    end


    """
    `JAC.Dielectronic.displayPathways(pathways::Array{Dielectronic.Pathway,1})`  ... to display a list of pathways and channels that have been 
         selected due to the prior settings. A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayPathways(pathways::Array{Dielectronic.Pathway,1})
        println(" ")
        println("  Selected dielectronic-recombination pathways:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(180))
        sa = "     ";   sb = "     "
        sa = sa * JAC.TableStrings.center(23, "Levels"; na=4);            sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=4);          
        sa = sa * JAC.TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=3);
        sa = sa * JAC.TableStrings.center(26, "Energies  " * JAC.TableStrings.inUnits("energy"); na=5);              
        sb = sb * JAC.TableStrings.center(26, "electron        photon "; na=5)
        sa = sa * JAC.TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
        sb = sb * JAC.TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(180)) 
        #   
        for  pathway in pathways
            sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=5)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", pathway.electronEnergy)) * "   "
            sa = sa * @sprintf("%.6e", JAC.convert("energy: from atomic", pathway.photonEnergy))   * "    "
            kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
            for  cChannel in pathway.captureChannels
                for  pChannel in pathway.photonChannels
                    push!( kappaMultipoleSymmetryList, (cChannel.kappa, pChannel.multipole, pChannel.gauge, cChannel.symmetry) )
                end
            end
            wa = JAC.TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
            if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
            for  i = 2:length(wa)
                sb = JAC.TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
        end
        println("  ", JAC.TableStrings.hLine(180))
        #
        return( nothing )
    end


    """
    `JAC.Dielectronic.displayResults(stream::IO, pathways::Array{Dielectronic.Pathway,1}, settings::Dielectronic.Settings)`  
         ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing 
         is returned otherwise.
    """
    function  displayResults(stream::IO, pathways::Array{Dielectronic.Pathway,1}, settings::Dielectronic.Settings)
        println(stream, " ")
        println(stream, "  Partial (Auger) capture and radiative decay rates:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(150))
        sa = "    ";   sb = "    "
        sa = sa * JAC.TableStrings.center(23, "Levels"; na=2);            sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=2);          
        sa = sa * JAC.TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=2);
        sa = sa * JAC.TableStrings.center(38, "Energies  " * JAC.TableStrings.inUnits("energy"); na=4);              
        sb = sb * JAC.TableStrings.center(38, "electron        m--i       photon  "; na=1)
        sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=6);        sb = sb * JAC.TableStrings.hBlank(17)
        sa = sa * JAC.TableStrings.center(36, "Rates  " * JAC.TableStrings.inUnits("rate"); na=2);   
        sb = sb * JAC.TableStrings.center(36, "(Auger) capture    Cou--photon--Bab";        na=2)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(150)) 
        #   
        for  pathway in pathways
            sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=3)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            en_mi = pathway.intermediateLevel.energy - pathway.initialLevel.energy
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", pathway.photonEnergy))   * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", en_mi))                  * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", pathway.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  pch in pathway.photonChannels
                multipoles = push!( multipoles, pch.multipole)
            end
            multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles)     * "                   "
            sa = sa * JAC.TableStrings.flushleft(16, mpString[1:16];  na=2)
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", pathway.captureRate))            * "     "
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", pathway.photonRate.Coulomb))     * "  "
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", pathway.photonRate.Babushkin))   * "  "
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(150))
        #
        #
        #
        println(stream, " ")
        println(stream, "  Partial resonance strength:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(150))
        sa = "    ";   sb = "    "
        sa = sa * JAC.TableStrings.center(23, "Levels"; na=2);            sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=2);          
        sa = sa * JAC.TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=2);
        sa = sa * JAC.TableStrings.center(38, "Energies  " * JAC.TableStrings.inUnits("energy"); na=4);              
        sb = sb * JAC.TableStrings.center(38, "electron        m--i       photon  "; na=1)
        sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=6);        sb = sb * JAC.TableStrings.hBlank(17)
        sa = sa * JAC.TableStrings.center(26, "resonance strength  " * JAC.TableStrings.inUnits("rate"); na=2);   
        sb = sb * JAC.TableStrings.center(26, " Cou -- photon -- Bab";        na=2)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(150)) 
        #   
        for  pathway in pathways
            sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=3)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            en_mi = pathway.intermediateLevel.energy - pathway.initialLevel.energy
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", pathway.photonEnergy))   * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", en_mi))                  * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", pathway.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  pch in pathway.photonChannels
                multipoles = push!( multipoles, pch.multipole)
            end
            multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles)     * "                   "
            sa = sa * JAC.TableStrings.flushleft(16, mpString[1:16];  na=1)
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", pathway.resonanceStrength.Coulomb))     * "     "
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", pathway.resonanceStrength.Babushkin))   * "     "
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(150))
        #
        return( nothing )
    end


    """
        + (stream::IO, resonances::Array{Dielectronic.Resonance,1}, settings::Dielectronic.Settings)`  
         ... to list all results for the resonances. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(stream::IO, resonances::Array{Dielectronic.Resonance,1}, settings::Dielectronic.Settings)
        println(stream, " ")
        println(stream, "  Total Auger rates, radiative rates and resonance strengths:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(150))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(14, "Energy"   ; na=4);               
        sb = sb * JAC.TableStrings.center(14,JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(44, "Auger rate    Cou -- rad. rates -- Bab"; na=2);       
        sb = sb * JAC.TableStrings.center(16, JAC.TableStrings.inUnits("rate"); na=2)
        sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("rate"); na=2)
        sb = sb * JAC.TableStrings.center(12, JAC.TableStrings.inUnits("rate"); na=2)
        sa = sa * JAC.TableStrings.center(32, "Cou -- res. strength -- Bab"; na=2);       
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("rate"); na=2)
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("rate"); na=2)
        sa = sa * JAC.TableStrings.center(18, "Widths Gamma_m"; na=2);       
        sb = sb * JAC.TableStrings.center(18, JAC.TableStrings.inUnits("energy"); na=2)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(150)) 
        #   
        for  resonance in resonances
            sa  = "";      isym = LevelSymmetry( resonance.initialLevel.J,      resonance.initialLevel.parity)
                           msym = LevelSymmetry( resonance.intermediateLevel.J, resonance.intermediateLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(resonance.initialLevel.index, resonance.intermediateLevel.index); na=4)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, msym);  na=4)
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", resonance.resonanceEnergy))          * "       "
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", resonance.augerRate))                  * "     "
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", resonance.photonRate.Coulomb))         * "  "
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", resonance.photonRate.Babushkin))       * "        "
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", resonance.resonanceStrength.Coulomb))  * "  "
            sa = sa * @sprintf("%.4e", JAC.convert("rate: from atomic", resonance.resonanceStrength.Coulomb))  * "     "
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(150))
        #
        return( nothing )
    end

end # module
