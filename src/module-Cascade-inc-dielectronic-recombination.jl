
# Functions and methods for cheme::Cascade.DielectronicRecombinationScheme computations


"""
`Cascade.computeSteps(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested capture & transition amplitudes as well as DielectronicCapture.Line's, AutoIonization.Line's, 
        etc. for all pre-specified decay steps of the cascade. When compared with standard computations of these atomic 
        processes, however, the amount of output is largely reduced and often just printed into the summary file. 
        A set of  data::Cascade.CaptureData  is returned.
"""
function computeSteps(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[];    cOrbitals = Dict{Subshell, Orbital}()
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nt = st = 0;   previousMeanEn = 0.
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        sa = "\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc decay lines (without selection rules): "
        println(sa);    if  printSummary   println(iostream, sa)   end 
                                                
        if      step.process == Basics.Auger()  &&   comp.approach == Cascade.AverageSCA()
            meanEn = 0.;    NoEn = 0
            for  p = 1:length(step.initialMultiplet.levels),  q = 1:length(step.finalMultiplet.levels)
                en = step.initialMultiplet.levels[p].energy - step.finalMultiplet.levels[q].energy
                if  en > 0.03    meanEn = meanEn + en;    NoEn = NoEn + 1    end
            end
            if  NoEn > 0     meanEn = meanEn/ NoEn     
            else     
                println(">> No transition with positive energy.")
                continue    
            end
        
            # Generate potential for continuum orbitals for this step
            nrContinuum  = Continuum.gridConsistency(meanEn + 0.1, comp.grid)
            contSettings = Continuum.Settings(false, nrContinuum);   
            npot         = Nuclear.nuclearPotential(comp.nuclearModel, comp.grid)
            ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
            ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
            wp           = Basics.computePotential(comp.asfSettings.scField, comp.grid, step.finalMultiplet.levels[1])         
            pot          = Basics.add(npot, wp)
            #  Generate continuum if not yet available
            generatedKappas = Int64[]
            for  p = 1:length(step.initialMultiplet.levels),  q = 1:length(step.finalMultiplet.levels)
                en = step.initialMultiplet.levels[p].energy - step.finalMultiplet.levels[q].energy
                if  en > 0.03    
                    symi      = LevelSymmetry(step.initialMultiplet.levels[p].J, step.initialMultiplet.levels[p].parity)
                    symf      = LevelSymmetry(step.finalMultiplet.levels[q].J, step.finalMultiplet.levels[q].parity)
                    kappaList = AngularMomentum.allowedKappaSymmetries(symi, symf)
                    for  kappa in kappaList
                        sh    = Subshell(101, kappa)
                        if      abs(kappa) > step.settings.maxKappa                               ## continuum orbital not needed
                        elseif  kappa in generatedKappas                                          ## already generated for this step
                        elseif  haskey(cOrbitals,sh)   &&   cOrbitals[sh].energy - en / en < 0.15 ## use previous one
                            println(">> No new continum orbital generated for $sh and energy $en ")
                        else                                                                      ## generate new continuum orbital
                            push!(generatedKappas, kappa)
                            cOrbital, phase, normF  = Continuum.generateOrbitalLocalPotential(en, sh, pot, contSettings)
                            cOrbitals[sh]           = cOrbital
                            println(">> New continum orbital generated for $sh and energy $en ")
                        end
                    end
                end
            end
        
            newLines = AutoIonization.computeLinesFromOrbitals(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                                step.settings, cOrbitals, output=true, printout=false) 
            append!(linesA, newLines);    nt = length(linesA)
        elseif  step.process == Basics.Auger() 
            # Compute continuum orbitals independently for all transitions in the given block.
            newLines = AutoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                            step.settings, output=true, printout=false) 
            append!(linesA, newLines);    nt = length(linesA)
        elseif  step.process == Basics.Radiative()
            newLines = PhotoEmission.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.grid, 
                                                            step.settings, output=true, printout=false) 
            append!(linesR, newLines);    nt = length(linesR)
        else   error("Unsupported atomic process for cascade computations.")
        end
        sa = "     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) decay lines."
        println(sa);    if  printSummary   println(iostream, sa)   end 
    end
    #
    data = [ Cascade.Data{PhotoEmission.Line}(linesR), Cascade.Data{AutoIonization.Line}(linesA) ]
end

"""
`Cascade.determineSteps(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, captureList::Array{Cascade.Block,1}, decayList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this dielectronic cascade. It considers the autoionization
        between the blocks from the captureList and initialList as well as the radiative stabilization between the blocks from the 
        captureList and decayList. It checks that at least on pair of levels supports either an `electron-capture' or 
        `radiative stabilization' within the step. A stepList::Array{Cascade.Step,1} is returned, and for which subsequently all 
        required transition amplitudes and rates/cross sections are computed.
"""
function determineSteps(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, capturedList::Array{Cascade.Block,1}, decayList::Array{Cascade.Block,1})
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        for  capturedBlock in capturedList
            for  initialBlock in initialList
                if  initialBlock.NoElectrons + 1 != capturedBlock.NoElectrons   error("stop a")     end
                # Check that at least one energy supports autoionization
                settings = AutoIonization.Settings(AutoIonization.Settings(), maxKappa=7)
                push!( stepList, Cascade.Step(Basics.Auger(), settings, capturedBlock.confs,     initialBlock.confs,
                                                                        capturedBlock.multiplet, initialBlock.multiplet) )
            end
        end
        #
        for  capturedBlock in capturedList
            for  decayBlock in decayList
                if  decayBlock.NoElectrons != capturedBlock.NoElectrons   error("stop b")     end
                # Check that at least one energy supports radiative stabilization
                if   Basics.determineMeanEnergy(capturedBlock.multiplet) - Basics.determineMeanEnergy(decayBlock.multiplet) > scheme.minPhotonEnergy
                    settings = PhotoEmission.Settings(PhotoEmission.Settings(), gauges=[UseCoulomb, UseBabushkin])
                    push!( stepList, Cascade.Step(Basics.Radiative(), settings, capturedBlock.confs,     decayBlock.confs,
                                                                                capturedBlock.multiplet, decayBlock.multiplet) )
                end
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


"""
`Cascade.generateBlocks(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this electron-capture and subsequent stabilization (DR) cascade, 
        and compute also the corresponding multiplets. The different cascade approches realized different strategies how these blocks are 
        selected and computed. A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[];    basis = Basis()
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        sa = "\n* Generate blocks for DR plasma rate coefficient computations: \n" *
                "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: " *
                "\n    + orbitals are generated independently for the first block in a Dirac-Fock-Slater potential; " *
                "\n    + these orbitals are re-used for all other block, together with hydrogenic orbitals for the outer part; " *
                "\n    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; " *
                "\n    + only the Coulomb interaction is considered for the electron capture. " *
                "\n    + only E1 excitations are considered for the stabilization. \n"
        if  printout       println(sa)              end
        if  printSummary   println(iostream, sa)    end
        #
        if  length(confs) > 1
            # Determine a list of hydrogenic orbitals for later use 
            relconfList = ConfigurationR[]
            for  confa in confs
                wa = Basics.generateConfigurationRs(confa)
                append!( relconfList, wa)
            end
            subshellList = Basics.generateSubshellList(relconfList)
            Defaults.setDefaults("relativistic subshell list", subshellList; printout=printout)
            wa                 = BsplinesN.generatePrimitives(comp.grid)
            ##x hydrogenicOrbitals = BsplinesN.generateOrbitalsHydrogenic(wa, comp.nuclearModel, subshellList; printout=printout)
            hydrogenicOrbitals = BsplinesN.generateOrbitalsHydrogenic(subshellList, comp.nuclearModel, wa; printout=printout)
        end
        
        for  (ia, confa)  in  enumerate(confs)
            sa = "  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... "
            print(sa);      if  printSummary   println(iostream, sa)   end
            # Now distinguish between the first and all other blocks; for the first block, a SCF is generated and the occupied orbital
            # used also for all other blocks. In addition, a set of hydrogenic orbitals generated for later use
            if  ia == 1
                ##x basis     = Basics.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                multiplet     = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                basis         = multiplet.levels[1].basis
            else
                # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
                relconfList  = ConfigurationR[]
                wa           = Basics.generateConfigurationRs(confa);    append!( relconfList, wa)
                subshellList = Basics.generateSubshellList(relconfList)
                Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
                # Generate the relativistic CSF's for the given subshell list
                csfList = CsfR[]
                for  relconf in relconfList
                    newCsfs = Basics.generateCsfRs(relconf, subshellList)
                    append!( csfList, newCsfs)
                end
                # Determine the list of coreSubshells
                coreSubshellList = Subshell[]
                for  k in 1:length(subshellList)
                    mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
                    for  csf in csfList
                        if  csf.occupation[k] != mocc    is_filled = false;    break   end
                    end
                    if   is_filled    push!( coreSubshellList, subshellList[k])    end
                end
                # Add all missing orbitals as hydrogenic
                orbitals      = copy(basis.orbitals)
                for  subsh in subshellList
                    if haskey(orbitals, subsh)   ## do nothing
                    else      orbitals[subsh] = hydrogenicOrbitals[subsh];   print("hydrogenic $subsh ...")
                    end
                end
                
                basis         = Basis(true, confa.NoElectrons, subshellList, csfList, coreSubshellList, orbitals)
            end
            ##x multiplet = Basics.perform("computation: mutiplet from orbitals, no CI, CSF diagonal", [confa],  basis.orbitals, 
            ##x                             comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            multiplet = Hamiltonian.performCIwithFrozenOrbitals([confa],  basis.orbitals, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    elseif    comp.approach == SCA()
        sa = "\n* Generate blocks for electron-capture & decay cascade computations: \n" *
                "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: " *
                "\n    + each single configuration forms an individual cascade block; " *
                "\n    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; " *
                "\n    + configuration mixing is included for each block, based on H^(DC); " *
                "\n    + all requested multipoles are considered for the stabilization. \n"
        if  printout       println(sa)              end
        if  printSummary   println(iostream, sa)    end
        #
        for  confa  in confs
            print("  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")
            if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")   end
                ##x basis     = Basics.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                ##x multiplet = Basics.performCI(basis,    comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                multiplet = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    else  error("Unsupported cascade approach.")
    end

    return( blockList )
end


"""
`Cascade.generateConfigurationsForDielectronicCapture(multiplets::Array{Multiplet,1},  scheme::DielectronicRecombinationScheme, 
                                                    nm::Nuclear.Model, grid::Radial.Grid)`  
    ... generates all possible doubly-excited configurations due to (dielectronic) electron capture into the given multiplets.
        The number and type of such doubly-generated configurations depend on (1) the maximum (electron) energy for capturing an electron
        that is closely related to the (maximum) temperature of the plasma; (2) the fromShells from which (and how many displacements)
        are accepted as well as (3) the maximum principle and orbital angular quantum number of the additional (to-) shellsthe fromShells
        into which electrons excited and/or captured. A Tuple(initialConfList::Array{Configuration,1}, confList::Array{Configuration,1}) 
        is returned.
"""
function generateConfigurationsForDielectronicCapture(multiplets::Array{Multiplet,1},  scheme::DielectronicRecombinationScheme, 
                                                    nm::Nuclear.Model, grid::Radial.Grid)
    # Determine all (reference) configurations from multiplets and generate the 'excited' configurations due to the specificed excitations
    initialConfList = Configuration[]
    for mp  in  multiplets   
        confList = Basics.extractNonrelativisticConfigurations(mp.levels[1].basis)
        for  conf in confList   if  conf in initialConfList   nothing   else   push!(initialConfList, conf)      end      end
    end
    coreConfList    = Basics.generateConfigurations(initialConfList, scheme.excitationFromShells, scheme.excitationToShells, 
                                                    scheme.NoExcitations)
    captureConfList = Basics.generateConfigurationsWithElectronCapture(coreConfList, scheme.excitationFromShells, scheme.intoShells, 0)
    decayConfList   = Basics.generateConfigurationsWithElectronCapture(coreConfList, scheme.excitationFromShells, scheme.decayShells, 0)
    @show initialConfList
    @show coreConfList
    @show captureConfList
    @show decayConfList
    #
    # Determine first a hydrogenic spectrum for all subshells of the initial and doubly-excited states
    allConfList   = Configuration[];      
    append!(allConfList, initialConfList);      append!(allConfList, captureConfList);      append!(allConfList, decayConfList)
    
    allSubshells  = Basics.extractRelativisticSubshellList(allConfList)
    primitives    = BsplinesN.generatePrimitives(grid)
    orbitals      = BsplinesN.generateOrbitalsHydrogenic(allSubshells, nm, primitives, printout=true)
    # Exclude configurations with too high mean energies
    en            = Float64[];   
    for conf in initialConfList
        wen = Basics.determineMeanEnergy(conf, orbitals, nm, grid)
        push!(en, wen)
    end
    initialMean = sum(en) / length(en)
    println(">>> initial configuration(s) have mena energies  $initialMean  [a.u.].")
    #
    newCaptureConfList = Configuration[];  meanEnergies = Float64[]
    for  conf  in  captureConfList    
        confMean = Basics.determineMeanEnergy(conf, orbitals, nm, grid)
        if  confMean - initialMean <= scheme.maxExcitationEnergy     push!(newCaptureConfList, conf)     end
    end

    return( (initialConfList, newCaptureConfList, decayConfList)  )
end


"""
`Cascade.perform(scheme::DielectronicRecombinationScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a dielectronic-recombination (DR) plasma rate coefficient computation that combines the electron capture
        and the subsequent radiative stabilization steps. Such a computation starts from a given set of initial configurations xor 
        initial multiplets and (1) generates all doubly-excited configurations due to the capture of an electron with a given maximum
        electron energy; (2) selects all electron capture (inverse Auger) and re-autoionization (Auger) steps and (3) selects
        all steps for radiative stabilization due to given parameters of the scheme::DielectronicRecombinationScheme. The results of 
        these DR plasma rate computation are comprised into (output) data::ExcitationData, while these data are only printed during 
        the generation and nothing is returned.

`Cascade.perform(scheme::DielectronicRecombinationScheme, comp::Cascade.Computation; output=true, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::DielectronicRecombinationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    # Perform the SCF and CI computation for the intial-state multiplets if initial configurations are given
    if  comp.initialConfigs != Configuration[]
        ##x basis      = Basics.performSCF(comp.initialConfigs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
        ##x multiplet  = Basics.performCI(basis, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
        multiplet  = SelfConsistent.performSCF(comp.initialConfigs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
        multiplets = [Multiplet("initial states", multiplet.levels)]
    else
        #== # Shift the initial level energy by -electronEnergyShift
        if  scheme.electronEnergyShift != 0.  
            multiplets = Multiplet
            for  multiplet  in  comp.initialMultiplets
                newMultiplet = Basics.shiftTotalEnergies(multiplet, Defaults.convertUnits("energy: to atomic", -scheme.electronEnergyShift)) 
                push!(multiplets, newMultiplet)
            end
            println(">> Shift all initial level energies by -$(scheme.electronEnergyShift) $(Defaults.getDefaults("unit: energy"))")
        else
            multiplets = comp.initialMultiplets
        end  ==#
        multiplets = comp.initialMultiplets
    end
    # Print out initial configurations and levels 
    Cascade.displayLevels(stdout, multiplets, sa="initial ")
    if  printSummary   Cascade.displayLevels(iostream, multiplets, sa="initial ")                            end
    #
    # Generate subsequent cascade configurations as well as display and group them together
    wa  = Cascade.generateConfigurationsForDielectronicCapture(multiplets, comp.scheme, comp.nuclearModel, comp.grid)
    wb1 = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa[1], sa="initial configurations of the DR cascade ")
    wb2 = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa[2], sa="doubly-excited capture configurations of the DR cascade ")
    wb3 = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa[3], sa="decay configurations of the DR cascade ")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc1 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb1)
    wc2 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb2, printout=false)
    wc3 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb3, printout=false)
    # Shift the initial level energy by -electronEnergyShift
    @show scheme.electronEnergyShift
    if  scheme.electronEnergyShift != 0. 
        wc1old = wc1;   wc1 = Cascade.Block[]
        for  block in wc1old
            newMultiplet = Basics.shiftTotalEnergies(multiplet, Defaults.convertUnits("energy: to atomic", -scheme.electronEnergyShift))
            push!(wc1, Cascade.Block(block.NoElectrons, block.confs, block.hasMultiplet, newMultiplet))
        end
        println(">> Shift all initial level energies by $(-scheme.electronEnergyShift) $(Defaults.getDefaults("unit: energy"))")
    end
    #
    Cascade.displayBlocks(stdout, wc1, sa="from the initial configurations of the DR cascade ");      
    Cascade.displayBlocks(stdout, wc2, sa="from the doubly-excited capture configurations of the DR cascade ")
    Cascade.displayBlocks(stdout, wc3, sa="from the decay configurations of the DR cascade ")
    if  printSummary   Cascade.displayBlocks(iostream, wc1, sa="from the initial configurations of the DR cascade ");      
                        Cascade.displayBlocks(iostream, wc2, sa="from the doubly-excited capture configurations of the DR cascade ")
                        Cascade.displayBlocks(iostream, wc3, sa="from the decay configurations of the DR cascade ")            end      
    #
    # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
    gMultiplets = Multiplet[];     
    for block in wc1  push!(gMultiplets, block.multiplet)    end
    for block in wc2  push!(gMultiplets, block.multiplet)    end
    for block in wc3  push!(gMultiplets, block.multiplet)    end
    we = Cascade.determineSteps(scheme, comp, wc1, wc2, wc3)
    Cascade.displaySteps(stdout, we, sa="electron capture and stabilization ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="electron capture and stabilization ")    end      
    wf   = Cascade.modifySteps(we)
    #
    data = Cascade.computeSteps(scheme, comp, wf)
    if output    
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => multiplets) )    
        results = Base.merge( results, Dict("generated multiplets:"         => gMultiplets) )    
        results = Base.merge( results, Dict("cascade data:"                 => data) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        filename = "zzz-cascade-dr-rate-computations-" * string(Dates.now())[1:13] * ".jld"
        println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
        if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
        JLD2.@save filename results
    end
    ## return( results )
    return( results )
end
