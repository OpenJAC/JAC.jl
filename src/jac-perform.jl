
"""
`Basics.perform()`  ... performs various computations.

  + `(computation::Atomic.Computation)`  
    ... to perform the computation as prescribed by comp. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

  + `(computation::Atomic.Computation; output=true)`  
    ... to perform the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the computations but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(computation::Atomic.Computation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel = computation.nuclearModel

    # Distinguish between the computation of level energies and properties and the simulation of atomic processes
    if   length(computation.configs) != 0
        ##x if  computation.asfSettings.breitCI
        ##x     @warn("perform: Test the Breit interaction; asfSettings = $(computation.asfSettings) ")
        ##x     wx = JAC.ManyElectron.Multiplet("from Ratip2012", "../test/testBreit/FeX-small-csl.inp",
        ##x                                     "../test/testBreit/FeX-small-scf.out","../test/testBreit/FeX-small-relci.mix-coulomb")  
        ##x     basis     = wx.levels[1].basis
        ##x     multiplet = perform("computation: CI", basis, nModel, computation.grid, computation.asfSettings)
        ##x     error("stop after Breit test")
        ##x end =#
        basis     = perform("computation: SCF", computation.configs, nModel, computation.grid, computation.asfSettings)
        multiplet = perform("computation: CI", basis, nModel, computation.grid, computation.asfSettings)
        if output    results = Base.merge( results, Dict("multiplet:" => multiplet) ) 
                     results = Base.merge( results, Dict("grid:"      => computation.grid) )  end
        
        # Now compute all requested properties
        if  JAC.EinsteinX      in computation.properties   
            outcome = JAC.Einstein.computeLines(multiplet,        computation.grid, computation.einsteinSettings)    
            if output    results = Base.merge( results, Dict("Einstein lines:" => outcome) )                  end
        end
        if  JAC.HFS            in computation.properties       && computation.hfsSettings.calcIJFexpansion  
            outcome = JAC.Hfs.computeHyperfineMultiplet(multiplet, nModel, computation.grid, computation.hfsSettings)         
            if output    results = Base.merge( results, Dict("IJF multiplet:" => outcome) )                   end
        end
        if  JAC.HFS            in computation.properties
            outcome = JAC.Hfs.computeOutcomes(multiplet, nModel,  computation.grid, computation.hfsSettings)         
            if output    results = Base.merge( results, Dict("HFS outcomes:" => outcome) )                    end
        end
        if  JAC.LandeJ         in computation.properties   
            outcome = JAC.LandeZeeman.computeOutcomes(multiplet, nModel,  computation.grid, computation.zeemanSettings)      
            if output    results = Base.merge( results, Dict("Zeeman parameter outcomes:" => outcome) )       end
        end
        if  JAC.Isotope        in computation.properties   
            outcome = JAC.IsotopeShift.computeOutcomes(multiplet, nModel, computation.grid, computation.isotopeSettings)         
            if output    results = Base.merge( results, Dict("Isotope parameter outcomes:" => outcome) )      end
        end
        if  JAC.AlphaX         in computation.properties   
            outcome = JAC.AlphaVariation.computeOutcomes(multiplet, nModel, computation.grid, computation.alphaSettings)         
            if output    results = Base.merge( results, Dict("alpha variation parameter outcomes:" => outcome) )      end
        end
        if  JAC.FormF          in computation.properties   
            outcome = JAC.FormFactor.computeOutcomes(multiplet, nModel, computation.grid, computation.formSettings)         
            if output    results = Base.merge( results, Dict("Form factor outcomes:" => outcome) )            end
        end
        if  JAC.Yields         in computation.properties   
            outcome = JAC.DecayYield.computeOutcomes(computation.configs, computation.asfSettings, 
                                                     multiplet, nModel, computation.grid, computation.yieldSettings)     
            if output    results = Base.merge( results, Dict("Fluorescence and AutoIonization yield outcomes:" => outcome) )   end
        end
        if  JAC.Green          in computation.properties   
            outcome = JAC.GreenFunction.computeRepresentation(multiplet, nModel, computation.grid, computation.greenSettings)         
            if output    results = Base.merge( results, Dict("Green function outcomes:" => outcome) )         end
        end
        if  JAC.Polarity       in computation.properties   
            outcome = JAC.MultipolePolarizibility.computeOutcomes(multiplet, nModel, computation.grid, computation.polaritySettings)         
            if output    results = Base.merge( results, Dict("Polarizibility outcomes:" => outcome) )         end
        end
        if  JAC.Plasma         in computation.properties   
            outcome = JAC.PlasmaShift.computeOutcomes(multiplet, nModel, computation.grid, computation.asfSettings, computation.plasmaSettings)         
            if output    results = Base.merge( results, Dict("Plasma shift outcomes:" => outcome) )            end
        end
    
    # Evaluate processes that need special SCF and CI procedures
    elseif  computation.process in [AugerInPlasma, PhotoInPlasma]
        pSettings        = computation.processSettings
        plasmaSettings   = PlasmaShift.Settings(pSettings.plasmaModel, pSettings.lambdaDebye, pSettings.ionSphereR0, pSettings.NoBoundElectrons)
        initialBasis     = perform("computation: SCF", computation.initialConfigs, nModel, computation.grid, 
                                                       computation.initialAsfSettings)
        initialMultiplet = perform("computation: CI for plasma",  initialBasis, nModel, computation.grid, computation.initialAsfSettings,
                                                                  plasmaSettings)
        finalBasis       = perform("computation: SCF", computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
        finalMultiplet   = perform("computation: CI for plasma",  finalBasis, nModel, computation.grid, computation.finalAsfSettings,
                                                                  plasmaSettings)
        #
        if      computation.process == JAC.AugerInPlasma   
            outcome = JAC.AutoIonization.computeLinesPlasma(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("AutoIonization lines in plasma:" => outcome) )         end
        elseif  computation.process == JAC.PhotoInPlasma   
            outcome = JAC.PhotoIonization.computeLinesPlasma(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Photo lines in plasma:" => outcome) )                  end
        else
            error("stop a")
        end
        
    else
        initialBasis     = perform("computation: SCF", computation.initialConfigs, nModel, computation.grid, computation.initialAsfSettings)
        initialMultiplet = perform("computation: CI",  initialBasis, nModel, computation.grid, computation.initialAsfSettings)
        finalBasis       = perform("computation: SCF", computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
        finalMultiplet   = perform("computation: CI",  finalBasis, nModel, computation.grid, computation.finalAsfSettings)
        #
        if computation.process in [PhotoExcFluor, PhotoExcAuto, PhotoIonFluor, PhotoIonAuto, ImpactExcAuto, Dierec]
            intermediateBasis     = perform("computation: SCF", computation.intermediateConfigs, nModel, computation.grid, 
                                                                computation.intermediateAsfSettings)
            intermediateMultiplet = perform("computation: CI",  intermediateBasis, nModel, computation.grid, 
                                                                computation.intermediateAsfSettings)
        end
        #
        if      computation.process == JAC.Auger   
            ##x println(" ")
            ##x for  i = 1:3   println("  perform-WARNING: Code modified to read in multiplets from Grasp computations and for testing integrals " *
            ##x                        "and Auger amplitudes !!") 
            ##x end
            ##x initialMultiplet = JAC.ManyElectron.Multiplet("from Grasp2013", "TestAuger/belike-resonance-a-csl.inp",
            ##x                                               "TestAuger/belike-resonance-a-scf.out","TestAuger/belike-resonance-a-relci.mix")  
            ##x finalMultiplet   = JAC.ManyElectron.Multiplet("from Grasp2013", "TestAuger/belike-ground-a-csl.inp",
            ##x                                               "TestAuger/belike-ground-a-scf.out","TestAuger/belike-ground-a-relci.mix") 
            outcome = JAC.AutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("AutoIonization lines:" => outcome) )                  end
        elseif  computation.process == JAC.Radiative  
            ##x println(" ")
            ##x for  i = 1:3   println("  perform-WARNING: Code modified to read in multiplets from Grasp computations and for testing integrals " *
            ##x                        "and radiative amplitudes !!") 
            ##x end
            ##x finalMultiplet   = JAC.ManyElectron.Multiplet("from Grasp2013", "TestRaditive/helike-resonance-a-csl.inp",
            ##x                                               "TestRaditive/helike-resonance-a-scf.out","TestRaditive/helike-resonance-a-relci.mix")  
            ##x initialMultiplet = JAC.ManyElectron.Multiplet("from Grasp2013", "TestRaditive/helike-resonance-a-csl.inp",
            ##x                                               "TestRaditive/helike-resonance-a-scf.out","TestRaditive/helike-resonance-a-relci.mix") 
            outcome = JAC.PhotoEmission.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("radiative lines:" => outcome) )              end
        elseif  computation.process == JAC.Photo   
            outcome = JAC.PhotoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photoionization lines:" => outcome) )        end
        elseif  computation.process == JAC.Rec   
            outcome = JAC.PhotoRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo recombination lines:" => outcome) )    end
        elseif  computation.process == JAC.Eimex   
            outcome = JAC.ImpactExcitation.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("impact-excitation lines:" => outcome) )      end
        elseif  computation.process == JAC.PhotoExc  
            outcome = JAC.PhotoExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation lines:" => outcome) )       end
        elseif  computation.process == JAC.PairA1P   
            outcome = JAC.PairAnnihilation1Photon.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("pair-annihilation-1-photon lines:" => outcome) )       end
        elseif  computation.process == JAC.MultiPhotonDE
            outcome = JAC.MultiPhotonDeExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon-de-excitation lines:" => outcome) )       end
        elseif  computation.process == JAC.RAuger
            outcome = JAC.RadiativeAuger.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("radiative Auger sharings:" => outcome) )               end
        elseif  computation.process == JAC.Coulex
            outcome = JAC.CoulombExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb excitation lines:" => outcome) )               end
        elseif  computation.process == JAC.Coulion
            outcome = JAC.CoulombIonization.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb ionization lines:" => outcome) )               end
        elseif  computation.process == JAC.PhotoExcAuto   
            outcome = JAC.PhotoExcitationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel, 
                                                                 computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation-autoionization pathways:" => outcome) )      end
        elseif  computation.process == JAC.PhotoExcFluor   
            outcome = JAC.PhotoExcitationFluores.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                 computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation-fluorescence pathways:" => outcome) )        end
        elseif  computation.process == JAC.PhotoIonAuto   
            outcome = JAC.PhotoIonizationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                 computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-ionization-autoionization pathways:" => outcome) )      end
        elseif  computation.process == JAC.PhotoIonFluor   
            outcome = JAC.PhotoIonizationFluores.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                 computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-ionizatiton-fluorescence pathways:" => outcome) )        end
        elseif  computation.process == JAC.ImpactExcAuto  
            outcome = JAC.ImpactExcitationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                  computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("impact-excitation-autoionization pathways:" => outcome) )     end
        elseif  computation.process == JAC.Compton   
            outcome = JAC.RayleighCompton.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Rayleigh-Compton lines:" => outcome) )        end
        elseif  computation.process == JAC.MultiPI   
            outcome = JAC.MultiPhotonIonization.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon single ionization:" => outcome) )        end
        elseif  computation.process == JAC.MultiPDI   
            outcome = JAC.MultiPhotonDoubleIon.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon double ionization:" => outcome) )        end
        elseif  computation.process == JAC.InternalConv   
            outcome = JAC.InternalConversion.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("internal conversion lines:" => outcome) )        end
        elseif  computation.process == JAC.Dierec  
            outcome = JAC.Dielectronic.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel, 
                                                                  computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("dielectronic recombination pathways:" => outcome) )           end
        else
            error("stop b")
        end
    end
       
    Defaults.warn(PrintWarnings)
    Defaults.warn(ResetWarnings)
    return( results )
end



"""
  + `(comp::Cascade.Computation)`  
    ... to set-up and perform a cascade computation that starts from a given set of initial configurations and proceeds via 
        various steps until a given number of electrons has been removed or the decay stops at some stable levels with regard 
        to the given atomic processes. The results of all individual steps are printed to screen but nothing is returned 
        otherwise.

  + `(comp::Cascade.Computation; output=true)`   
    ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the type 
        and specifications of the cascade but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(comp::Cascade.Computation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    # Perform the SCF and CI computation for the intial-state multiplet and print them out with their relative occupation
    basis     = perform("computation: SCF", comp.initialConfs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
    multiplet = perform("computation: CI", basis, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
    multiplet = Multiplet("initial states", multiplet.levels)
    JAC.Cascade.displayInitialLevels(stdout, multiplet, comp.initialLevels)
    if  printSummary   JAC.Cascade.displayInitialLevels(iostream, multiplet, comp.initialLevels)    end      
    if output    results = Base.merge( results, Dict("initial multiplet:" => multiplet) )           end
    #
    # Generate subsequent cascade configurations as well as display and group them together
    wa = JAC.Cascade.generateConfigurationList(comp.initialConfs, comp.maxElectronLoss, comp.NoShakeDisplacements)
    wb = JAC.Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa)
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc = JAC.Cascade.generateBlocks(comp, wb, basis.orbitals)
    JAC.Cascade.displayBlocks(stdout, wc)
    if  printSummary   JAC.Cascade.displayBlocks(iostream, wc)    end      
    # Determine, modify and compute the transition data for all steps, ie. the PhotoEmission.Line's, the AutoIonization.Line's, etc.
    wd = JAC.Cascade.determineSteps(comp, wc)
    JAC.Cascade.displaySteps(stdout, wd)
    if  printSummary   JAC.Cascade.displaySteps(iostream, wd)    end      
    we   = JAC.Cascade.modifySteps(wd)
    data = JAC.Cascade.computeSteps(comp, we)
    if output    
        results = Base.merge( results, Dict("cascade data:" => data) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        filename = "zzz-Cascade-" * string(now())[1:13] * ".jld"
        println("\nWrite all results to disk; use:\n   save($filename, results) \n   using JLD " *
                "\n   results = load($filename)    ... to load the results back from file.")
        if  printSummary   println(iostream, "\nWrite all results to disk; use:\n   save($filename, results) \n   using JLD " *
                                             "\n   results = load($filename)    ... to load the results back from file." )      end      
        save(filename, results)
    end
    ## return( results )
    return( data )
end



"""
  + `(simulation::Cascade.Simulation, data::Cascade.Data)`  
    ... to simulate a cascade decay (and excitation) by using the given data::Cascade.Data. Different computational methods and 
        different properties of the ionic system, such as the ion distribution or final-level distribution can be derived and 
        displayed from these simulations. Of course, the details of these simulations strongly depend on the atomic processes 
        and data that have been generated before by performing a computation::Cascade.Computation. The results of all individual 
        steps are printed to screen but nothing is returned otherwise.

  + `(simulation::Cascade.Simulation, data::Cascade.Data; output=true)`   
    ... to perform the same but to return the complete output in a dictionary; the particular output depends on the method and 
        specifications of the cascade but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(simulation::Cascade.Simulation, data::Cascade.Data; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    #
    # Distinguish between the different computational methods for running the simulations
    if  simulation.method == Cascade.ProbPropagation()
        if   JAC.Cascade.FinalDist()  in  simulation.properties   ||   JAC.Cascade.IonDist()  in  simulation.properties
            wa = Cascade.simulateLevelDistribution(simulation, data) 
        end
        if   JAC.Cascade.ElectronIntensity()    in simulation.properties    ||   JAC.Cascade.PhotonIntensity()  in simulation.properties   ||
             JAC.Cascade.ElectronCoincidence()  in simulation.properties 
            error("stop a: Not yet implemented")    
        end
    else  error("stop b")
    end

    if output 
        wx = 0.   
        results = Base.merge( results, Dict("ion distribution:" => wx) )
        #  Write out the result to file to later continue with simulations on the cascade data
        filename = "xxx-simulation-" * string(now())[1:13] * ".jld"
        println("\nWrite all results to disk; use:\n   save($filename, results) \n   using JLD " *
                "\n   results = load($filename)    ... to load the results back from file.")
        save(filename, results)
    end
    return( results )
end



"""
  + `(computation::Atomic.CasComputation)`  
    ... to perform the computation ... . Nothing is returned.  **Not yet implemented !**

  + `(computation::Atomic.CasComputation; output=true)`  
    ... to perform the same ....  **Not yet implemented !**
"""
function Basics.perform(computation::Atomic.CasComputation; output::Bool=false)
    error("Not yet implemented")
end



"""
  + `("computation: SCF", configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings;
                          printout::Bool=true)`  
    ... to generate an atomic basis and to compute the self-consistent field (SCF) for this basis due to the given settings; 
        a basis::Basis is returned.  
"""
function Basics.perform(sa::String, configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings;
                 printout::Bool=true)
    !(sa == "computation: SCF")   &&   error("Unsupported keystring = $sa")
    if  printout    println("\n... in perform('computation: SCF', ...")    end
    
    # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generate("configuration list: relativistic", conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println("perform-aa: ", relconfList[i])    end   end
    subshellList = Basics.generate("subshells: ordered list for relativistic configurations", relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generate("CSF list: from single ConfigurationR", relconf, subshellList)
        append!( csfList, newCsfs)
    end

    # Determine the number of electrons and the list of coreSubshells
    NoElectrons      = sum( csfList[1].occupation )
    coreSubshellList = Subshell[]
    for  k in 1:length(subshellList)
        mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    is_filled = false;    break   end
        end
        if   is_filled    push!( coreSubshellList, subshellList[k])    end
    end
        
    ##x waL = JAC.Bsplines.generatePrimitives(7, 7, grid);    nsL = length(waL.bsplines) - 3
    ##x waS = JAC.Bsplines.generatePrimitives(8, 7, grid);    nsS = length(waS.bsplines) - 3 
    wa = Bsplines.generatePrimitives(grid)

    # Generate start orbitals
    if  settings.startScf == "hydrogenic"
        # Generate start orbitals for the SCF field by using B-splines
        ##x orbitals  = JAC.Bsplines.generateOrbitalsHydrogenic(waL, nsL, waS, nsS, nuclearModel, subshellList)
        orbitals  = JAC.Bsplines.generateOrbitalsHydrogenic(wa, nuclearModel, subshellList; printout=printout)
    elseif  settings.startScf == "fromNRorbitals"
        # Generate starting orbitals for this csfList by adapting non-relativistic orbitals with a proper nuclear charge
        orbitals = Dict{Subshell, Orbital}()
        for  subsh in subshellList
            orb      = JAC.HydrogenicIon.radialOrbital(subsh, nuclearModel.Z, grid)
            orbitals = Base.merge( orbitals, Dict( subsh => orb) ) 
        end
    else  error("stop a")
    end

    ##x error("stop ** for test reasons ** ")
    
    # Perform the SCF computations for the orbitals due to the given settings
    basis    = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, orbitals)
    if   settings.methodScf in ["meanDFS", "meanHS"]
        ##x newBasis = JAC.Bsplines.solveSelfConsistentFieldMean(waL, nsL, waS, nsS, nuclearModel, basis, settings) 
        newBasis = JAC.Bsplines.solveSelfConsistentFieldMean(wa, nuclearModel, basis, settings; printout=printout) 
    elseif   settings.methodScf in ["AL", "OL"] 
        ##x newBasis = JAC.Bsplines.solveSelfConsistentFieldOptimized(waL, nsL, waS, nsS, nuclearModel, basis, settings)
        newBasis = JAC.Bsplines.solveSelfConsistentFieldOptimized(wa, nuclearModel, basis, settings; printout=printout)
    else  error("stop b")
    end
    
    return( newBasis )  
end



"""
  + `("computation: mutiplet from orbitals, no CI, CSF diagonal", configs::Array{Configuration,1}, 
        initalOrbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)` 
    ... to generate from the given initial orbitals a multiplet of single-CSF levels by just using the diagonal 
        part of the Hamiltonian matrix; a multiplet::Multiplet is returned.  
"""
function Basics.perform(sa::String, configs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, 
                 grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)
    !(sa == "computation: mutiplet from orbitals, no CI, CSF diagonal")   &&   error("Unsupported keystring = $sa")
    if  printout    println("\n... in perform('computation: mutiplet from orbitals, no CI, CSF diagonal', ...")    end
    
    # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generate("configuration list: relativistic", conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println("perform-aa: ", relconfList[i])    end   end
    subshellList = Basics.generate("subshells: ordered list for relativistic configurations", relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generate("CSF list: from single ConfigurationR", relconf, subshellList)
        append!( csfList, newCsfs)
    end

    # Determine the number of electrons and the list of coreSubshells
    NoElectrons      = sum( csfList[1].occupation )
    coreSubshellList = Subshell[]
    for  k in 1:length(subshellList)
        mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    is_filled = false;    break   end
        end
        if   is_filled    push!( coreSubshellList, subshellList[k])    end
    end
    
    # Set-up a basis for calculating the Hamiltonian matrix
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, initalOrbitals)
   
    # Calculate the diagonal part of the Hamiltonian matrix and define a multiplet for these approximate ASf;
    # Generate first an effective nuclear charge Z(r) on the given grid
    potential = JAC.Nuclear.nuclearPotential(nuclearModel, grid)
    
    n = length(csfList);    matrix = zeros(Float64, n, n)
    for  r = 1:n
        wa = compute("angular coefficients: e-e, Ratip2013", basis.csfs[r], basis.csfs[r])
        me = 0.
        for  coeff in wa[1]
            jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
            me = me + coeff.T * sqrt( jj + 1) * JAC.RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
        end

        for  coeff in wa[2]
            if  settings.coulombCI    
                me = me + coeff.V * JAC.InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                 basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)   end
            if  settings.breitCI
                me = me + coeff.V * JAC.InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                               basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)     end
        end
        matrix[r,r] = me
    end 
    #
    eigen  = Basics.diagonalize("matrix: Julia, eigfact", matrix)
    levels = Level[]
    for  ev = 1:length(eigen.values)
        vector = eigen.vectors[ev];   ns = 0
        for  r = 1:length(vector)   
            if      (vector[r])^2 > 0.99    ns = r;   break    
            elseif  (vector[r])^2 > 0.01    error("stop a; unexpected mixing coefficieint; vector = $vector ")    
            end
        end
        J = csfList[ns].J
        newlevel = Level( J, AngularM64(J.num//J.den), csfList[ns].parity, 0, eigen.values[ev], 0., true, basis, vector ) 
        push!( levels, newlevel)
    end
    mp = Multiplet("noName", levels)
    mp = Basics.sort("multiplet: by energy", mp)
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate("multiplet: energies", mp)
        Basics.tabulate("multiplet: energy relative to immediately lower level",    mp)
        Basics.tabulate("multiplet: energy of each level relative to lowest level", mp)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary     
        Basics.tabulate("multiplet: energies", mp, stream=iostream)
        Basics.tabulate("multiplet: energy relative to immediately lower level",    mp, stream=iostream)
        Basics.tabulate("multiplet: energy of each level relative to lowest level", mp, stream=iostream)
    end
  
    return( mp )
end


"""
  + `("computation: CI", basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)`  
    ... to  set-up and diagonalize from the (SCF) basis the configuration-interaction matrix and to derive and display the 
        level structure of the corresponding multiplet due to the given settings; a multiplet::Multiplet is returned.   
"""
function Basics.perform(sa::String, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)
    !(sa == "computation: CI")   &&   error("Unsupported keystring = $sa")
    
    # Determine the J^P symmetry blocks
    symmetries = Dict{JAC.LevelSymmetry,Int64}()
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if     haskey(symmetries, sym)    symmetries[sym] = symmetries[sym] + 1
        else   symmetries = Base.merge( symmetries, Dict( sym => 1 ) )
        end
    end

    # Test the total number of CSF
    NoCsf = 0
    for (sym,v) in symmetries   NoCsf = NoCsf + v   end
   
    # Calculate for each symmetry block the corresponding CI matrix, diagonalize it and append a Multiplet for this block
    multiplets = Multiplet[]
    for  (sym,v) in  symmetries
        matrix = compute("matrix: CI, J^P symmetry", sym, basis, nuclearModel, grid, settings; printout=printout)
        eigen  = Basics.diagonalize("matrix: Julia, eigfact", matrix)
        levels = Level[]
        for  ev = 1:length(eigen.values)
            # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
            evSym = eigen.vectors[ev];    vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs) 
                if LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = evSym[ns]   end
            end
            newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector ) 
            push!( levels, newlevel)
        end
        wa = Multiplet("noName", levels)
        push!( multiplets, wa)
    end
    
    # Merge all multiplets into a single one
    mp = Basics.merge("multiplets", multiplets)
    mp = Basics.sort("multiplet: by energy", mp)
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate("multiplet: energies", mp)
        Basics.tabulate("multiplet: energy relative to immediately lower level",    mp)
        Basics.tabulate("multiplet: energy of each level relative to lowest level", mp)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary     
        Basics.tabulate("multiplet: energies", mp, stream=iostream)
        Basics.tabulate("multiplet: energy relative to immediately lower level",    mp, stream=iostream)
        Basics.tabulate("multiplet: energy of each level relative to lowest level", mp, stream=iostream)
    end
  
    return( mp )
end


"""
  + `("computation: CI for plasma", basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, 
    settings::AsfSettings, plasmaSettings::PlasmaShift.Settings; printout::Bool=true)`  
    ... to  set-up and diagonalize from the given (SCF) basis the configuration-interaction matrix and to derive and
        display the level structure of the corresponding multiplet due to the given settings. Here, the CI matrix
        includes the modifications of the Hamiltonian due to the given plasmaSettings; a multiplet::Multiplet is returned.   
"""
function Basics.perform(sa::String, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, 
                 settings::AsfSettings, plasmaSettings::PlasmaShift.Settings; printout::Bool=true)
    !(sa == "computation: CI for plasma")   &&   error("Unsupported keystring = $sa")
    
    # Determine the J^P symmetry blocks
    symmetries = Dict{JAC.LevelSymmetry,Int64}()
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if     haskey(symmetries, sym)    symmetries[sym] = symmetries[sym] + 1
        else   symmetries = Base.merge( symmetries, Dict( sym => 1 ) )
        end
    end

    # Test the total number of CSF
    NoCsf = 0
    for (sym,v) in symmetries   NoCsf = NoCsf + v   end
   
    # Calculate for each symmetry block the corresponding CI matrix for the given plasma model and parameters as specified by
    # plasmaSettings, diagonalize it and append a Multiplet for this block
    multiplets = Multiplet[]
    for  (sym,v) in  symmetries
        matrix = compute("matrix: CI for plasma, J^P symmetry", sym, basis, nuclearModel, grid, settings, plasmaSettings; printout=printout)
        eigen  = Basics.diagonalize("matrix: Julia, eigfact", matrix)
        levels = Level[]
        for  ev = 1:length(eigen.values)
            # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
            evSym = eigen.vectors[ev];    vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs) 
                if LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = evSym[ns]   end
            end
            newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector ) 
            push!( levels, newlevel)
        end
        wa = Multiplet("noName", levels)
        push!( multiplets, wa)
    end
    
    # Merge all multiplets into a single one
    mp = Basics.merge("multiplets", multiplets)
    mp = Basics.sort("multiplet: by energy", mp)
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate("multiplet: energies", mp)
        Basics.tabulate("multiplet: energy relative to immediately lower level",    mp)
        Basics.tabulate("multiplet: energy of each level relative to lowest level", mp)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary     
        Basics.tabulate("multiplet: energies", mp, stream=iostream)
        Basics.tabulate("multiplet: energy relative to immediately lower level",    mp, stream=iostream)
        Basics.tabulate("multiplet: energy of each level relative to lowest level", mp, stream=iostream)
    end
  
    return( mp )
end

