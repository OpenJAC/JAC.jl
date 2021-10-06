
    # Functions and methods for cascade computation


    """
    `Cascade.computeSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
        ... computes in turn all the requested transition amplitudes as well as PhotoEmission.Line's, AutoIonization.Line's, 
            etc. for all pre-specified decay steps of the cascade. When compared with standard computations of these atomic 
            processes, however, the amount of output is largely reduced and often just printed into the summary file. 
            A set of  data::Cascade.DecayData  is returned.
    """
    function computeSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
        linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[];    cOrbitals = Dict{Subshell, Orbital}()    
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        nt = 0;   st = 0;   previousMeanEn = 0.
        for  step  in  stepList
            st = st + 1
            nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
            println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc decay lines (without selection rules): ")
            if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                 "up to $nc decay lines (without selection rules): ")   end 
                                                 
            if      step.process == Basics.Auger()  &&   comp.approach == Cascade.AverageSCA()
                # First determine the `mean' free-electron energy for this Auger block and calculate a common set of continuum orbital
                meanEn = 0.;    NoEn = 0
                for  p = 1:length(step.initialMultiplet.levels),  q = 1:length(step.finalMultiplet.levels)
                    en = step.initialMultiplet.levels[p].energy - step.finalMultiplet.levels[q].energy
                    if  en > 0.1    meanEn = meanEn + en;    NoEn = NoEn + 1    end
                end
                if  NoEn > 0     meanEn = meanEn/ NoEn     else     meanEn = 0.1    end
                if  abs(meanEn - previousMeanEn) / meanEn  < 0.15   # no new continuum orbitals
                    println(">> No new continum orbitals are generated for $(keys(cOrbitals)) and for the energy $meanEn ")
                else
                    previousMeanEn = meanEn;    cOrbitals = Dict{Subshell, Orbital}()
                    for kappa = -step.settings.maxKappa-1:step.settings.maxKappa        if   kappa == 0     continue    end
                        sh           = Subshell(101, kappa);     nrContinuum = Continuum.gridConsistency(meanEn, comp.grid)
                        contSettings = Continuum.Settings(false, nrContinuum);   
                        npot         = Nuclear.nuclearPotential(comp.nuclearModel, comp.grid)
                        ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
                        ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
                        ## wp3 = compute("radial potential: Kohn-Sham", grid, wLevel)
                        ## wp           = Basics.compute("radial potential: Dirac-Fock-Slater", comp.grid, step.finalMultiplet.levels[1].basis)
                        wp           = Basics.computePotentialDFS(comp.grid, step.finalMultiplet.levels[1])
                        pot          = Basics.add(npot, wp)
                        cOrbital, phase, normF  = Continuum.generateOrbitalLocalPotential(meanEn, sh, pot, contSettings)
                        cOrbitals[sh] = cOrbital
                    end
                    println(">> Generate continum orbitals in DFS potential for $(keys(cOrbitals)) and for the energy $meanEn ")
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
            println("     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                    "to a total of $nt $(string(step.process)) decay lines." )
            if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                 "giving now rise to a total of $nt $(string(step.process)) decay lines." )   end      
        end
        #
        data = Cascade.DecayData(linesR, linesA)
    end


    """
    `Cascade.determineSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, blockList::Array{Cascade.Block,1})`  
        ... determines all step::Cascade.Step's that need to be computed for this decay cascade. It cycles through all processes of the given
            decay scheme and selects all pairs of blocks due to the selected cascade approach. It is checked that the (averaged) energies 
            each block or level supports a `decay' within the step. A stepList::Array{Cascade.Step,1} is returned, and for which subsequently
            all required transition amplitudes and rates are computed.
    """
    function determineSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, blockList::Array{Cascade.Block,1})
        stepList = Cascade.Step[]
        if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
            for  a = 1:length(blockList)
                for  b = 1:length(blockList)
                    minEn = 100000.;   maxEn = -100000.
                    for  p = 1:length(blockList[a].multiplet.levels),  q = 1:length(blockList[b].multiplet.levels)
                        minEn = min(minEn, blockList[a].multiplet.levels[p].energy - blockList[b].multiplet.levels[q].energy)
                        maxEn = max(maxEn, blockList[a].multiplet.levels[p].energy - blockList[b].multiplet.levels[q].energy)
                    end
                    for  process  in  comp.scheme.processes
                        if      process == Basics.Radiative()   
                            if  a == b   ||   maxEn < 0.    continue   end
                            if  blockList[a].NoElectrons == blockList[b].NoElectrons
                                settings = PhotoEmission.Settings([E1], [UseBabushkin], false, false, LineSelection(), 0., 0., 1.0e6)
                                push!( stepList, Cascade.Step(process, settings, blockList[a].confs, blockList[b].confs, 
                                                              blockList[a].multiplet, blockList[b].multiplet) )
                            end
                        elseif  process == Basics.Auger()       
                            if  a == b   ||   maxEn < 0.    continue   end
                            if  blockList[a].NoElectrons == blockList[b].NoElectrons + 1
                                settings = AutoIonization.Settings(false, false, LineSelection(), 0., 1.0e6, 3, CoulombInteraction())
                                push!( stepList, Cascade.Step(process, settings, blockList[a].confs, blockList[b].confs, 
                                                              blockList[a].multiplet, blockList[b].multiplet) )
                            end
                        else    error("stop a")
                        end
                    end
                end
            end
            #
        else  error("Unsupported cascade approach.")
        end
        return( stepList )
    end


    """
    `Cascade.perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation)`  
        ... to set-up and perform a cascade computation that starts from a given set of initial configurations and proceeds via 
            various steps until a given number of electrons has been removed or the decay stops at some stable levels with regard 
            to the given atomic processes. The results of all individual steps are printed to screen but nothing is returned 
            otherwise.

    `Cascade.perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)`   
        ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the type 
            and specifications of the cascade but can easily accessed by the keys of this dictionary.
    """
    function perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)
        if  output    results = Dict{String, Any}()    else    results = nothing    end
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        # Perform the SCF and CI computation for the intial-state multiplets if initial configurations are given
        if  comp.initialConfigs != Configuration[]
            basis      = Basics.performSCF(comp.initialConfigs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            multiplet  = Basics.performCI(basis, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            multiplets = [Multiplet("initial states", multiplet.levels)]
        else
            basis      = comp.initialMultiplets[1].levels[1].basis
            multiplets = comp.initialMultiplets
        end
        # Print out initial configurations and levels 
        Cascade.displayLevels(stdout, multiplets, sa="initial ")
        if  printSummary   Cascade.displayLevels(iostream, multiplets, sa="initial ")          end
        #
        # Generate subsequent cascade configurations as well as display and group them together
        wa = Cascade.generateConfigurationList(multiplets, comp.scheme.maxElectronLoss, comp.scheme.NoShakeDisplacements)
        wb = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa, sa="decay ")
        #
        # Determine first all configuration 'blocks' and from them the individual steps of the cascade
        wc = Cascade.generateBlocks(comp, wb, basis.orbitals, sa="for the decay cascade:")
        Cascade.displayBlocks(stdout, wc, sa="for the decay cascade ")
        if  printSummary   Cascade.displayBlocks(iostream, wc, sa="for the decay cascade ")    end 
        gMultiplets = Multiplet[];     for block in wc  push!(gMultiplets, block.multiplet)    end
        # Determine, modify and compute the transition data for all steps, ie. the PhotoEmission.Line's, the AutoIonization.Line's, etc.
        wd = Cascade.determineSteps(scheme, comp, wc)
        Cascade.displaySteps(stdout, wd, sa="decay ")
        if  printSummary   Cascade.displaySteps(iostream, wd, sa="decay ")    end      
        we   = Cascade.modifySteps(wd)
        @time data = Cascade.computeSteps(scheme, comp, we)
        if output    
            results = Base.merge( results, Dict("name"                  => comp.name) ) 
            results = Base.merge( results, Dict("cascade scheme"        => comp.scheme) ) 
            results = Base.merge( results, Dict("initial multiplets:"   => multiplets) )    
            results = Base.merge( results, Dict("generated multiplets:" => gMultiplets) )    
            results = Base.merge( results, Dict("decay line data:"      => data) )
            #
            #  Write out the result to file to later continue with simulations on the cascade data
            if outputToFile
                filename = "zzz-cascade-decay-computations-" * string(Dates.now())[1:13] * ".jld"
                println("\n* Write all results to disk; use:\n   JLD.save(''$filename'', results) \n   using JLD " *
                        "\n   results = JLD.load(''$filename'')    ... to load the results back from file.")
                if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD.save(''$filename'', results) \n   using JLD " *
                                                     "\n   results = JLD.load(''$filename'')    ... to load the results back from file." )      end      
                ## JLD.save(filename, results)
                JLD2.@save filename results
            end
        else    results = nothing
        end
        
        return( results )
    end
