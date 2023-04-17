
    # Functions and methods for cascade computation


    """
    `Cascade.computeSteps(scheme::Cascade.PhotoAbsorptionScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
        ... computes in turn all the requested transition amplitudes as well as PhotoExcitation.Line's, etc. for all pre-specified 
            excitation steps of the cascade. When compared with standard excitation computations, however, the amount of output is 
            largely reduced and often just printed into the summary file. A set of  data::Cascade.ExcitationData  is returned.
    """
    function computeSteps(scheme::Cascade.PhotoAbsorptionScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
        linesE = PhotoExcitation.Line[];    linesP = PhotoIonization.Line[]    
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        nt = 0;   st = 0;   previousMeanEn = 0.
        for  step  in  stepList
            st = st + 1
            nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
            println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc excitation lines (without selection rules): ")
            if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                 "up to $nc excitation lines (without selection rules): ")   end 
                                                 
            if  step.process == Basics.Photo() 
                newLines = PhotoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                               step.settings, output=true, printout=false) 
                append!(linesP, newLines);    nt = length(linesP)
            elseif  step.process == Basics.PhotoExc() 
                newLines = PhotoExcitation.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.grid, 
                                                               step.settings, output=true, printout=false) 
                append!(linesE, newLines);    nt = length(linesE)
            else   error("Unsupported atomic process for excitation computations.")
            end
            println("     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                    "to a total of $nt $(string(step.process)) decay lines." )
            if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                 "giving now rise to a total of $nt $(string(step.process)) decay lines." )   end      
        end
        #
        data = [ Cascade.PhotoIonizationData(linesP),  Cascade.ExcitationData(linesE) ]
    end


    """
    `Cascade.determineSteps(scheme::Cascade.PhotoAbsorptionScheme, comp::Cascade.Computation, 
                            initialList::Array{Cascade.Block,1}, ionizedList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1})`  
        ... determines all step::Cascade.Step's that need to be computed for this decay cascade. It cycles through all processes of the given
            decay scheme and selects all pairs of blocks due to the selected processes and cascade approach. It checks that at least on 
            pair of levels supports either a `photo-ionization' or `photo-excitation' within the step. A stepList::Array{Cascade.Step,1} 
            is returned, and for which subsequently all required transition amplitudes and rates/cross sections are computed.
    """
    function determineSteps(scheme::Cascade.PhotoAbsorptionScheme, comp::Cascade.Computation, 
                            initialList::Array{Cascade.Block,1}, ionizedList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1})
        #
        stepList = Cascade.Step[]
        if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
            for  initialBlock in initialList
                for  ionizedBlock in ionizedList
                    for  process  in  comp.scheme.processes
                        if      process == Basics.Photo()  
                            if  initialBlock.NoElectrons == ionizedBlock.NoElectrons + 1
                                photonEnergies = scheme.photonEnergies
                                println(">>> Photon energies must still be given in user-selected units: $(photonEnergies)")
                                settings = PhotoIonization.Settings(scheme.multipoles, [UseCoulomb, UseBabushkin], photonEnergies, Float64[],
                                                                    false, false, false, false, LineSelection(), Basics.ExpStokes() )
                                push!( stepList, Cascade.Step(process, settings, initialBlock.confs, ionizedBlock.confs, 
                                                                                 initialBlock.multiplet, ionizedBlock.multiplet) )
                            end
                        elseif  process == Basics.PhotoExc()
                        else    error("stop a")
                        end
                    end
                end
            end
            #
            for  initialBlock in initialList
                for  excitedBlock in excitedList
                    for  process  in  comp.scheme.processes
                        if      process == Basics.PhotoExc()  
                            if  initialBlock.NoElectrons == excitedBlock.NoElectrons
                                settings = PhotoExcitation.Settings(scheme.multipoles, [UseCoulomb, UseBabushkin], 
                                                           false, false, false, false, LineSelection(), 0., 0., 1.0e6, Basics.ExpStokes() )
                                push!( stepList, Cascade.Step(process, settings, initialBlock.confs, excitedBlock.confs, 
                                                                                 initialBlock.multiplet, excitedBlock.multiplet) )
                            end
                        elseif  process == Basics.Photo()
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
    `Cascade.generateBlocks(scheme::Cascade.PhotoAbsorptionScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
        ... generate all block::Cascade.Block's, that need to be computed for this excitation cascade, and compute also the corresponding multiplets.
            The different cascade approaches realizes different strategies how these block are selected and computed. 
            A blockList::Array{Cascade.Block,1} is returned.
    """
    function generateBlocks(scheme::Cascade.PhotoAbsorptionScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
        blockList = Cascade.Block[]
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        if    comp.approach == AverageSCA()
            if  printout
            println("\n* Generate blocks for photoabsorption computations:")
            println("\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println("    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; ")
            println("    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
            println("    + only E1 excitations are considered. \n")
            if  printSummary   
            println(iostream, "\n* Generate blocks for photoabsorption computations:")
            println(iostream, "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println(iostream, "    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; ")
            println(iostream, "    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
            println(iostream, "    + only E1 excitations are considered. \n")
            end
            end     # printout
            #
            for  confa  in confs
                print("  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")
                if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")   end
                    basis     = Basics.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                    multiplet = Basics.perform("computation: mutiplet from orbitals, no CI, CSF diagonal", [confa],  basis.orbitals, 
                                               comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
                println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
            end
        elseif    comp.approach == SCA()  error("Not yet implemented.")
        else  error("Unsupported cascade approach.")
        end

        return( blockList )
    end


    """
    `Cascade.generateConfigurationsForPhotoabsorption(multiplets::Array{Multiplet,1},  scheme::PhotoAbsorptionScheme, nm::Nuclear.Model)`  
        ... generates all possible (photo-absorption) configurations with a single displacements of an electron from
            scheme.excitationFromShells  to either scheme.excitationToShells or into partial waves with scheme.lValues.
            A Tuple(initialConfList::Array{Configuration,1}, confList::Array{Configuration,1}) is returned.
    """
    function generateConfigurationsForPhotoabsorption(multiplets::Array{Multiplet,1},  scheme::PhotoAbsorptionScheme, nm::Nuclear.Model)
        # Determine all (reference) configurations from multiplets and generate the 'excited/ionized' configurations due to the 
        # specificed excitation/ionization
        initialConfList = Configuration[];   ionConfList = Configuration[];   excConfList = Configuration[];
        for mp  in  multiplets   
            confList = Basics.extractNonrelativisticConfigurations(mp.levels[1].basis)
            for  conf in confList   if  conf in initialConfList   nothing   else   push!(initialConfList, conf)      end      end
        end
        # 
        # Generate all photoionized configurations if photoionization is to be considered; no configuration need to be excluded
        # since parity is give by the partial waves.
        if  Photo() in scheme.processes
            ionConfList = Basics.generateConfigurationsWithElectronLoss(initialConfList, scheme.excitationFromShells)
            ionConfList = unique(ionConfList)
        end
        # 
        # Generate all photoexcited configurations if photoexcitation is to be considered
        if  PhotoExc() in scheme.processes
            excConfList = Basics.generateConfigurations(initialConfList, scheme.excitationFromShells, scheme.excitationToShells, 1)
            #==
            # Exclude configurations with too low or too high mean energies as well as those that are parity forbidden for the given multipoles
            hasPlus = false;   hasMinus = false
            for  conf  in  initialConfList   
                if  Basics.determineParity(conf) == Basics.plus   hasPlus = true    else   hasMinus = true      end
            end
            en     = Float64[];   for conf in initialConfList    push!(en, -Semiempirical.estimate("binding energy", round(Int64, nm.Z), conf))   end
            maxen  = maximum(en);    minen  = minimum(en);  
            println(">>> initial configuration(s) have energies from $minen  to  $maxen  [a.u.].")
            #
            newExcConfList = Configuration[]
            for  conf  in  excConfList    meanEnergy = -Semiempirical.estimate("binding energy", round(Int64, nm.Z), conf) 
                if  minen + scheme.photonEnergies[1]  <= meanEnergy <= maxen + scheme.photonEnergies[end]
                    if      hasPlus  &&  hasMinus                                                                    push!(newExcConfList, conf)   
                    elseif  M1 in scheme.multipoles  ||  E2 in scheme.multipoles  ||  M2 in scheme.multipoles        push!(newExcConfList, conf)  
                    elseif  hasPlus  &&  E1 in scheme.multipoles  &&   Basics.determineParity(conf) == Basics.minus  push!(newExcConfList, conf)  
                    elseif  hasMinus &&  E1 in scheme.multipoles  &&   Basics.determineParity(conf) == Basics.plus   push!(newExcConfList, conf)  
                    else    println(">>> exclude $conf because of parity reasons.")
                    end
                else        println(">>> exclude $conf with energy $meanEnergy [a.u.] because of energy reasons.")
                end
            end  ==#
            #
            excConfList = unique(excConfList)
        end
        #
        return( initialConfList, ionConfList, excConfList )
    end


    """
    `Cascade.perform(scheme::PhotoAbsorptionScheme, comp::Cascade.Computation)`  
        ... to set-up and perform a photo-excitation computation that starts from a given set of initial configurations xor initial multiplets
            and comprises (various) photoexcitation processes into configurations with up-to NoExcitations single-electron excitations with 
            regard to the initial multiplets. The results of these excitation are comprised into (output) data::PhotoExcData, while these 
            data are only printed during the generation and nothing is returned.

    `Cascade.perform(scheme::PhotoAbsorptionScheme, comp::Cascade.Computation; output=true, outputToFile::Bool=true)`   
        ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
            cascade simulation. The particular output depends on the specifications of the cascade.
    """
    function perform(scheme::PhotoAbsorptionScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)
        ##
        ##
        if  output    results = Dict{String, Any}()    else    results = nothing    end
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        # Perform the SCF and CI computation for the intial-state multiplets if initial configurations are given
        if  comp.initialConfigs != Configuration[]
            basis      = Basics.performSCF(comp.initialConfigs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            multiplet  = Basics.performCI(basis, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            multiplets = [Multiplet("initial states", multiplet.levels)]
        else
            multiplets = comp.initialMultiplets
        end
        # Print out initial configurations and levels 
        Cascade.displayLevels(stdout, multiplets, sa="initial ")
        if  printSummary   Cascade.displayLevels(iostream, multiplets, sa="initial ")                            end
        #
        # Generate subsequent cascade configurations as well as display and group them together
        wa  = Cascade.generateConfigurationsForPhotoabsorption(multiplets, comp.scheme, comp.nuclearModel)
        wb1 = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa[1], sa="(initial part of the) photoabsorption ")
        wb2 = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa[2], sa="(photo-ionized part of the) photoabsorption ")
        wb3 = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa[3], sa="(photo-excited part of the) photoabsorption ")
        #
        # Determine first all configuration 'blocks' and from them the individual steps of the cascade
        wc1 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb1)
        wc2 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb2, printout=false)
        wc3 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb3, printout=false)
        Cascade.displayBlocks(stdout, wc1, sa="for the (initial part of the) photoabsorption cascade ");      
        Cascade.displayBlocks(stdout, wc2, sa="for the (photo-ionized part of the) photoabsorption cascade ")
        Cascade.displayBlocks(stdout, wc3, sa="for the (photo-excited part of the) photoabsorption cascade ")
        if  printSummary    Cascade.displayBlocks(iostream, wc1, sa="for the (initial part of the) photoabsorption cascade ");      
                            Cascade.displayBlocks(iostream, wc2, sa="for the (photo-ionized part of the) photoabsorption cascade ")
                            Cascade.displayBlocks(iostream, wc3, sa="for the (photo-excited part of the) photoabsorption cascade ")    end      
        # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
        gMultiplets = Multiplet[];     for block in wc2  push!(gMultiplets, block.multiplet)    end
                                       for block in wc3  push!(gMultiplets, block.multiplet)    end
        we = Cascade.determineSteps(scheme, comp, wc1, wc2, wc3)
        Cascade.displaySteps(stdout, we, sa="ionized & excited ")
        if  printSummary   Cascade.displaySteps(iostream, we, sa="ionized & excited ")    end      
        wf   = Cascade.modifySteps(we)
        data = Cascade.computeSteps(scheme, comp, wf)
        if output    
            results = Base.merge( results, Dict("name"                          => comp.name) ) 
            results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
            results = Base.merge( results, Dict("initial multiplets:"           => multiplets) )    
            results = Base.merge( results, Dict("generated multiplets:"         => gMultiplets) )    
            results = Base.merge( results, Dict("photo-ionization line data:"   => data[1]) )
            results = Base.merge( results, Dict("photo-excitation line data:"   => data[2]) )
            #
            #  Write out the result to file to later continue with simulations on the cascade data
            filename = "zzz-cascade-excitation-computations-" * string(Dates.now())[1:13] * ".jld"
            println("\n* Write all results to disk; use:\n   JLD.save(''$filename'', results) \n   using JLD " *
                    "\n   results = JLD.load(''$filename'')    ... to load the results back from file.")
            if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD.save(''$filename'', results) \n   using JLD " *
                                                 "\n   results = JLD.load(''$filename'')    ... to load the results back from file." )      end      
            JLD2.@save filename results
        end
        ## return( results )
        return( results )
    end
