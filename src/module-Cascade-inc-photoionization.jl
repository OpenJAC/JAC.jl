
    # Functions and methods for cascade computation


    """
    `Cascade.computeSteps(scheme::Cascade.PhotonIonizationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
        ... computes in turn all the requested transition amplitudes as well as PhotoIonization.Line's, etc. for all pre-specified 
            (photo-) ionizing steps of the cascade. When compared with standard computations of these atomic processes, however, 
            the amount of output is largely reduced and often just printed into the summary file. A set of  
            data::Cascade.PhotoIonData  is returned.
    """
    function computeSteps(scheme::Cascade.PhotonIonizationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
        data = Cascade.PhotoIonData[]
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        
        # Perform computations in turn for all photon energies
        for  photonEnergy in  scheme.photonEnergies
            nt = 0;   st = 0;   linesP = PhotoIonization.Line[]  
            settings = PhotoIonization.Settings(PhotoIonization.Settings(), photonEnergies=[Defaults.convertUnits("energy: from atomic", photonEnergy)])
            println("\n\n* Perform ionizing computations for the photon energy $photonEnergy  [a.u.]")
            if  printSummary   println(iostream, "\n\n* Perform ionizing computations for the photon energy $photonEnergy  [a.u.]")   end 
            #
            for  step  in  stepList
                st = st + 1
                nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
                println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc ionizing lines (without selection rules): ")
                if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                     "up to $nc ionizing lines (without selection rules): ")   end 
          
                if      step.process == Basics.Photo()
                    newLines = PhotoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                                   settings, output=true, printout=true) 
                    append!(linesP, newLines);    nt = length(linesP)
                else   error("Unsupported atomic process for cascade computations.")
                end
                println("     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                        "to a total of $nt $(string(step.process)) ionizing lines." )
                if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                     "giving now rise to a total of $nt $(string(step.process)) ionizing lines." )   end      
            end
            #
            push!(data, Cascade.PhotoIonData(photonEnergy, linesP))
        end
        
        return( data )
    end


    """
    `Cascade.determineSteps(scheme::Cascade.PhotonIonizationScheme, comp::Cascade.Computation, 
                            initialList::Array{Cascade.Block,1}, ionizedList::Array{Cascade.Block,1})`  
        ... determines all steps::Cascade.Step's that need to be computed for this ionizing cascade. It cycles through all processes 
            of the given ionizing scheme and selects all pairs of initial and ionized blocks due to the selected photon energies and 
            cascade approaches.             It also checks that the (averaged) energies of the configuration support a photoionization 
            energetically.             A stepList::Array{Cascade.Step,1} is returned
    """
    function determineSteps(scheme::Cascade.PhotonIonizationScheme, comp::Cascade.Computation, 
                            initialList::Array{Cascade.Block,1}, ionizedList::Array{Cascade.Block,1})
        stepList = Cascade.Step[]
        if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
            for  initialBlock in initialList
                for  ionizedBlock in ionizedList
                    for  process  in  comp.scheme.processes
                        if      process == Basics.Photo()  
                            if  initialBlock.NoElectrons == ionizedBlock.NoElectrons + 1
                                settings = PhotoIonization.Settings()   ##x settings, photonEnergies=[1.0], printBefore=false)
                                push!( stepList, Cascade.Step(process, settings, initialBlock.confs, ionizedBlock.confs, 
                                                                                 initialBlock.multiplet, ionizedBlock.multiplet) )
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
    `Cascade.generateConfigurationsForPhotoionization(multiplets::Array{Multiplet,1},  maxPhotoElectrons::Int64, 
                                                      maxPhotonEnergy::Float64, nm::Nuclear.Model)`  
        ... generates all possible (photoionized) configurations with upto maxPhotoElectrons less electrons and whose averaged energy
            is NOT higher than maxPhotoEnergy above of the given configuration. In the present version, no shake contributions are taken
            into account. A Tuple(initialConfList::Array{Configuration,1}, confList::Array{Configuration,1}) is returned.
    """
    function generateConfigurationsForPhotoionization(multiplets::Array{Multiplet,1}, maxPhotoElectrons::Int64, maxPhotonEnergy::Float64, 
                                                      nm::Nuclear.Model)
        # Determine all (different) configurations from multiplets
        newconfList = Configuration[]
        for mp  in  multiplets   
            confList = Basics.extractNonrelativisticConfigurations(mp.levels[1].basis)
            for  conf in confList   if  conf in newconfList   nothing   else   push!(newconfList, conf)      end      end
        end
        en     = Float64[];   for conf in confList    push!(en, -Semiempirical.estimate("binding energy", round(Int64, nm.Z), conf))   end
        maxen  = maximum(en)
        shList = Basics.generate("shells: ordered list for NR configurations", newconfList)
        
        # Collect all configuration with up to maxPhotoElectrons less electrons
        initialConfList = deepcopy(newconfList);      
        allnewconfList  = Configuration[];    confList = deepcopy(newconfList);    newconfList = Configuration[]
        for ne = 1:maxPhotoElectrons
            for  sh in shList
                for conf in confList
                    shells = deepcopy(conf.shells);     
                    if shells[sh] > 0   shells[sh] = shells[sh] - 1;     push!(newconfList, Configuration(shells, conf.NoElectrons - 1))  end
                end
            end
            append!(allnewconfList, newconfList)
            confList = deepcopy(newconfList);    newconfList = Configuration[]
        end
        
        # Exclude double configurations as well as those with too high average energy
        ##x confList = Basics.excludeDoubles(allnewconfList)
        confList = unique(allnewconfList)
        newconfList = Configuration[]
        for  conf in confList  
            enconf = -Semiempirical.estimate("binding energy", round(Int64, nm.Z), conf)
            if enconf - maxen < maxPhotonEnergy   push!(newconfList, conf)    end
        end

        return( (initialConfList, newconfList)  )
    end


    """
    `Cascade.perform(scheme::PhotonIonizationScheme, comp::Cascade.Computation)`  
        ... to set-up and perform a cascade computation that starts from a given set of initial configurations xor initial multiplets
            and proceeds via photoionizing processes to configurations with a given No of photoelectrons. The results of all individual 
            photoionization steps are comprised into (output) data::PhotoIonData, although all data are only printed and nothing 
            is returned.

    `Cascade.perform(scheme::PhotonIonizationScheme, comp::Cascade.Computation; output=true, outputToFile::Bool=true)`   
        ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
            cascade simulation. The particular output also depends on the specifications of the cascade.
    """
    function perform(scheme::PhotonIonizationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)
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
        maxen = maximum(comp.scheme.photonEnergies) + 20
        wa  = Cascade.generateConfigurationsForPhotoionization(multiplets, comp.scheme.maxPhotoElectrons, maxen, comp.nuclearModel)
        wb1 = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa[1], sa="(initial part of the) ionizing ")
        wb2 = Cascade.groupDisplayConfigurationList(comp.nuclearModel.Z, wa[2], sa="(generated part of the) ionizing ")
        #
        # Determine first all configuration 'blocks' and from them the individual steps of the cascade
        wc1 = Cascade.generateBlocks(comp, wb1, basis.orbitals, sa="for the ionizing cascade ")
        wc2 = Cascade.generateBlocks(comp, wb2, basis.orbitals, sa="for the ionizing cascade ", printout=false)
        Cascade.displayBlocks(stdout, wc1, sa="for the (initial part of the) ionizing cascade ");      
        Cascade.displayBlocks(stdout, wc2, sa="for the (generated part of the) ionizing cascade ")
        if  printSummary   Cascade.displayBlocks(iostream, wc1, sa="for the (initial part of the) ionizing cascade ")
                           Cascade.displayBlocks(iostream, wc2, sa="for the (generated part of the) ionizing cascade ")    end      
        # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
        gMultiplets = Multiplet[];     for block in wc2  push!(gMultiplets, block.multiplet)    end
        we = Cascade.determineSteps(scheme, comp, wc1, wc2)
        Cascade.displaySteps(stdout, we, sa="ionizing ")
        if  printSummary   Cascade.displaySteps(iostream, we, sa="ionizing ")    end      
        wf   = Cascade.modifySteps(we)
        data = Cascade.computeSteps(scheme, comp, wf)
        if output    
            results = Base.merge( results, Dict("name"                       => comp.name) ) 
            results = Base.merge( results, Dict("cascade scheme"             => comp.scheme) ) 
            results = Base.merge( results, Dict("initial multiplets:"        => multiplets) )    
            results = Base.merge( results, Dict("generated multiplets:"      => gMultiplets) )    
            results = Base.merge( results, Dict("photo-ionizing line data:"  => data) )
            #
            #  Write out the result to file to later continue with simulations on the cascade data
            filename = "zzz-cascade-ionizing-computations-" * string(Dates.now())[1:13] * ".jld"
            println("\n* Write all results to disk; use:\n   JLD.save(''$filename'', results) \n   using JLD " *
                    "\n   results = JLD.load(''$filename'')    ... to load the results back from file.")
            if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD.save(''$filename'', results) \n   using JLD " *
                                                 "\n   results = JLD.load(''$filename'')    ... to load the results back from file." )      end      
            JLD2.@save filename results
        end
        ## return( results )
        return( results )
    end
