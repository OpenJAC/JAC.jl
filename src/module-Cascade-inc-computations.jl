
    # Functions and methods for cascade computation

    """
    `Cascade.computeDecayYieldOutcome(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                      linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)` 
        ... to compute the flourescence and Auger yields for a single outcome as specified by its level; an 
            outcome::DecayYield.Outcome is returned in which all physical parameters are now specified for the given
            decay-yield.
    """
    function computeDecayYieldOutcome(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                      linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)
        # Identify the level key of the given level also in the lists of radiative and Auger lines
        level = outcome.level;    levelKey = LevelKey( LevelSymmetry(level.J, level.parity), level.index, level.energy, 0.)
        similarKey = LevelKey();  rateR = 0.;    rateA = 0.;   NoPhotonLines = 0;   NoAugerLines = 0
        for  line in linesR
            compareKey = LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
            if   Basics.isSimilar(levelKey, compareKey, 1.0e-3)    println("** compareKey = $compareKey");   similarKey = deepcopy(compareKey)    end
        end
        if   similarKey == LevelKey()    error("No similar level found !")   end
        
        for  line in linesR
            if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
                rateR = rateR + line.photonRate;   NoPhotonLines = NoPhotonLines + 1    
            end
        end
        for  line in linesA
            if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
                rateA = rateA + line.totalRate;   NoAugerLines = NoAugerLines + 1    
            end
        end
        
        omegaR     = rateR / (rateR + rateA);   omegaA = rateA / (rateR + rateA)
        newOutcome = DecayYield.Outcome(level, NoPhotonLines, NoAugerLines, rateR, rateA, omegaR, omegaA)
        return( newOutcome )
    end


    """
    `Cascade.computeSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
        ... computes in turn all the requested transition amplitudes and PhotoEmission.Line's, AutoIonization.Line's, 
            etc. for all pre-specified decay steps of the cascade. When compared with standard computations of these atomic 
            processes, however, the amount of output is largely reduced and often just printed into the summary file. 
            A set of  data::Cascade.DecayData  is returned.
    """
    function computeSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
        linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[]    
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        nt = 0;   st = 0
        for  step  in  stepList
            st = st + 1
            nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
            println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc decay lines (without selection rules): ")
            if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                 "up to $nc decay lines (without selection rules): ")   end 
          
            if      step.process == Basics.Auger 
                newLines = AutoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                              step.settings, output=true, printout=false) 
                append!(linesA, newLines);    nt = length(linesA)
            elseif  step.process == Basics.Radiative
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
    `Cascade.computeSteps(scheme::Cascade.PhotonIonizationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
        ... computes in turn all the requested transition amplitudes and PhotoIonization.Line's, etc. for all pre-specified 
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
          
                if      step.process == Basics.Photo
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
                        if      process == Basics.Radiative   
                            if  a == b   ||   minEn < 0.    continue   end
                            if  blockList[a].NoElectrons == blockList[b].NoElectrons
                                settings = PhotoEmission.Settings([E1], [UseBabushkin], false, false, false, Tuple{Int64,Int64}[], 0., 0., 1.0e6)
                                push!( stepList, Cascade.Step(process, settings, blockList[a].confs, blockList[b].confs, 
                                                              blockList[a].multiplet, blockList[b].multiplet) )
                            end
                        elseif  process == Basics.Auger       
                            if  a == b   ||   minEn < 0.    continue   end
                            if  blockList[a].NoElectrons == blockList[b].NoElectrons + 1
                                settings = AutoIonization.Settings(false, false, false, Tuple{Int64,Int64}[], 0., 1.0e6, 2, "Coulomb")
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
                        if      process == Basics.Photo   
                            if  initialBlock.NoElectrons == ionizedBlock.NoElectrons + 1
                                ##x settings = PhotoIonization.Settings()
                                ##x for  en  in scheme.photonEnergies   push!(phEnergies, Defaults.convertUnits("energy: from atomic", en))   end
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
    `Cascade.displayBlocks(stream::IO, blockList::Array{Cascade.Block,1}; sa::String="")` 
        ... group & display the blocks of the cascade with same No. of electrons; this blocks are displayed with the
            minimum and maximum energy of each multiplet. The optional sa::String can be used to display some details
            about the given blocks. nothing is returned.
    """
    function displayBlocks(stream::IO, blockList::Array{Cascade.Block,1}; sa::String="")
        #
        println(stream, "\n* Configuration 'blocks' (multiplets) " * sa * "in the given (cascade) model: \n")
        println(stream, "  ", TableStrings.hLine(134))
        println(stream, "      No.   Configurations                                                                       " *
                        "      Range of total energies " * TableStrings.inUnits("energy") ) 
        println(stream, "  ", TableStrings.hLine(134))
        i = 0
        for  block  in  blockList
            i = i + 1;    
            sa = "   " * TableStrings.flushright( 6, string(i); na=2)
            sb = " ";         for conf  in blockList[i].confs   sb = sb * string(conf) * ", "    end
            en = Float64[];   for level in  block.multiplet.levels    push!(en, level.energy)    end
            minEn = minimum(en);   minEn = Defaults.convertUnits("energy: from atomic", minEn)
            maxEn = maximum(en);   maxEn = Defaults.convertUnits("energy: from atomic", maxEn)
            sa = sa * TableStrings.flushleft(90, sb[1:end-2]; na=2) 
            sa = sa * TableStrings.flushleft(30, string( round(minEn)) * " ... " * string( round(maxEn)); na=2)
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(134))

        return( nothing )
    end
   

    """
    `Cascade.displayLevels(stream::IO, multiplets::Array{Multiplet,1}; sa::String="")`  
        ... display on stream the initial configurations as well as the calculated levels for all initial multiplets.
    """
    function displayLevels(stream::IO, multiplets::Array{Multiplet,1}; sa::String="")
        println(stream, " ")
        println(stream, "* Configurations and levels for all given " * sa * "multiplets of the cascade, relative to the lowest:")
        for  multiplet  in multiplets
            println(stream, "  ")
            confList = Basics.extractNonrelativisticConfigurations(multiplet.levels[1].basis)
            for  conf in confList
                println(stream, "  $conf")
            end
            println(stream, "  ", TableStrings.hLine(44))
            println(stream, "    Level  J Parity          Energy " * TableStrings.inUnits("energy") ) 
            println(stream, "  ", TableStrings.hLine(44))
            for  i = 1:length(multiplet.levels)
                lev = multiplet.levels[i]
                en  = lev.energy - multiplet.levels[1].energy;    en_requested = Defaults.convertUnits("energy: from atomic", en)
                sc  = "   "  * TableStrings.level(i) * "     " * string(LevelSymmetry(lev.J, lev.parity)) * "     "
                @printf(stream, "%s %.15e %s", sc, en_requested, "\n")
            end
            println(stream, "  ", TableStrings.hLine(44))
        end
        return( nothing )
    end
   

    """
    `Cascade.displaySteps(stream::IO, steps::Array{Cascade.Step,1}; sa::String="")` 
        ... displays all predefined steps in a neat table and supports to delete individual steps from the list.
            sa::String can be used to display details about the given steps
    """
    function displaySteps(stream::IO, steps::Array{Cascade.Step,1}; sa::String="")
        println(stream, " ")
        println(stream, "* Steps that are defined for the current " * sa * "cascade due to the given approach:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(170))
        sa = "  "
        sa = sa * TableStrings.center( 9, "Step-No"; na=2)
        sa = sa * TableStrings.flushleft(11, "Process"; na=1)
        sa = sa * TableStrings.flushleft(55, "Initial:  No CSF, configuration(s)"; na=4)
        sa = sa * TableStrings.flushleft(55, "Final:  No CSF, configuration(s)"; na=4)
        sa = sa * TableStrings.flushleft(40, "Energies from ... to in " * TableStrings.inUnits("energy"); na=4)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(170))
        #
        for  i = 1:length(steps)
            sa = " " * TableStrings.flushright( 7, string(i); na=5)
            sa = sa  * TableStrings.flushleft( 11, string(steps[i].process); na=1)
            sb = "";   for conf in steps[i].initialConfigs   sb = sb * string(conf) * ", "    end
            sa = sa  * TableStrings.flushright( 5, string( length(steps[i].initialMultiplet.levels[1].basis.csfs) )*", "; na=0) 
            sa = sa  * TableStrings.flushleft( 50, sb[1:end-2]; na=4)
            sb = "";   for conf in steps[i].finalConfigs     sb = sb * string(conf) * ", "    end
            sa = sa  * TableStrings.flushright( 5, string( length(steps[i].finalMultiplet.levels[1].basis.csfs) )*", "; na=0) 
            sa = sa  * TableStrings.flushleft( 50, sb[1:end-2]; na=4)
            minEn = 1000.;   maxEn = -1000.;
            for  p = 1:length(steps[i].initialMultiplet.levels),  q = 1:length(steps[i].finalMultiplet.levels)
                minEn = min(minEn, steps[i].initialMultiplet.levels[p].energy - steps[i].finalMultiplet.levels[q].energy)
                maxEn = max(maxEn, steps[i].initialMultiplet.levels[p].energy - steps[i].finalMultiplet.levels[q].energy)
            end
            minEn = Defaults.convertUnits("energy: from atomic", minEn);   maxEn = Defaults.convertUnits("energy: from atomic", maxEn)
            sa = sa * string( round(minEn)) * " ... " * string( round(maxEn))
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(170))
    end
    
    
    
    """
    `Cascade.generateBlocks(comp::Cascade.Computation, confs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}; 
                            sa::String="", printout::Bool=true)`  
        ... generate all block::Cascade.Block's, that need to be computed for this cascade, and compute also the corresponding multiplets.
            The different cascade approches enables one to realized follow different strategies how these block are selected and computed. 
            A blockList::Array{Cascade.Block,1} is returned.
    """
    function generateBlocks(comp::Cascade.Computation, confs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}; 
                            sa::String="", printout::Bool=true)
        blockList = Cascade.Block[]
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        if    comp.approach == AverageSCA()
            if  printout
            println("\n* Generate blocks " * sa)
            println("\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println("    + orbitals from the initial multiplet are applied throughout; ")
            println("    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
            println("    + only E1 dipole transitions are applied in all radiative decay stets; ")
            println("    + for each decay step, a (single) set of continuum orbitals with an configuration averaged energy is applied " *
                           "for all transitions of the same step. \n")
            if  printSummary   
            println(iostream, "\n* Generate blocks " * sa)
            println(iostream, "\n* In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println(iostream, "    + orbitals from the initial multiplet are applied throughout; ")
            println(iostream, "    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
            println(iostream, "    + only E1 dipole transitions are applied in all radiative decay stets; ")
            println(iostream, "    + for each decay step, a (single) set of continuum orbitals with an configuration averaged energy is applied " *
                              "for all transitions of the same step. \n")
            end
            end     # printout
            #
            for  confa  in confs
                print("  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")
                if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")   end
                multiplet = Basics.perform("computation: mutiplet from orbitals, no CI, CSF diagonal", [confa],  initalOrbitals, 
                                           comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
                println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
                ##x if  printSummary   println(iostream, "* ... and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")   end
            end
        elseif    comp.approach == SCA()
            if printout
            println("\n* Generate blocks " * sa)
            println("\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println("    + orbitals are generated independently for each multiplet (block); ")
            println("    + configuration interaction is included for each block; ")
            println("    + only E1 dipole transitions are applied in all radiative decay stets; \n")
            if  printSummary   
            println(iostream, "\n* Generate blocks " * sa)
            println(iostream, "\n* In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println(iostream, "    + orbitals are generated independently for each multiplet (block); ")
            println(iostream, "    + configuration interaction is included for each block; ")
            println(iostream, "    + only E1 dipole transitions are applied in all radiative decay stets; \n")
            end
            end     # printout
            #
            i = 0
            for  confa  in confs
                ## i = i + 1;    if   i in [1,2, 4,5,6,7,8,9,10,11,12,13,14]  ||  i > 15   println("  Block $i omitted.");    continue    end
                ## i = i + 1;    if   i < 11  ||  i > 11   println("  Block $i omitted.");    continue    end
                print("  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")
                if  printSummary   print(iostream, "* Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")   end
                basis     = Basics.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                multiplet = Basics.performCI(basis, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
                println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
                ##x if  printSummary   println(iostream, "and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")   end
            end
        else  error("Unsupported cascade approach.")
        end

        return( blockList )
    end
    

    """
    `Cascade.generateConfigurationList(multiplets::Array{Multiplet,1}, further::Int64, NoShake::Int64)`  
        ... generates all possible (decay) configurations with up to further holes and with NoShake displacements with regard
            to the given multiplets. First, all configuratons are generated for which the hole is either moved 'outwards' or 
            is moved and a second 'outer' hole is created; this step is repated further + 2 times to make sure that all relevant
            configurations are met. From the generated list, however, only those configurations are kept eventually with 
            up to further holes, when compared to the configurations of the given multiplets. A confList::Array{Configuration,1} 
            is returned.
    """
    function generateConfigurationList(multiplets::Array{Multiplet,1}, further::Int64, NoShake::Int64)
        # Determine all (different) configurations from multiplets
        confList = Configuration[]
        for mp  in  multiplets   
            cfList = Basics.extractNonrelativisticConfigurations(mp.levels[1].basis)
            for  cf in cfList   if  cf in confList   nothing   else   push!(confList, cf)      end      end
        end
        cList = copy(confList);   initialNoElectrons = multiplets[1].levels[1].basis.NoElectrons
        # First, move and generate new 'outer' hole without displacements
        for  fur = 1:further+1
            newConfList = Configuration[]
            for conf  in cList
                holeList = Basics.determineHoleShells(conf)
                for  holeShell in holeList
                    wa = generateConfigurationsWith1OuterHole(conf,  holeShell);   append!(newConfList, wa)
                    wa = generateConfigurationsWith2OuterHoles(conf, holeShell);   append!(newConfList, wa)
                end
            end
            if  length(newConfList) > 0    newConfList = Basics.excludeDoubles(newConfList)    end
            cList = newConfList
            append!(confList, newConfList)
        end
        # Make sure that only configurations with up to further holes are returned
        newConfList = Configuration[]
        for   conf in confList   
            if  conf.NoElectrons + further >= initialNoElectrons   push!(newConfList, conf)    end
        end
        # Add further shake-displacements if appropriate
        newConfList = Basics.excludeDoubles(newConfList)
        return( newConfList )
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
        ##x println("xx  $(multiplets)   $(maxen)  $(shList)    $(newconfList)")
        
        # Collect all configuration with up to maxPhotoElectrons less electrons
        initialConfList = deepcopy(newconfList);      
        allnewconfList  = Configuration[];    confList = deepcopy(newconfList);    newconfList = Configuration[]
        ##x println("aa  $(length(confList))")
        for ne = 1:maxPhotoElectrons
            for  sh in shList
                for conf in confList
                    shells = deepcopy(conf.shells);     
                    if shells[sh] > 0   shells[sh] = shells[sh] - 1;     push!(newconfList, Configuration(shells, conf.NoElectrons - 1))  end
                end
            end
            append!(allnewconfList, newconfList)
            confList = deepcopy(newconfList);    newconfList = Configuration[]
            ##x println("bb  $(length(confList))")
        end
        
        # Exclude double configurations as well as those with too high average energy
        confList = Basics.excludeDoubles(allnewconfList)
        newconfList = Configuration[]
        for  conf in confList  
            enconf = -Semiempirical.estimate("binding energy", round(Int64, nm.Z), conf)
            ##x @show enconf, maxen, maxPhotonEnergy
            if enconf - maxen < maxPhotonEnergy   push!(newconfList, conf)    end
        end
        ##x println("cc  $(length(newconfList))")

        ##x @show length(initialConfList), length(newconfList)
        return( (initialConfList, newconfList)  )
    end


    """
    `Cascade.generateConfigurationsWith1OuterHole(conf,  holeShell)`  
        ... generates all possible (decay) configurations where the hole in holeShell is moved 'outwards'. 
            A confList::Array{Configuration,1} is returned.
    """
    function generateConfigurationsWith1OuterHole(conf::Configuration,  holeShell::Shell)
         shList = Basics.generate("shells: ordered list for NR configurations", [conf]);   i0 = 0
         for  i = 1:length(shList)
             if   holeShell == shList[i]    i0 = i;    break    end
         end
         if  i0 == 0   error("stop a")   end
         #
         # Now move the hole 'outwards'
         confList = Configuration[]
         for  i = i0+1:length(shList)
             if  haskey(conf.shells, shList[i])  &&  conf.shells[ shList[i] ] >= 1  
                 newshells = copy( conf.shells )
                 newshells[ shList[i] ] = newshells[ shList[i] ] - 1
                 newshells[ holeShell ] = newshells[ holeShell ] + 1
                 push!(confList, Configuration( newshells, conf.NoElectrons ) )
             end
         end
         return( confList )
    end


    """
    `Cascade.generateConfigurationsWith2OuterHoles(conf,  holeShell)`  
        ... generates all possible (decay) configurations where the hole in holeShell is moved 'outwards'. 
            A confList::Array{Configuration,1} is returned.
    """
    function generateConfigurationsWith2OuterHoles(conf::Configuration,  holeShell::Shell)
         shList = Basics.generate("shells: ordered list for NR configurations", [conf]);   i0 = 0
         for  i = 1:length(shList)
             if   holeShell == shList[i]    i0 = i;    break    end
         end
         if  i0 == 0   error("stop a")   end
         #
         # Now move the hole 'outwards'
         confList = Configuration[]
         for  i = i0+1:length(shList)
             if  haskey(conf.shells, shList[i])  &&  conf.shells[ shList[i] ] >= 2  
                 newshells = copy( conf.shells )
                 newshells[ shList[i] ] = newshells[ shList[i] ] - 2
                 newshells[ holeShell ] = newshells[ holeShell ] + 1
                 push!(confList, Configuration( newshells, conf.NoElectrons - 1 ) )
             end
             #
             for  j = i0+1:length(shList)
                 if  i != j   &&   haskey(conf.shells, shList[i])  &&  conf.shells[ shList[i] ] >= 1   &&
                                   haskey(conf.shells, shList[j])  &&  conf.shells[ shList[j] ] >= 1 
                     newshells = copy( conf.shells )
                     newshells[ shList[i] ] = newshells[ shList[i] ] - 1
                     newshells[ shList[j] ] = newshells[ shList[j] ] - 1
                     newshells[ holeShell ] = newshells[ holeShell ] + 1
                     push!(confList, Configuration( newshells, conf.NoElectrons - 1 ) )
                 end
             end
         end
         return( confList )
    end



    """
    `Cascade.groupDisplayConfigurationList(Z::Float64, confs::Array{Configuration,1}; sa::String="")` 
        ... group & display the configuration list into sublists with the same No. of electrons; this lists are displayed together 
            with an estimated total energy.
    """
    function groupDisplayConfigurationList(Z::Float64, confs::Array{Configuration,1}; sa::String="")
        minNoElectrons = 1000;   maxNoElectrons = 0  
        for  conf in confs
            minNoElectrons = min(minNoElectrons, conf.NoElectrons)
            maxNoElectrons = max(maxNoElectrons, conf.NoElectrons)
        end
        #
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        println("\n* Electron configuration used in the " * sa * "cascade:")
        @warn "*** Limit to just two configurations for each No. of electrons. ***"                       ## delete nxx
        if  printSummary   println(iostream, "\n* Electron configuration used in the cascade:")    end
        confList = Configuration[];   nc = 0
        for  n = maxNoElectrons:-1:minNoElectrons
            nxx = 0                                                                                        ## delete nxx
            println("\n  Configuration(s) with $n electrons:")
            if  printSummary   println(iostream, "\n    Configuration(s) with $n electrons:")      end
            for  conf in confs
                if n == conf.NoElectrons  
                    nxx = nxx + 1;    if nxx > 2   break    end                                            ## delete nxx
                    nc = nc + 1
                    push!(confList, conf ) 
                    wa = Semiempirical.estimate("binding energy", round(Int64, Z), conf);    wa = Defaults.convertUnits("energy: from atomic", wa)
                    sb = "   av. BE = "  * string( round(-wa) ) * "  " * TableStrings.inUnits("energy")
                    println("      " * string(conf) * sb * "      ($nc)" )
                    if  printSummary   println(iostream, "      " * string(conf) * sb * "      ($nc)")      end
                end  
            end
        end
        
        println("\n  A total of $nc configuration have been defined for this " * sa * "cascade, and selected configurations could be " *
                "removed here:  [currently not supported]")
        if  printSummary   println(iostream, "\n* A total of $nc configuration have been defined for this cascade, and selected " *
                                             "configurations could be removed here:  [currently not supported]")      end
        return( confList )
    end


    """
    `Cascade.modifySteps(stepList::Array{Cascade.Step,1})` 
        ... allows the user to modify the steps, for instance, by deleting selected steps of the cascade or by modifying the settings of
            one or several steps. A newStepList::Array{Cascade.Step,1} for which the transition data are eventually computed.
    """
    function modifySteps(stepList::Array{Cascade.Step,1})
        #
        newStepList = Cascade.Step[]
        #
        println("\n* Here, modify the individual steps explicitly in the code, if needed, ...... and just do it !!")
        # 
        #  Delete individual steps from stepList
        #  if  i in [1,2,5, ...] modify the particular settings, etc.
        for  i = 1:length(stepList)
            step = stepList[i]
            #
            if  i in []
                println("  Modify step $i :")
                newStep = Cascade.Step(step.process, step.settings, step.initialConfs, step.finalConfs, step.initialMultiplet, step.initialMultiplet)
                push!(newStepList, newStep)
            else
                push!(newStepList, step)
            end
        end
        #
        # wa = [1,2,3]
        # delete from list
        #
        println("\n  A total of $(length(newStepList)) steps are still defined in the cascade.")
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   println(iostream, "\n* A total of $(length(newStepList)) steps are still defined in the cascade.")    end      
        
        return( newStepList )
    end


    """
    `Cascade.perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation)`  
        ... to set-up and perform a cascade computation that starts from a given set of initial configurations and proceeds via 
            various steps until a given number of electrons has been removed or the decay stops at some stable levels with regard 
            to the given atomic processes. The results of all individual steps are printed to screen but nothing is returned 
            otherwise.

    `Cascade.perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation; output=true)`   
        ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the type 
            and specifications of the cascade but can easily accessed by the keys of this dictionary.
    """
    function perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation; output::Bool=false)
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
        data = Cascade.computeSteps(scheme, comp, we)
        if output    
            results = Base.merge( results, Dict("name"                  => comp.name) ) 
            results = Base.merge( results, Dict("cascade scheme"        => comp.scheme) ) 
            results = Base.merge( results, Dict("initial multiplets:"   => multiplets) )    
            results = Base.merge( results, Dict("generated multiplets:" => gMultiplets) )    
            results = Base.merge( results, Dict("decay line data:"      => data) )
            #
            #  Write out the result to file to later continue with simulations on the cascade data
            filename = "zzz-cascade-decay-computations-" * string(Dates.now())[1:13] * ".jld"
            println("\n* Write all results to disk; use:\n   JLD.save(''$filename'', results) \n   using JLD " *
                    "\n   results = JLD.load(''$filename'')    ... to load the results back from file.")
            if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD.save(''$filename'', results) \n   using JLD " *
                                                 "\n   results = JLD.load(''$filename'')    ... to load the results back from file." )      end      
            JLD.save(filename, results)
        else    results = nothing
        end
        
        return( results )
    end


    """
    `Cascade.perform(scheme::PhotonIonizationScheme, comp::Cascade.Computation)`  
        ... to set-up and perform a cascade computation that starts from a given set of initial configurations xor initial multiplets
            and proceeds via photoionizing processes to configurations with a given No of photoelectrons. The results of all individual 
            photoionization steps are comprised into (output) data::PhotoIonData, although all data are only printed and nothing 
            is returned.

    `Cascade.perform(scheme::PhotonIonizationScheme, comp::Cascade.Computation; output=true)`   
        ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
            cascade simulation. The particular output also depends on the specifications of the cascade.
    """
    function perform(scheme::PhotonIonizationScheme, comp::Cascade.Computation; output::Bool=false)
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
        ##x @show length(wa[1]), length(wa[2])
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
            JLD.save(filename, results)
        end
        ## return( results )
        return( data )
    end
