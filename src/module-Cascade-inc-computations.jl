
    # Functions and methods for cascade computation

    """
    `Cascade.computeDecayYieldOutcome(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                      linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)` 
        ... to compute the flourescence and Auger yields for a single outcome as specified by its level; an 
            outcome::DecayYield.Outcome is returned in which are physical parameters are now specified.
    """
    function computeDecayYieldOutcome(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                      linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)
        # Identify the level key of the given level also in the lists of radiative and Auger lines
        level = outcome.level;    levelKey = LevelKey( LevelSymmetry(level.J, level.parity), level.index, level.energy, 0.)
        similarKey = LevelKey();  rateR = 0.;    rateA = 0.;   NoPhotonLines = 0;   NoAugerLines = 0
        ##x println("** levelKey = $levelKey")
        for  line in linesR
            compareKey = LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
            if   Basics.isSimilar(levelKey, compareKey, 1.0e-3)    println("** compareKey = $compareKey");   similarKey = deepcopy(compareKey)    end
        end
        if   similarKey == LevelKey()    error("No similar level found !")   end
        
        for  line in linesR
            ##x println("** Radiative: similarKey =$similarKey     $(line.initialLevel)")
            if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
                rateR = rateR + line.photonRate;   NoPhotonLines = NoPhotonLines + 1    
                ##x println("NoPhotonLines = $NoPhotonLines")
            end
        end
        for  line in linesA
            ##x println("** Auger: similarKey =$similarKey     $(line.initialLevel)")
            if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
                rateA = rateA + line.totalRate;   NoAugerLines = NoAugerLines + 1    
            end
        end
        
        omegaR = rateR / (rateR + rateA);   omegaA = rateA / (rateR + rateA)
        newOutcome = DecayYield.Outcome(level, NoPhotonLines, NoAugerLines, rateR, rateA, omegaR, omegaA)
        return( newOutcome )
    end


    """
    `Cascade.computeSteps(comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
        ... computes in turn all the requested transition amplitudes and PhotoEmission.Line's, AutoIonization.Line's, 
            etc. for all pre-specified decay steps of the cascade. When compared with the standard atomic process computations, 
            however, the amount of output is largely reduced. A set of  data::Cascade.Data  is returned.
    """
    function computeSteps(comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
        linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[];    linesP = PhotoIonization.Line[]    
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        nt = 0;   st = 0
        for  step  in  stepList
            st = st + 1
            nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels) 
            println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc lines (without selection rules): ")
            if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                 "up to $nc lines (without selection rules): ")   end      
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
                    "to a total of $nt $(string(step.process)) lines." )
            if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                 "giving now rise to a total of $nt $(string(step.process)) lines." )   end      
        end
        #
        data = Cascade.Data(comp.name, linesR, linesA, linesP)
    end


    """
    `Cascade.determineSteps(comp::Cascade.Computation, blockList::Array{Cascade.Block,1})`  
        ... determines all steps::Cascade.Step that need to be computed for this cascade. It cycles through the given processes and 
            distinguished between the different cascade approaches. It also checks that the averaged energies of the configuration 
            allows such a step energetically. A stepList::Array{Cascade.Step,1} is returned
    """
    function determineSteps(comp::Cascade.Computation, blockList::Array{Cascade.Block,1})
        stepList = Cascade.Step[]
        ##x println("comp.approach = $(comp.approach)")
        if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
            for  a = 1:length(blockList)
                for  b = 1:length(blockList)
                    minEn = 100000.;   maxEn = -100000.
                    for  p = 1:length(blockList[a].multiplet.levels),  q = 1:length(blockList[b].multiplet.levels)
                        minEn = min(minEn, blockList[a].multiplet.levels[p].energy - blockList[b].multiplet.levels[q].energy)
                        maxEn = max(maxEn, blockList[a].multiplet.levels[p].energy - blockList[b].multiplet.levels[q].energy)
                    end
                    for  process  in  comp.processes
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
    `Cascade.displayBlocks(stream::IO, blockList::Array{Cascade.Block,1})` 
        ... group & display the blocks of the cascade with same No. of electrons; this blocks are displayed with the
            minimum and maximum energy of each multiplet.
    """
    function displayBlocks(stream::IO, blockList::Array{Cascade.Block,1})
        #
        println(stream, "\n* Configuration 'blocks' (multiplets) in the given cascade model: \n")
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
    `Cascade.displayInitialLevels(stream::IO, multiplet::Multiplet)`  
        ... display the calculated initial levels to screen together with their given relative occupation.
    """
    function displayInitialLevels(stream::IO, multiplet::Multiplet)
        println(stream, " ")
        println(stream, "* Initial levels of the given cascade, relative to the lowest:")
        ##x println(stream, " ")
        println(stream, "  ", TableStrings.hLine(44))
        println(stream, "    Level  J Parity          Energy " * TableStrings.inUnits("energy") ) 
        println(stream, "  ", TableStrings.hLine(44))
        for  i = 1:length(multiplet.levels)
            lev = multiplet.levels[i]
            en  = lev.energy - multiplet.levels[1].energy;    en_requested = Defaults.convertUnits("energy: from atomic", en)
            ##x wx = 0.
            ##x for  ilevel in initialLevels
            ##x if  i in ilevel   wx = ilevel[2];    break    end
            ##x end
            ##x sb  = "          "   * string(wx)
            sc  = "   "  * TableStrings.level(i) * "     " * string(LevelSymmetry(lev.J, lev.parity)) * "     "
            @printf(stream, "%s %.15e %s", sc, en_requested, "\n")
        end
        println(stream, "  ", TableStrings.hLine(44))
        return( nothing )
    end
   

    """
    `Cascade.displaySteps(stream::IO, steps::Array{Cascade.Step,1})` 
        ... displays all predefined steps in a neat table and supports to delete individual steps from the list.
    """
    function displaySteps(stream::IO, steps::Array{Cascade.Step,1})
        println(stream, " ")
        println(stream, "* Steps that are defined for the current cascade due to the given approach:")
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
    `Cascade.generateBlocks(comp::Cascade.Computation, confs::Array{Configuration,1})`  
        ... generate all block::Cascade.Block that need to be computed for this cascade and compute the corresponding multiplets.
            The different cascade approches follow different strategies in defining and computing these blocks. 
            A blockList::Array{Cascade.Block,1} is returned.
    """
    function generateBlocks(comp::Cascade.Computation, confs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital})
        blockList = Cascade.Block[]
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        if    comp.approach == AverageSCA()
            println("\n* In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println("    + orbitals from the initial multiplet are applied throughout; ")
            println("    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
            println("    + only E1 dipole transitions are applied in all radiative decay stets; ")
            println("    + for each decay step, a (single) set of continuum orbitals with an configuration averaged energy is applied " *
                           "for all transitions of the same step. \n")
            if  printSummary   
            println(iostream, "\n* In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println(iostream, "    + orbitals from the initial multiplet are applied throughout; ")
            println(iostream, "    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
            println(iostream, "    + only E1 dipole transitions are applied in all radiative decay stets; ")
            println(iostream, "    + for each decay step, a (single) set of continuum orbitals with an configuration averaged energy is applied " *
                              "for all transitions of the same step. \n")
            end
            #
            for  confa  in confs
                print("  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")
                if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")   end
                multiplet = perform("computation: mutiplet from orbitals, no CI, CSF diagonal", [confa],  initalOrbitals, 
                                    comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
                println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
                ##x if  printSummary   println(iostream, "* ... and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")   end
            end
        elseif    comp.approach == SCA()
            println("\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println("    + orbitals are generated independently for each multiplet (block); ")
            println("    + configuration interaction is included for each block; ")
            println("    + only E1 dipole transitions are applied in all radiative decay stets; \n")
            if  printSummary   
            println(iostream, "\n* In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println(iostream, "    + orbitals are generated independently for each multiplet (block); ")
            println(iostream, "    + configuration interaction is included for each block; ")
            println(iostream, "    + only E1 dipole transitions are applied in all radiative decay stets; \n")
            end
            #
            i = 0
            for  confa  in confs
                ## i = i + 1;    if   i in [1,2, 4,5,6,7,8,9,10,11,12,13,14]  ||  i > 15   println("  Block $i omitted.");    continue    end
                ## i = i + 1;    if   i < 11  ||  i > 11   println("  Block $i omitted.");    continue    end
                print("  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")
                if  printSummary   print(iostream, "* Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")   end
                basis     = Basics.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                ##x multiplet = perform("computation: CI",  basis, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
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
    `Cascade.generateConfigurationList(initialConfs::Array{Configuration,1}, further::Int64, NoShake::Int64)`  
        ... generates all possible (decay) configurations with up to further holes and with NoShake displacements. First, all 
            configuratons are generated for which the hole is either moved 'outwards' or is moved and a second 'outer' hole is 
            created; this step is repated further + 2 times. From the generated list, only those configurations are kept with 
            up to further holes, when compared with the initial configuration. A confList::Array{Configuration,1} is returned.
    """
    function generateConfigurationList(initialConfs::Array{Configuration,1}, further::Int64, NoShake::Int64)
        confList = copy(initialConfs);    cList = copy(initialConfs);   initialNoElectrons = initialConfs[1].NoElectrons
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
    `Cascade.groupDisplayConfigurationList(Z::Float64, confs::Array{Configuration,1})` 
        ... group & display the configuration list into sublists with the same No. of electrons; this lists are displayed together 
            with an estimated total energy.
    """
    function groupDisplayConfigurationList(Z::Float64, confs::Array{Configuration,1})
        minNoElectrons = 1000;   maxNoElectrons = 0  
        for  conf in confs
            minNoElectrons = min(minNoElectrons, conf.NoElectrons)
            maxNoElectrons = max(maxNoElectrons, conf.NoElectrons)
        end
        #
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        println("\n* Electron configuration used in the cascade:")
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
                    sa = "   av. BE = "  * string( round(-wa) ) * "  " * TableStrings.inUnits("energy")
                    println("      " * string(conf) * sa * "      ($nc)" )
                    if  printSummary   println(iostream, "      " * string(conf) * sa * "      ($nc)")      end
                end  
            end
        end
        
        println("\n  A total of $nc configuration have been defined for this cascade, and selected configurations could be " *
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
