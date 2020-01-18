
    # Functions and methods for cascade simulations

    """
    `Cascade.displayIonDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})` 
        ... displays the (current or final) ion distribution in a neat table. Nothing is returned.
    """
    function displayIonDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})
        minElectrons = 1000;   maxElectrons = 0
        for  level in levels   minElectrons = min(minElectrons, level.NoElectrons);   maxElectrons = max(maxElectrons, level.NoElectrons)   end
        println(stream, " ")
        println(stream, "  (Final) Ion distribution for the cascade:  $sc ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(31))
        sa = "  "
        sa = sa * TableStrings.center(14, "No. electrons"; na=4)        
        sa = sa * TableStrings.center(10,"Rel. occ.";      na=2)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(31))
        for n = maxElectrons:-1:minElectrons
            sa = "             " * string(n);   sa = sa[end-10:end];   prob = 0.
            for  level in levels    if  n == level.NoElectrons   prob = prob + level.relativeOcc    end    end
            sa = sa * "         " * @sprintf("%.5e", prob)
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(31))

        return( nothing )
    end


    """
    `Cascade.displayLevelDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})` 
        ... displays the (current or final) level distribution in a neat table. Only those levels with a non-zero 
            occupation are displayed here. Nothing is returned.
    """
    function displayLevelDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})
        minElectrons = 1000;   maxElectrons = 0;   energies = zeros(length(levels))
        for  i = 1:length(levels)
            minElectrons = min(minElectrons, levels[i].NoElectrons);   maxElectrons = max(maxElectrons, levels[i].NoElectrons)
            energies[i]  = levels[i].energy   
        end
        enIndices = sortperm(energies, rev=true)
        # Now printout the results
        println(stream, " ")
        println(stream, "  (Final) Level distribution for the cascade:  $sc")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(69))
        sa = "  "
        sa = sa * TableStrings.center(14, "No. electrons"; na=2)        
        sa = sa * TableStrings.center( 8, "Lev-No"; na=2)        
        sa = sa * TableStrings.center( 6, "J^P"          ; na=3);               
        sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(10, "Rel. occ.";                                    na=2)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(69))
        for n = maxElectrons:-1:minElectrons
            sa = "            " * string(n);        sa  = sa[end-10:end]
            for  en in enIndices
                saa = "            " * string(en);  saa = saa[end-12:end]
                if  n == levels[en].NoElectrons  ##    &&  levels[en].relativeOcc > 0
                    sx = "    " * string( LevelSymmetry(levels[en].J, levels[en].parity) )                       * "           "
                    sb = sa * saa * sx[1:15]
                    sb = sb * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", levels[en].energy))  * "      "
                    sb = sb * @sprintf("%.5e", levels[en].relativeOcc) 
                    sa = "           "
                    println(stream, sb)
                end
            end
        end
        println(stream, "  ", TableStrings.hLine(69))

        return( nothing )
    end


    """
    `Cascade.displayLevelTree(stream::IO, levels::Array{Cascade.Level,1}, data::Cascade.Data; extended::Bool=false)` 
        ... displays all defined levels  in a neat table, together with their No. of electrons, symmetry, level energy, 
            current (relative) population as well as analogue information about their parents and daugther levels. This 
            enables one to recognize (and perhaps later add) missing parent and daughter levels. Nothing is returned.
    """
    function displayLevelTree(stream::IO, levels::Array{Cascade.Level,1}, data::Cascade.Data; extended::Bool=false)
        minElectrons = 1000;   maxElectrons = 0;   energies = zeros(length(levels))
        for  i = 1:length(levels)
            minElectrons = min(minElectrons, levels[i].NoElectrons);   maxElectrons = max(maxElectrons, levels[i].NoElectrons)
            energies[i]  = levels[i].energy   
        end
        enIndices = sortperm(energies, rev=true)
        # Now printout the results
        println(stream, " ")
        println(stream, "  Level tree of this cascade:  $(data.name)")
        println(stream, " ")
        if  extended    println(stream, "  ", TableStrings.hLine(179))  else    println(stream, "  ", TableStrings.hLine(65))  end
        sa = " "
        sa = sa * TableStrings.center( 6, "No. e-"; na=2)        
        sa = sa * TableStrings.center( 6, "Lev-No"; na=2)        
        sa = sa * TableStrings.center( 6, "J^P"          ; na=2);               
        sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(10, "Rel. occ.";                                na=3)
        if  extended
            sb = "Parents P(A/R: No_e, sym, energy) and Daughters D(A/R: No_e, sym, energy);  all energies in " * TableStrings.inUnits("energy")
            sa = sa * TableStrings.flushleft(100, sb; na=2)
        end
        #
        println(stream, sa)
        if  extended    println(stream, "  ", TableStrings.hLine(179))  else    println(stream, "  ", TableStrings.hLine(65))  end
        for n = maxElectrons:-1:minElectrons
            sa = "            " * string(n);     sa  = sa[end-5:end]
            for  en in enIndices
                saa = "         " * string(en);  saa = saa[end-8:end]
                if  n == levels[en].NoElectrons
                    sx = "    " * string( LevelSymmetry(levels[en].J, levels[en].parity) )                       * "           "
                    sb = sa * saa * sx[1:15]
                    sb = sb * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", levels[en].energy))  * "    "
                    sb = sb * @sprintf("%.4e", levels[en].relativeOcc)                                           * "  "
                    if extended
                    pProcessSymmetryEnergyList = Tuple{AtomicProcess,Int64,LevelSymmetry,Float64}[]
                    dProcessSymmetryEnergyList = Tuple{Basics.AtomicProcess,Int64,LevelSymmetry,Float64}[]
                    for  p in levels[en].parents
                        idx = p.index
                        if      p.process == Basics.Auger         lev = data.linesA[idx].initialLevel
                        elseif  p.process == Basics.Radiative     lev = data.linesR[idx].initialLevel
                        elseif  p.process == Basics.Photo         lev = data.linesP[idx].initialLevel
                        else    error("stop a")    end
                        push!( pProcessSymmetryEnergyList, (p.process, lev.basis.NoElectrons, LevelSymmetry(lev.J, lev.parity), lev.energy) )
                    end
                    for  d in levels[en].daugthers
                        idx = d.index
                        if      d.process == Basics.Auger         lev = data.linesA[idx].finalLevel
                        elseif  d.process == Basics.Radiative     lev = data.linesR[idx].finalLevel
                        elseif  d.process == Basics.Photo         lev = data.linesP[idx].finalLevel
                        else    error("stop b")    end
                        push!( dProcessSymmetryEnergyList, (d.process, lev.basis.NoElectrons, LevelSymmetry(lev.J, lev.parity), lev.energy) )
                    end
                    wa = TableStrings.processSymmetryEnergyTupels(120, pProcessSymmetryEnergyList, "P")
                    if  length(wa) > 0    sc = sb * wa[1];    println(stream,  sc )    else    println(stream,  sb )   end  
                    for  i = 2:length(wa)
                        sc = TableStrings.hBlank( length(sb) ) * wa[i];    println(stream,  sc )
                    end
                    wa = TableStrings.processSymmetryEnergyTupels(120, dProcessSymmetryEnergyList, "D")
                    for  i = 1:length(wa)
                        sc = TableStrings.hBlank( length(sb) ) * wa[i];    println(stream,  sc )
                    end
                    else    println(stream,  sb )
                    end  ## extended
                    sa = "      "
                end
            end
        end
        if  extended    println(stream, "  ", TableStrings.hLine(179))  else    println(stream, "  ", TableStrings.hLine(65))  end

        return( nothing )
    end


    """
    `Cascade.extractLevels(data::Cascade.Data, settings::Cascade.SimulationSettings)` 
        ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
            propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
            with its parent lines (which may populate the level) and the daugther lines (to which the pobability may decay). 
            A levelList::Array{Cascade.Level,1} is returned.
    """
    function extractLevels(data::Cascade.Data, settings::Cascade.SimulationSettings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        levels = Cascade.Level[]
        print("\n  Extract, sort and unify the list of levels of the cascade ... ")
        if printSummary     print(iostream, "\n  Extract, sort and unify the list of levels of the cascade ... ")     end
        
        for  i = 1:length(data.linesR)
            line = data.linesR[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(Basics.Radiative, i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(Basics.Radiative, i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end

        for  i = 1:length(data.linesA)
            line = data.linesA[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(Basics.Auger, i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(Basics.Auger, i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end

        for  i = 1:length(data.linesP)
            line = data.linesP[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(Basics.Photo, i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(Basics.Photo, i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end
        
        # Sort the levels by energy in reversed order
        energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
        enIndices = sortperm(energies, rev=true)
        newlevels = Cascade.Level[]
        for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
        
        # Assign the relative occupation of the levels in this list due to the given settings
        for  pair  in  settings.initialOccupations
            newlevels[pair[1]].relativeOcc = pair[2]
        end

        println("a total of $(length(newlevels)) levels were found.")
        println("\n    + all level numbers below refer and are sorted with regard to the overall cascade (i.e. all charged states involved);")
        if printSummary     
            println(iostream, "a total of $(length(newlevels)) levels were found.")
            println(iostream, "\n    + all level numbers below refer and are sorted with regard to the overall cascade (i.e. all charged states involved);")
        end
        
        return( newlevels )
    end


    """
    `Cascade.findLevelIndex(level::Cascade.Level, levels::Array{Cascade.Level,1})` 
        ... find the index of the given level within the given list of levels; an idx::Int64 is returned and an error message is 
            issued if the level is not found in the list.
    """
    function findLevelIndex(level::Cascade.Level, levels::Array{Cascade.Level,1})
        for  k = 1:length(levels)
            if  level.energy == levels[k].energy  &&   level.J == levels[k].J   &&   level.parity == levels[k].parity   &&
                level.NoElectrons == levels[k].NoElectrons
                kk = k;   return( kk )
            end
        end
        error("findLevelIndex():  No index was found;\n   level = $(level) ")
    end


    """
    `Cascade.propagateProbability!(levels::Array{Cascade.Level,1}, data::Cascade.Data)` 
        ... propagates the relative level occupation through the levels of the cascade until no further change occur in the 
            relative level occupation. The argument levels is modified during the propagation, but nothing is returned otherwise.
    """
    function propagateProbability!(levels::Array{Cascade.Level,1}, data::Cascade.Data)
        n = 0
        println("\n  Probability propagation through $(length(levels)) levels of the cascade:")
        while true
            n = n + 1;    totalProbability = 0.
            print("    $n-th round ... ")
            for  level in levels
                if   level.relativeOcc > 0.   && length(level.daugthers) > 0
                    # A level with relative occupation > 0 has still 'daugther' levels; collect all decay rates for this level
                    prob  = level.relativeOcc;   totalProbability = totalProbability + prob;   rates = zeros(length(level.daugthers))
                    level.relativeOcc = 0.
                    for  i = 1:length(level.daugthers)
                        idx = level.daugthers[i].index
                        if      level.daugthers[i].process == Basics.Radiative     rates[i] = data.linesR[idx].photonRate.Babushkin
                        elseif  level.daugthers[i].process == Basics.Auger         rates[i] = data.linesA[idx].totalRate
                        else    error("stop a; process = $(level.daugthers[i].process) ")
                        end
                    end
                    totalRate = sum(rates)
                    # Shift the relative occupation to the 'daugther' levels due to the different decay pathes
                    for  i = 1:length(level.daugthers)
                        idx = level.daugthers[i].index
                        if      level.daugthers[i].process == Basics.Radiative     line = data.linesR[idx]
                        elseif  level.daugthers[i].process == Basics.Auger         line = data.linesA[idx]
                        else    error("stop b; process = $(level.daugthers[i].process) ")
                        end
                        newLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, 
                                               line.finalLevel.basis.NoElectrons, 0., Cascade.LineIndex[], Cascade.LineIndex[] )
                        kk    = Cascade.findLevelIndex(newLevel, levels)
                        levels[kk].relativeOcc = levels[kk].relativeOcc + prob * rates[i] / totalRate
                    end
                end
            end
            println("has propagated a total of $totalProbability level occupation.")
            # Cycle once more if the relative occupation has still changed
            if  totalProbability == 0.    break    end
        end

        return( nothing )
    end


    """
    `Cascade.pushLevels!(levels::Array{Cascade.Level,1}, newLevel::Cascade.Level)` 
        ... push's the information of newLevel of levels. This is the standard 'push!(levels, newLevel)' if newLevel is not yet 
            including in levels, and the proper modification of the parent and daugther lines of this level otherwise. The argument 
            levels::Array{Cascade.Level,1} is modified and nothing is returned otherwise.
    """
    function pushLevels!(levels::Array{Cascade.Level,1}, newLevel::Cascade.Level)
        for  i = 1:length(levels)
            if  newLevel.energy == levels[i].energy  &&  newLevel.J == levels[i].J  &&  newLevel.parity == levels[i].parity
                append!(levels[i].parents,   newLevel.parents)
                append!(levels[i].daugthers, newLevel.daugthers)
                return( nothing )
            end
        end
        push!( levels, newLevel)
        ##x info("... one level added, n = $(length(levels)) ")
        return( nothing )
    end
    

    """
    `Cascade.simulateLevelDistribution(simulation::Cascade.Simulation)` 
        ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
            propagation, it runs through all levels and propagates the probabilty until no level probability changes anymore. The final 
            level distribution is then used to derive the ion distribution or the level distribution, if appropriate. Nothing is returned.
    """
    function simulateLevelDistribution(simulation::Cascade.Simulation)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        levels = Cascade.extractLevels(simulation.cascadeData, simulation.settings)
        Cascade.displayLevelTree(stdout, levels, simulation.cascadeData, extended=false)
        if  printSummary   Cascade.displayLevelTree(iostream, levels, simulation.cascadeData, extended=false)
                           Cascade.displayLevelTree(iostream, levels, simulation.cascadeData, extended=true)
        end
        Cascade.propagateProbability!(levels, simulation.cascadeData)
        #
        if  Cascade.IonDistribution()        in simulation.properties    
            Cascade.displayIonDistribution(stdout, simulation.cascadeData.name, levels)     
            if  printSummary   Cascade.displayIonDistribution(iostream, simulation.cascadeData.name, levels)      end
        end
        if  Cascade.FinalLevelDistribution() in simulation.properties    
            Cascade.displayLevelDistribution(stdout, simulation.cascadeData.name, levels)   
            if  printSummary   Cascade.displayLevelDistribution(iostream, simulation.cascadeData.name, levels)    end
        end

        return( nothing )
    end
