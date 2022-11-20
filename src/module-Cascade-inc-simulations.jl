
    # Functions and methods for cascade simulations

    """
    `Cascade.addLevels(levelsA::Array{Cascade.Level,1}, levelsB::Array{Cascade.Level,1})` 
        ... adds two sets of levels so that each levels occurs only 'once' in the list; in practice, however, this 'addition' requires also 
            that the parent and daughter processes are added properly so that all information is later available for the simulations.
            It is assumed here that all daugther and parent (processes) appear only once if levels from different data sets
            (Cascade.DecayData, Cascade.PhotoIonData) are added to each other. A message is issued about the number of levels before and 
            after this 'addition', and how many of the levels have been modified by this method. Note that all relative occucations are 
            set to zero in this addition; a newlevels::Array{Cascade.Level,1} is returned.
    """
    function  addLevels(levelsA::Array{Cascade.Level,1}, levelsB::Array{Cascade.Level,1})
        nA = length(levelsA);   nB = length(levelsB);    nmod = 0;    nnew = 0;    newlevels = Cascade.Level[];  appendedB = falses(nB)
        
        # First all levels from levels A but taking additional parents and daugthers into accout
        for  levA in levelsA
            parents   = deepcopy(levA.parents);     daugthers = deepcopy(levA.daugthers);   
            for  (i,levB) in enumerate(levelsB)
                if levA == levB     appendedB[i] = true
                    for p in levB.parents     push!(parents,   p)   end
                    for d in levB.daugthers   push!(daugthers, d)   end
                    if  length(parents) > length(levA.parents)  ||  length(daugthers) > length(levA.daugthers)  nmod = nmod + 1   end
                    break
                end
            end
            push!(newlevels, Cascade.Level(levA.energy, levA.J, levA.parity, levA.NoElectrons, 0., parents, daugthers) )
        end
        
        # Append those levels from levelsB that are not yet appended
        for  (i,levB) in enumerate(levelsB)
            if appendedB[i]
            else    push!(newlevels, levB);     nnew = nnew + 1
            end 
        end
        nN = length(newlevels)
        println("> Append $nnew (new) levels to $nA levels results in a total of $nN levels (with $nmod modified levels) in the list.")
        return( newlevels )
    end
    

    """
    `Cascade.assignOccupation!(levels::Array{Cascade.Level,1}, property::AbstractSimulationProperty)` 
        ... assigns the occupation due to the given property
    """
    function assignOccupation!(levels::Array{Cascade.Level,1}, property::AbstractSimulationProperty)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        
        if  typeof(property) in [Cascade.IonDistribution, Cascade.PhotonIntensities]
            # Assign the relative occupation of the levels in this list due to the given initial occupation
            for  pair  in  property.initialOccupations
                levels[pair[1]].relativeOcc = pair[2]
            end
            println("> Assign an initial occupation for given level numbers.")
            if printSummary  println(iostream, "> Assign an initial occupation for given level numbers.")   end
        end

        return( nothing )
    end
    

    """
    `Cascade.combineEnergiesIntensities(w1::Float64, w1enInts::Array{Tuple{Float64,Float64},1}, 
                                        w2::Float64, w2enInts::Array{Tuple{Float64,Float64},1})` 
        ... combines w1 * w1enInts + w2 * w2enInts; a newEnergiesInts::Array{Tuple{Float64,Float64},1} is returned.
    """
    function combineEnergiesIntensities(w1::Float64, w1enInts::Array{Tuple{Float64,Float64},1}, 
                                        w2::Float64, w2enInts::Array{Tuple{Float64,Float64},1})
        newEnergiesInts = Tuple{Float64,Float64}[]
        for  enInt in w1enInts   push!(newEnergiesInts, (enInt[1], w1*enInt[2]))    end
        for  enInt in w2enInts   push!(newEnergiesInts, (enInt[1], w2*enInt[2]))    end
        
        return( newEnergiesInts )
    end
    
    
    """
    `Cascade.displayExpansionOpacities(stream::IO, sc::String, property::Cascade.ExpansionOpacities, 
                                       energyInterval::Tuple{Float64, Float64}, kappas::Array{Basics.EmProperty,1})` 
        ... displays the expansion opacities in a neat table. Nothing is returned.
    """
    function displayExpansionOpacities(stream::IO, sc::String, property::Cascade.ExpansionOpacities, 
                                       energyInterval::Tuple{Float64, Float64}, kappas::Array{Basics.EmProperty,1})
        nx = 63
        println(stream, " ")
        sa = "  Expansion opacities:  $sc       ... are evaluated for the following parameters: \n" *
            "\n    + level population                   = $(property.levelPopulation)    " *
            "\n    + opacityDependence                  = $(property.opacityDependence)    " *
            "\n    + ion density [ions/cm^3]            = $(property.ionDensity)    " *
            "\n    + plasma temperature [K]             = $(property.temperature)   " *
            "\n    + expansion/observation time [sec]   = $(property.expansionTime) " *
            "\n    + binning [Hartree]                  = $(property.opacityDependence.binning) " *
            "\n    + energy interval [Hartree]          = $(energyInterval) " *
            "\n    + energy shift  [Hartree]            = $(property.transitionEnergyShift) \n"
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        sb = TableStrings.inUnits("energy")
        if  typeof(property.opacityDependence) == Cascade.TemperatureOpacityDependence    sb = "[dim-less]"     end
        sa = "  "
        sa = sa * TableStrings.center(20, "Values " * sb; na=1)        
        sa = sa * TableStrings.center(36, "Cou -- kappa^(expansion) [cm^2/g] -- Bab";      na=2)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        for (i,value) in enumerate(property.dependencyValues)
            sa = "       " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", value)) * 
                 "         " * @sprintf("%.6e", kappas[i].Coulomb) * "        " * @sprintf("%.6e", kappas[i].Babushkin)
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))

        return( nothing )
    end
    
    
    """
    `Cascade.displayIonDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})` 
        ... displays the (current or final) ion distribution in a neat table. Nothing is returned.
    """
    function displayIonDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})
        minElectrons = 1000;   maxElectrons = 0;   totalProb = 0.;    nx = 31
        for  level in levels   minElectrons = min(minElectrons, level.NoElectrons);   maxElectrons = max(maxElectrons, level.NoElectrons)   end
        println(stream, " ")
        println(stream, "  (Final) Ion distribution for the cascade:  $sc ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(14, "No. electrons"; na=4)        
        sa = sa * TableStrings.center(10,"Rel. occ.";      na=2)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        for n = maxElectrons:-1:minElectrons
            sa = "             " * string(n);   sa = sa[end-10:end];   prob = 0.
            for  level in levels    if  n == level.NoElectrons   prob = prob + level.relativeOcc    end    end
            sa = sa * "         " * @sprintf("%.5e", prob)
            println(stream, sa)
            totalProb = totalProb + prob
        end
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  Total distributed probability:  " * @sprintf("%.5e", totalProb)
        println(stream, sa)

        return( nothing )
    end


    """
    `Cascade.displayLevelDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})` 
        ... displays the (current or final) level distribution in a neat table. Only those levels with a non-zero 
            occupation are displayed here. Nothing is returned.
    """
    function displayLevelDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})
        minElectrons = 1000;   maxElectrons = 0;   energies = zeros(length(levels));    nx = 69
        for  i = 1:length(levels)
            minElectrons = min(minElectrons, levels[i].NoElectrons);   maxElectrons = max(maxElectrons, levels[i].NoElectrons)
            energies[i]  = levels[i].energy   
        end
        enIndices = sortperm(energies, rev=true)
        # Now printout the results
        println(stream, " ")
        println(stream, "  (Final) Level distribution for the cascade:  $sc")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(14, "No. electrons"; na=2)        
        sa = sa * TableStrings.center( 8, "Lev-No"; na=2)        
        sa = sa * TableStrings.center( 6, "J^P"          ; na=3);               
        sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(10, "Rel. occ.";                                    na=2)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
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
        println(stream, "  ", TableStrings.hLine(nx))

        return( nothing )
    end


    """
    `Cascade.displayLevelTree(stream::IO, levels::Array{Cascade.Level,1}; extended::Bool=false)` 
        ... displays all defined levels  in a neat table, together with their No. of electrons, symmetry, level energy, 
            current (relative) population as well as analogue information about their parents and daugther levels. This 
            enables one to recognize (and perhaps later add) missing parent and daughter levels. Nothing is returned.
    """
    function displayLevelTree(stream::IO, levels::Array{Cascade.Level,1}; extended::Bool=false)
        minElectrons = 1000;   maxElectrons = 0;   energies = zeros(length(levels));    nx = 179;    ny = 65
        for  i = 1:length(levels)
            minElectrons = min(minElectrons, levels[i].NoElectrons);   maxElectrons = max(maxElectrons, levels[i].NoElectrons)
            energies[i]  = levels[i].energy   
        end
        enIndices = sortperm(energies, rev=true)
        # Now printout the results
        println(stream, " ")
        println(stream, "* Level tree of this cascade:")
        println(stream, " ")
        if  extended    println(stream, "  ", TableStrings.hLine(nx))  else    println(stream, "  ", TableStrings.hLine(ny))  end
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
        if  extended    println(stream, "  ", TableStrings.hLine(nx))  else    println(stream, "  ", TableStrings.hLine(ny))  end
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
                        if      p.process == Basics.Auger()         lev = p.lineSet.linesA[idx].initialLevel
                        elseif  p.process == Basics.Radiative()     lev = p.lineSet.linesR[idx].initialLevel
                        elseif  p.process == Basics.Photo()         lev = p.lineSet.linesP[idx].initialLevel
                        else    error("stop a")    end
                        push!( pProcessSymmetryEnergyList, (p.process, lev.basis.NoElectrons, LevelSymmetry(lev.J, lev.parity), lev.energy) )
                    end
                    for  d in levels[en].daugthers
                        idx = d.index
                        if      d.process == Basics.Auger()         lev = d.lineSet.linesA[idx].finalLevel
                        elseif  d.process == Basics.Radiative()     lev = d.lineSet.linesR[idx].finalLevel
                        elseif  d.process == Basics.Photo()         lev = d.lineSet.linesP[idx].finalLevel
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
        if  extended    println(stream, "  ", TableStrings.hLine(nx))  else    println(stream, "  ", TableStrings.hLine(ny))  end

        return( nothing )
    end
    

    """
    `Cascade.displayIntensities(stream::IO, property::PhotonIntensities, energiesIntensities::Array{Tuple{Float64,Float64},1})` 
        ... displays the (tuples of) energiesIntensities in a neat table. Nothing is returned.
    """
    function displayIntensities(stream::IO, property::PhotonIntensities, energiesIntensities::Array{Tuple{Float64,Float64},1})
        nx = 40
        sMinEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", property.minPhotonEnergy))
        sMaxEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", property.maxPhotonEnergy))
        println(stream, " ")
        println(stream, "* Energies & (relative) photon intensities between " * sMinEn * " and "  * sMaxEn * 
                           TableStrings.inUnits("energy") * ":  ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(16, "Energy "        * TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(10, "Rel. Intensity"; na=2)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        #
        totalIntensity = 0.
        for  enInt in  energiesIntensities
            sa = "     "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", enInt[1])) * "         " * @sprintf("%.3e", enInt[2])
            totalIntensity = totalIntensity + enInt[2]
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "  Total (relative) intensity:  " * @sprintf("%.3e", totalIntensity))

        return( nothing )
    end
    

    """
    `Cascade.displayPhotoAbsorptionSpectrum(stream::IO, crossSections::Array{Cascade.AbsorptionCrossSection,1}, settings::Cascade.SimulationSettings)` 
        ... displays the photoabsorption cross sections a neat table. Nothing is returned.
    """
    function displayPhotoAbsorptionSpectrum(stream::IO, crossSections::Array{Cascade.AbsorptionCrossSection,1}, settings::Cascade.SimulationSettings)
        nx = 83
        println(stream, " ")
        println(stream, "* Absorption cross sections:  ")
        println(stream, " ")
        sMinEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", settings.minPhotonEnergy))
        sMaxEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", settings.maxPhotonEnergy))
        println(stream, "  Absorption cross sections are determined for photon energies between " * sMinEn * " and "  *
                        sMaxEn * TableStrings.inUnits("energy") * " as well as for levels \n  with the initial population " *
                        "$(settings.initialOccupations) \n")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(16, "Energy "        * TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(10, "Ionization CS " * TableStrings.inUnits("cross section"); na=11)
        sa = sa * TableStrings.center(10, "Excitation CS " * TableStrings.inUnits("cross section"); na=11)
        println(stream, sa)
        sa = "                      Coulomb       Babushkin         Coulomb       Babushkin"
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        #
        for  cs in  crossSections
            sa = "     "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", cs.photonEnergy)) * "     "
            if     cs.ionizationCS == Basics.EmProperty(0.)       sa = sa * "                                         "
            else   sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.ionizationCS.Coulomb))   * "   " 
                   sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.ionizationCS.Babushkin)) * "         "
            end
            #
            if     cs.ionizationCS == Basics.EmProperty(0.)       sa = sa * "                                         "
            else   sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.excitationCS.Coulomb))   * "   " 
                   sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.excitationCS.Babushkin)) * "    "
            end
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))

        return( nothing )
    end
    

    """
    `Cascade.displayRelativeOccupation(stream::IO, levels::Array{Cascade.Level,1}, settings::Cascade.SimulationSettings)` 
        ... displays the (initial) relative occupation of the levels in a neat table; an error message is issued if the population is
            given for those levels in the settings, which do not exist in the present simulation. Nothing is returned.
    """
    function displayRelativeOccupation(stream::IO, levels::Array{Cascade.Level,1}, settings::Cascade.SimulationSettings)
       nx = 69
        println(stream, " ")
        println(stream, "* Initial level occupation:  ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(14, "No. electrons"; na=2)        
        sa = sa * TableStrings.center( 8, "Lev-No"; na=2)        
        sa = sa * TableStrings.center( 6, "J^P"          ; na=3);               
        sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(10, "Rel. occ.";                                    na=2)
        println(stream, sa)
        #
        for  initialOcc in  settings.initialOccupations
            idx = initialOcc[1];   occ = initialOcc[2]
            if   idx < 1   ||   idx > length(levels)       error("In appropriate choice of initial occupation; idx = $idx")    end
            level = levels[idx]
            sa = "            " * string(level.NoElectrons);                                  sa  = sa[end-10:end]
            sb = "            " * string(idx);                                                sb  = sb[end-12:end]
            sc = "    " * string( LevelSymmetry(level.J, level.parity) )  * "           ";    sc  = sc[1:15]
            sd = sa * sb * sc  * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", level.energy))
            sd = sd * "      " * @sprintf("%.5e", level.relativeOcc) 
            println(stream, sd)
        end
        println(stream, "  ", TableStrings.hLine(nx))

        return( nothing )
    end
    

    """
    `Cascade.displayRelativeOccupation(stream::IO, levels::Array{Cascade.Level,1})` 
        ... displays the (initial) relative occupation of the levels in a neat table. Nothing is returned.
    """
    function displayRelativeOccupation(stream::IO, levels::Array{Cascade.Level,1})
       nx = 69
        println(stream, " ")
        println(stream, "* Initial level occupation:  ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(14, "No. electrons"; na=2)        
        sa = sa * TableStrings.center( 8, "Lev-No"; na=2)        
        sa = sa * TableStrings.center( 6, "J^P"          ; na=3);               
        sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(10, "Rel. occ.";                                    na=2)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        #
        for  (levi, level) in enumerate(levels)
            sa = "            " * string(level.NoElectrons);                                  sa  = sa[end-10:end]
            sb = "            " * string(levi);                                               sb  = sb[end-12:end]
            sc = "    " * string( LevelSymmetry(level.J, level.parity) )  * "           ";    sc  = sc[1:15]
            sd = sa * sb * sc  * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", level.energy))
            sd = sd * "      " * @sprintf("%.5e", level.relativeOcc) 
            println(stream, sd)
        end
        println(stream, "  ", TableStrings.hLine(nx))

        return( nothing )
    end
    

    """
    `Cascade.extractOccupation(levels::Array{Cascade.Level,1}, groundConfigs::Array{Configuration,1})` 
        ... determines the total occupation of the levels in (one of) the groundConfigs. A occ::Float64 is returned.
    """
    function extractOccupation(levels::Array{Cascade.Level,1}, groundConfigs::Array{Configuration,1})
        #
        wocc = 0.
        for level in levels
            ##x @show level
            ##x conf = Basics.extractLeadingConfiguration(level)
            ##x @show conf, groundConfigs
            if  Basics.extractLeadingConfiguration(level) in groundConfigs   wocc = wocc + level.relativeOcc     end
        end

        return( wocc )
    end
    


    """
    `Cascade.extractPhotoExcitationData(dataDicts::Array{Dict{String,Any},1})` 
        ... returns the available photoexcitation data.
    """
    function extractPhotoExcitationData(dataDicts::Array{Dict{String,Any},1})
        photoexcitationData = Cascade.ExcitationData[]
        for data  in  dataDicts       results = data["results"]
            if  haskey(results, "photoexcitation line data:") 
                push!(photoexcitationData, results["photoexcitation line data:"])   
            end
        end
        
        return( photoexcitationData )
    end
    


    """
    `Cascade.extractPhotoIonizationData(dataDicts::Array{Dict{String,Any},1})` 
        ... returns the available photoionization data.
    """
    function extractPhotoIonizationData(dataDicts::Array{Dict{String,Any},1})
        photoionizationData = Cascade.PhotoIonData[]
        for data  in  dataDicts       results = data["results"]
            if  haskey(results, "photo-ionizing line data:")  photoionizationData = results["photo-ionizing line data:"]    end
        end
        
        return( photoionizationData )
    end


    """
    `Cascade.extractLevels(data::Cascade.DecayData, settings::Cascade.SimulationSettings)` 
        ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
            propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
            with its parent lines (which may populate the level) and the daugther lines (to which the pobability may decay). 
            A levelList::Array{Cascade.Level,1} is returned.
    """
    function extractLevels(data::Cascade.DecayData, settings::Cascade.SimulationSettings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        levels = Cascade.Level[]
        print("> Extract and sort the list of levels for the given decay data ... ")
        if printSummary     print(iostream, "> Extract and sort the list of levels for the given decay data ... ")     end
        
        for  i = 1:length(data.linesR)
            line = data.linesR[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(data, Basics.Radiative(), i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(data, Basics.Radiative(), i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end

        for  i = 1:length(data.linesA)
            line = data.linesA[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(data, Basics.Auger(), i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(data, Basics.Auger(), i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end
        
        # Sort the levels by energy in reversed order
        energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
        enIndices = sortperm(energies, rev=true)
        newlevels = Cascade.Level[]
        for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
        
        println("a total of $(length(newlevels)) levels were found.")
        if printSummary     println(iostream, "a total of $(length(newlevels)) levels were found.")     end
        
        return( newlevels )
    end


    """
    `Cascade.extractLevels(data::Array{Cascade.PhotoIonData,1}, settings::Cascade.SimulationSettings)` 
        ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
            propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
            with its parent lines (which may populate the level) and the daugther lines (to which the pobability may decay). 
            A levelList::Array{Cascade.Level,1} is returned.
    """
    function extractLevels(data::Array{Cascade.PhotoIonData,1}, settings::Cascade.SimulationSettings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream");   found = false
        photonEnergy = settings.initialPhotonEnergy;   newlevels = Cascade.Level[]
        for dataset in data
            println(">>> photonEnergy = $(dataset.photonEnergy) ")
            if  dataset.photonEnergy == photonEnergy       found     = true
                newlevels = Cascade.extractLevels(dataset, settings);    break
            end
        end
        
        if  !found  error("No proper photo-ionizing data set (Cascade.PhotoIonData) found for the photon energy $photonEnergy ")  end
        return( newlevels )
    end


    """
    `Cascade.extractLevels(data::Cascade.ExcitationData, settings::Cascade.SimulationSettings)` 
        ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
            propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
            with its parent lines (which may populate the level) and the daugther lines (to which the pobability may decay). 
            A levelList::Array{Cascade.Level,1} is returned.
    """
    function extractLevels(data::Cascade.ExcitationData, settings::Cascade.SimulationSettings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        levels = Cascade.Level[]
        print("> Extract and sort the list of levels for the given excitation data ... ")
        if printSummary     print(iostream, "> Extract and sort the list of levels for the given excitation data ... ")     end
        
        for  i = 1:length(data.linesE)
            line = data.linesE[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(data, Basics.Radiative(), i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(data, Basics.Radiative(), i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end
        
        # Sort the levels by energy in reversed order
        energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
        enIndices = sortperm(energies, rev=true)
        newlevels = Cascade.Level[]
        for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
        
        println("a total of $(length(newlevels)) levels were found.")
        if printSummary     println(iostream, "a total of $(length(newlevels)) levels were found.")     end
        
        return( newlevels )
    end


    """
    `Cascade.extractLevels(data::Cascade.PhotoIonData, settings::Cascade.SimulationSettings)` 
        ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
            propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
            with its parent lines (which may populate the level) and the daugther lines (to which the pobability may decay). 
            A levelList::Array{Cascade.Level,1} is returned.
    """
    function extractLevels(data::Cascade.PhotoIonData, settings::Cascade.SimulationSettings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        levels = Cascade.Level[]
        print("> Extract and sort the list of levels for the given photo-ionization data ... ")
        if printSummary     print(iostream, "> Extract and sort the list of levels for the given photo-ionization data ... ")     end

        for  i = 1:length(data.linesP)
            line = data.linesP[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(data, Basics.Photo(), i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(data, Basics.Photo(), i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end
        
        # Sort the levels by energy in reversed order
        energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
        enIndices = sortperm(energies, rev=true)
        newlevels = Cascade.Level[]
        for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
        
        println("a total of $(length(newlevels)) levels were found.")
        if printSummary     println(iostream, "a total of $(length(newlevels)) levels were found.")     end
        
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
    `Cascade.generateAbsorptionCrossSections(excitationData::Array{Cascade.ExcitationData,1}, 
                                             ionizationData::Array{Cascade.PhotoIonData,1}, settings::Cascade.SimulationSettings)` 
        ... generates a list of photo-absorption cross sections by cycling through all (photo-) excitationData and ionizationData.
            The photoionization contributions to the absorption cross section is obtained by averaging the photoionization cross sections
            for the occupation of the given initial levels; for the photoionization contributions to the absorption cross section,
            the (relative) cross section is taken at neighboured energies of the excitation line of interest, and the associated 
            ionization contribution is obtained from a linear interpolation of two neighboured photon energies or the (constant)
            extrapolation of just one neighboured photon energies of the previously determined photoionization cross sections.
            An list of crossSections::Array{AbsorptionCrossSection,1} is returned which is ordered by energy.
    """
    function generateAbsorptionCrossSections(excitationData::Array{Cascade.ExcitationData,1}, 
                                             ionizationData::Array{Cascade.PhotoIonData,1}, settings::Cascade.SimulationSettings)
        ionizationCS = Cascade.AbsorptionCrossSection[];    excitationCS = Cascade.AbsorptionCrossSection[]
        # First collect the contributions of the initially occupied levels from the photoionization cross sections, and
        # weigthed by the relative population; data are collected for all available photon energies as they might be needed for
        # interpolation
        for  data  in  ionizationData
            cs = Basics.EmProperty(0.)
            for  (levelIdx, occ)  in  settings.initialOccupations
                for  line  in  data.linesP
                    ##x @show  line.photonEnergy, data.photonEnergy
                    if  abs( (line.photonEnergy - data.photonEnergy)/line.photonEnergy ) > 1.0e-5     error("stop a")      end
                    if  line.initialLevel.index ==  levelIdx            cs = cs + occ * line.crossSection                  end
                end
            end
            push!(ionizationCS, Cascade.AbsorptionCrossSection(data.photonEnergy, Basics.EmProperty(0.), cs))
        end
        #
        # Next determine all excitation cross sections for the requested interval of photon energies; here the photon energies are
        # simply given by the photon excitation energies, and which must be in the requested interval. There will be one
        # cross section for each photoexcitation independent of possible blendings.
        for  data  in  excitationData
            ics = Basics.EmProperty(0.);   ecs = Basics.EmProperty(0.)
            for  (levelIdx, occ)  in  settings.initialOccupations
                for  line  in  data.linesE
                    addEcs = false
                    if  settings.minPhotonEnergy <= line.omega <= settings.maxPhotonEnergy     &&
                        line.initialLevel.index ==  levelIdx            
                        ecs = occ * line.crossSection;      addEcs = true 
                    end
                    if  addEcs      ics = Cascade.interpolateIonizationCS(line.omega, ionizationCS)
                                    push!(excitationCS, Cascade.AbsorptionCrossSection(line.omega, ecs, ics))
                    end
                end
            end
        end
        
        # Append those (pure) ionization cross sections which belong the requested interval
        for  cs in  ionizationCS
            if  settings.minPhotonEnergy <= cs.photonEnergy <= settings.maxPhotonEnergy     push!(excitationCS, cs)     end
        end
        sortedCS = Base.sort( excitationCS , lt=Base.isless)
        return( sortedCS )
    end
    


    """
    `Cascade.interpolateIonizationCS(photonEnergy::Float64, ionizationCS::Array{Cascade.AbsorptionCrossSection,1})` 
        ... interpolates (or extrapolates) the ionization cross sections as defined by ionizationCS for the given photonEnergy.
            If photonEnergy is outside the photon energies from ionizationCS, simply the cross section from the nearest energy
            is returned; if photonEnergy lays between two photon energies from ionizationCS, a simple linear interpolation
            rules is applied here. A cs::Basics.EmProperty is returned.
    """
    function interpolateIonizationCS(photonEnergy::Float64, ionizationCS::Array{Cascade.AbsorptionCrossSection,1})
        imin = imax = 0
        for  (i, ionCS)  in  enumerate(ionizationCS)
           if  ionCS.photonEnergy <= photonEnergy   imin = i    end
        end
        for  (i, ionCS)  in  enumerate(ionizationCS)
           if  ionCS.photonEnergy >  photonEnergy   imax = i;   break    end
        end
        #
        if       imin == 0  &&  imax == 1        return(ionizationCS[1].ionizationCS)
        elseif   imax == 0                       return(ionizationCS[end].ionizationCS)
        elseif   imax - imin == 1
            deltaEnergy = photonEnergy - ionizationCS[imin].photonEnergy
            totalEnergy = ionizationCS[imax].photonEnergy - ionizationCS[imin].photonEnergy  
            csCou       = ionizationCS[imin].ionizationCS.Coulomb + 
                          deltaEnergy/totalEnergy * (ionizationCS[imax].ionizationCS.Coulomb - ionizationCS[imin].ionizationCS.Coulomb)
            csBab       = ionizationCS[imin].ionizationCS.Babushkin + 
                          deltaEnergy/totalEnergy * (ionizationCS[imax].ionizationCS.Babushkin - ionizationCS[imin].ionizationCS.Babushkin)
            return( EmProperty(csCou, csBab) )
        else  error("stop b")    
        end
    end



    """
    `Cascade.perform(simulation::Cascade.Simulation`  
        ... to simulate a cascade decay (and excitation) from the given data. Different computational methods and different properties of 
            the ionic system, such as the ion distribution or final-level distribution can be derived and displayed from these simulations. 
            Of course, the details of these simulations strongly depend on the atomic processes and data that have been generated before by 
            performing a computation::Cascade.Computation. The results of all individual steps are printed to screen but nothing is 
            returned otherwise.

    `Cascade.perform(simulation::Cascade.Simulation; output=true)`   
        ... to perform the same but to return the complete output in a dictionary; the particular output depends on the method and 
            specifications of the cascade but can easily accessed by the keys of this dictionary.
    """
    function perform(simulation::Cascade.Simulation; output::Bool=false)
        if  output    results = Dict{String, Any}()    
            results = Base.merge( results, Dict("name:"         => simulation.name) ) 
            results = Base.merge( results, Dict("property:"     => simulation.property) )
        
        else    results = nothing    end
        
        # First review and display the computation data for this simulation; this enables the reader to check the consistency of data.
        # It also returns the level tree to be used in the simulations
        if      typeof(simulation.property) == Cascade.PhotoAbsorption
                                             # -----------------------
            levels = Cascade.reviewData(simulation, ascendingOrder=false)
            wa     = Cascade.simulatePhotoAbsorptionSpectrum(levels, simulation) 
            #
        elseif  typeof(simulation.property) == Cascade.IonDistribution         &&   simulation.method == Cascade.ProbPropagation()
                                             # -----------------------
            levels = Cascade.reviewData(simulation, ascendingOrder=true)
            wa     = Cascade.simulateIonDistribution(levels, simulation) 
            #
        elseif  typeof(simulation.property) == Cascade.FinalLevelDistribution  &&   simulation.method == Cascade.ProbPropagation()
                                             # ------------------------------
            levels = Cascade.reviewData(simulation, ascendingOrder=true)
            wa     = Cascade.simulateLevelDistribution(levels, simulation) 
            #
        elseif  typeof(simulation.property) == Cascade.PhotonIntensities       &&   simulation.method == Cascade.ProbPropagation()
                                             # -------------------------
            levels = Cascade.reviewData(simulation, ascendingOrder=true)
            wa     = Cascade.simulatePhotonIntensities(levels, simulation) 
            if output   results = Base.merge( results, Dict("energies/intensities:"  => wa) )   end
            #
        elseif  typeof(simulation.property) == Cascade.DrRateCoefficients
                                             # --------------------------
            levels = Cascade.reviewData(simulation, ascendingOrder=true)
            wa     = Cascade.simulateDrRateCoefficients(levels, simulation) 
            #
        elseif  typeof(simulation.property) == Cascade.RrRateCoefficients
                                             # --------------------------
            linesR = simulation.computationData[1]["results"]["photo-recombination line data:"].linesR
            @show typeof(linesR)
            wa     = Cascade.simulateRrRateCoefficients(linesR, simulation) 
            if output   results = Base.merge( results, Dict("alpha^RR:"  => wa) )   end
            #
        elseif  typeof(simulation.property) == Cascade.MeanRelaxationTime
                                             # --------------------------
            levels = Cascade.reviewData(simulation, ascendingOrder=true)
            wa     = Cascade.simulateMeanRelaxationTime(levels, simulation) 
            if output   results = Base.merge( results, Dict("relaxPercentage:"  => wa[1]) )   
                        results = Base.merge( results, Dict("relaxTimes:"       => wa[2]) )   end
            #
        elseif  typeof(simulation.property) == Cascade.ExpansionOpacities
                                             # --------------------------
            photoExcData = Cascade.extractPhotoExcitationData(simulation.computationData)
            wa           = Cascade.simulateExpansionOpacities(photoExcData, simulation.name, simulation.property, printout=true) 
            #
        elseif  typeof(simulation.property) == Cascade.RosselandOpacities
                                             # --------------------------
            photoExcData = Cascade.extractPhotoExcitationData(simulation.computationData)
            wa           = Cascade.simulateRosselandOpacities(photoExcData, simulation) 
            #
        else         error("stop b")
        end

        return( results )
    end
    

    """
    `Cascade.propagateOccupationInTime!(levels::Array{Cascade.Level,1}, dt::Float64)` 
        ... propagates the occupation of the levels by dt in time. 
    """
    function propagateOccupationInTime!(levels::Array{Cascade.Level,1}, dt::Float64)
        #
        relativeOcc = zeros(length(levels));    relativeLoss = zeros(length(levels));    down = 0.0;    loss = 0.0
        for (i, level) in  enumerate(levels)
            # Cycle through all daugthers of level and 'shift' the part of occupation down the daugther levels
            lossFactor = 0.;    occ = level.relativeOcc
            for  daugther in level.daugthers
                idx = daugther.index
                if      daugther.process == Basics.Radiative()     line = daugther.lineSet.linesR[idx];  rate = line.photonRate.Coulomb
                elseif  daugther.process == Basics.Auger()         line = daugther.lineSet.linesA[idx];  rate = line.totalRate
                else    error("stop b; process = $(daugther.process) ")
                end
                downFactor = (1.0 - exp(-rate*dt));         
                ##x @show downFactor, occ
                newLevel   = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, 
                                            line.finalLevel.basis.NoElectrons, 0., Cascade.LineIndex[], Cascade.LineIndex[] )
                kk         = Cascade.findLevelIndex(newLevel, levels)
                if lossFactor + downFactor < 1.0    lossFactor = lossFactor + downFactor
                else                                downFactor = 1.0 - lossFactor;      lossFactor = 1.0
                end
                relativeOcc[kk] = relativeOcc[kk] + downFactor * occ;   down = down + downFactor * occ
            end
            ##x @show lossFactor
            levels[i].relativeOcc = levels[i].relativeOcc - lossFactor * occ;   loss = loss + lossFactor * occ
        end
        for i = 1:length(levels)  levels[i].relativeOcc = levels[i].relativeOcc + relativeOcc[i]    end
        ##x @show down, loss
        ##x @show relativeOcc
        
        return( nothing )
    end


    """
    `Cascade.propagateProbability!(levels::Array{Cascade.Level,1}; 
                                   collectPhotonIntensities::Bool=false, collectElectronIntensities::Bool=false)` 
        ... propagates the relative level occupation through the levels of the cascade until no further change occur in the 
            relative level occupation. The argument levels is modified during the propagation, but nothing is returned otherwise.
    """
    function propagateProbability!(levels::Array{Cascade.Level,1}; 
                                   collectPhotonIntensities::Bool=false, collectElectronIntensities::Bool=false)
        # Inititalize arrays to be returned
        if      collectPhotonIntensities  && collectElectronIntensities      error("stop a")
        elseif  collectPhotonIntensities  || collectElectronIntensities
                energiesIntensities = Tuple{Float64,Float64}[]
        end
        
        n = 0
        println("\n*  Probability propagation through $(length(levels)) levels of the cascade:")
        while true
            n = n + 1;    totalProbability = 0.
            print("    $n-th round ... ")
            relativeOcc = zeros(length(levels));    relativeLoss = zeros(length(levels))
            for  level in levels
                ##x @show level.relativeOcc, length(level.daugthers)
                if   level.relativeOcc > 0.   && length(level.daugthers) > 0
                    # A level with relative occupation > 0 has still 'daugther' levels; collect all excitation/decay rates for this level
                    # Here, an excitation cross section is formally treated as a rate as it is assumed that the initial levels of
                    # the photoionization process cannot decay by photon emission or autoionization processes.
                    prob  = level.relativeOcc;   rates = zeros(length(level.daugthers))
                    for  (i,daugther) in  enumerate(level.daugthers)
                        idx = daugther.index
                        if      daugther.process == Basics.Radiative()     rates[i] = daugther.lineSet.linesR[idx].photonRate.Coulomb
                        elseif  daugther.process == Basics.Auger()         rates[i] = daugther.lineSet.linesA[idx].totalRate
                        elseif  daugther.process == Basics.Photo()         rates[i] = daugther.lineSet.linesP[idx].crossSection.Coulomb
                        else    error("stop a; process = $(daugther.process) ")
                        end
                        ##x @show daugther.process, rates[i]
                    end
                    totalRate = sum(rates)
                    if      totalRate <  0.    error("stop b")
                    elseif  totalRate == 0.    # do nothing
                    else
                        # Shift the relative occupation to the 'daugther' levels due to the different ionization and decay pathes
                        for  (i,daugther) in  enumerate(level.daugthers)
                            idx = daugther.index
                            if      daugther.process == Basics.Radiative()     line = daugther.lineSet.linesR[idx]
                            elseif  daugther.process == Basics.Auger()         line = daugther.lineSet.linesA[idx]
                            elseif  daugther.process == Basics.Photo()         line = daugther.lineSet.linesP[idx]
                            else    error("stop b; process = $(daugther.process) ")
                            end
                            newLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, 
                                                      line.finalLevel.basis.NoElectrons, 0., Cascade.LineIndex[], Cascade.LineIndex[] )
                            kk    = Cascade.findLevelIndex(newLevel, levels)
                            relativeOcc[kk] = relativeOcc[kk] + prob * rates[i] / totalRate
                            #
                            if  collectPhotonIntensities   && daugther.process == Basics.Radiative()
                                push!(energiesIntensities, (level.energy - levels[kk].energy, prob * rates[i] / totalRate))     end
                            if  collectElectronIntensities && daugther.process == Basics.Auger()
                                push!(energiesIntensities, (level.energy - levels[kk].energy, prob * rates[i] / totalRate))     end
                        end
                        level.relativeOcc = 0.;   totalProbability = totalProbability + prob
                    end
                end
            end
            for i = 1:length(levels)  levels[i].relativeOcc = levels[i].relativeOcc + relativeOcc[i]    end
            println("has propagated a total of $totalProbability level occupation.")
            # Cycle once more if the relative occupation has still changed
            if  totalProbability == 0.    break    end
        end

        if  collectPhotonIntensities   ||   collectElectronIntensities   return(energiesIntensities)
        else                                                             return( nothing )
        end
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
        return( nothing )
    end


    """
    `Cascade.reviewData(simulation::Cascade.Simulation; ascendingOrder::Bool=false)` 
        ... reviews and displays the (computation) data for the given simulation; these data contains the name of the data set, 
            its initial and generated multiplets for the various blocks of (some part of the ionization and/or decay) cascade as 
            well as all the line data [lineR, linesA, lineP, ...]. From these data, this function also generates and returns
            the level tree that is to be used in the subsequent simulations, and where levels are odered in `ascending' order
            if selected.
    """
    function reviewData(simulation::Cascade.Simulation; ascendingOrder::Bool=false)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        dataDicts = simulation.computationData;     settings = simulation.settings;     allLevels = Cascade.Level[]
        
        # Loop through all (computation) data set and display the major results
        for  (i,data) in  enumerate(dataDicts)
            results     = data["results"]
            multiplets  = results["initial multiplets:"]
            gMultiplets = results["generated multiplets:"]
            nlev = 0;    for multiplet in multiplets     nlev  = nlev  + length(multiplet.levels)     end
            nglev = 0;   for multiplet in gMultiplets    nglev = nglev + length(multiplet.levels)     end
            println("\n* $i) Data dictionary for cascade computation:   $(results["name"])  with  $nlev initial and  $nglev generated levels") 
            println(  "  ===========================================")
            
            Cascade.displayLevels(stdout, multiplets, sa="initial ")
            if  printSummary 
                println(iostream, "\n* $i) Data dictionary for cascade computation:   $(results["name"])  with  $nlev initial and  $nglev generated levels") 
                println(iostream,   "  ===========================================")
                Cascade.displayLevels(iostream, multiplets,  sa="initial ")
                Cascade.displayLevels(iostream, gMultiplets, sa="generated ")        
            end
            #
            if      haskey(results,"decay line data:")                          lineData = results["decay line data:"]
            elseif  haskey(results,"photo-ionizing line data:")                 lineData = results["photo-ionizing line data:"]
            elseif  haskey(results,"photo-excited line data:")                  lineData = results["photo-excited line data:"]
            elseif  haskey(results,"hollow-ion line data:")                     lineData = results["hollow-ion line data:"]
            else    error("stop a")
            end
            levels = Cascade.extractLevels(lineData, settings)
            allLevels = Cascade.addLevels(allLevels, levels)
        end
        
        allLevels = Cascade.sortByEnergy(allLevels, ascendingOrder=ascendingOrder)
        Cascade.assignOccupation!(allLevels, simulation.property)
        if  simulation.settings.printTree      Cascade.displayLevelTree(stdout, allLevels, extended=false)     end
        if  simulation.settings.printLongTree  Cascade.displayLevelTree(stdout, allLevels, extended=true)      end
        
        return( allLevels )
    end
    

    """
    `Cascade.simulateDrRateCoefficients(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
        ... Determines and prints the DR resonance strength and (plasma) rate coefficients for all resonance levels.
            Nothing is returned.
    """
    function simulateDrRateCoefficients(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        resonances = Dielectronic.Resonance[]
        @show length(levels)
        #
        # Collect the information about all resonances
        for  level in levels
            for daugther in level.daugthers
                if  daugther.process != Basics.Auger();                            continue   end
                aLine   = daugther.lineSet.linesA[daugther.index]
                if  aLine.finalLevel.index != simulation.property.initialLevelNo   continue   end
                #
                dJ      = aLine.initialLevel.J;         iJ      = aLine.finalLevel.J
                dM      = aLine.initialLevel.M;         iM      = aLine.finalLevel.M
                dParity = aLine.initialLevel.parity;    iParity = aLine.finalLevel.parity
                dIndex  = aLine.initialLevel.index;     iIndex  = aLine.finalLevel.index
                dEnergy = aLine.initialLevel.energy;    iEnergy = aLine.finalLevel.energy
                captureRate = aLine.totalRate
                augerRate   = Cascade.computeTotalAugerRate(level)
                ##x @show level.energy, captureRate, augerRate
                photonRate  = Cascade.computeTotalPhotonRate(level)
                strength    = Basics.EmProperty(0.)
                #
                iLevel      = ManyElectron.Level(iJ, iM, iParity, iIndex, iEnergy, 0., false, ManyElectron.Basis(), Float64[])
                dLevel      = ManyElectron.Level(dJ, dM, dParity, dIndex, dEnergy, 0., false, ManyElectron.Basis(), Float64[])
                es          = Defaults.convertUnits("energy: to atomic", simulation.property.electronEnergyShift)
                en          = dLevel.energy-iLevel.energy + es;    if  en < 0.  continue   end
                #
                wa          = Defaults.convertUnits("kinetic energy to wave number: atomic units", en)
                wa          = pi*pi / (wa*wa) * captureRate  * 2 * # factor 2 is not really clear.
                              ((Basics.twice(dJ) + 1) / (Basics.twice(iJ) + 1)) /
                              (augerRate + photonRate.Babushkin)
                strength     = EmProperty(wa * photonRate.Coulomb, wa * photonRate.Babushkin)
                newResonance = Dielectronic.Resonance(iLevel, dLevel, en, strength, 0., augerRate, photonRate)
                ## newResonance = Dielectronic.Resonance(iLevel, dLevel, en, strength, captureRate, augerRate, photonRate)
                push!(resonances, newResonance)
                @show augerRate, photonRate, strength
            end
        end 
        
        # Printout the resonance strength
        settings = Dielectronic.Settings(Dielectronic.Settings(), calcRateAlpha=true, temperatures=simulation.property.temperatures)
        Dielectronic.displayResults(stdout, resonances, settings)
        Dielectronic.displayRateCoefficients(stdout, resonances, settings)

        return( nothing )
    end
    

    """
    `Cascade.simulateExpansionOpacities(photoexcitationData::Array{Cascade.ExcitationData,1}, name::String, 
                                        property::Cascade.ExpansionOpacities; printout::Bool=true)` 
        ... runs through all excitation lines, sums up their contributions and form a (list of) expansion opacities for the given 
            parameters. Nothing is returned.
    """
    function simulateExpansionOpacities(photoexcitationData::Array{Cascade.ExcitationData,1}, name::String, 
                                        property::Cascade.ExpansionOpacities; printout::Bool=true)
        #
        function lambda_over_dlambda(opacityDependence::AbstractOpacityDependence, lineOmega::Float64, kT::Float64, depValue::Float64)
            # Calculates the value lambda / Delta lambda for the given binning and depValue, for which the opacity need to be
            # determined; the binning is assumed in Hartree (for frequency- and temperature-normalized dependence) and in
            # nm for wavelength-dependent opacities; an value = 0. is returned if the lineOmega does not fall into the (binning)
            # interval
            wa = 0.;  halfBinning = opacityDependence.binning / 2.
            #
            if     typeof(opacityDependence) == FrequencyOpacityDependence
                # Binning, depValue and lineOmega are all in Hartree and readily to compare; ratio need to be inverted, when compared
                # with wavelength.
                if  depValue - halfBinning < lineOmega < depValue + halfBinning             wa = 2* halfBinning / lineOmega        end
            elseif typeof(opacityDependence) == WavelengthOpacityDependence
                # Binning in nm, lineOmega & depValue in Hartree are first converted into nm ... and ratio is determined in nm as well
                lineOmega_nm = convertUnits("energy: from atomic to Angstrom", lineOmega) / 10.
                depValue_nm  = convertUnits("energy: from atomic to Angstrom", depValue)  / 10.
                if  depValue_nm - halfBinning < lineOmega_nm < depValue_nm + halfBinning    wa = lineOmega_nm / (2* halfBinning)   end
            elseif typeof(opacityDependence) == TemperatureOpacityDependence
                # Binning and lineOmega are in Hartree, depValue in [u] and need to be converted; 
                # ratio still need to be inverted, when compared with wavelength.
                if  depValue * kT - halfBinning < lineOmega < depValue + halfBinning * kT   wa = 2* halfBinning / lineOmega        end
            else   error("stop a")
            end
            
            return( wa )
        end
        #
        #
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        values     = property.dependencyValues;                
        exptime    = property.expansionTime;                   
        dependence = property.opacityDependence;                   
        exptime_au = exptime / convertUnits("time: from atomic to sec", 1.0)
        T          = property.temperature;                                
        rho        = property.ionDensity
        eshift     = property.transitionEnergyShift
        ne         = 1.0 ## number density [1/a_o^3] ?? 
        NoValues   = length(values);                                      
        kappas     = Basics.EmProperty[];    for  i = 1:NoValues     push!(kappas, Basics.EmProperty(0.))   end
        #
        # Determine c in cm/s
        alpha   = Defaults.getDefaults("alpha")
        c_in_SI = Defaults.getDefaults("speed of light: c") * convertUnits("length: from atomic to fm", 1.0) * 1.0e-13 /
                  convertUnits("time: from atomic to sec", 1.0)
        factor  = 1.0 / (rho * c_in_SI * exptime)
        kT      = convertUnits("temperature: from Kelvin to (Hartree) units", T)
        A_au    = convertUnits("length: from fm to atomic", 1.0e5)
        #
        if length(photoexcitationData) == 0     error("No photoexcitationData provided.")     end
        #
        minEnergy = 1000.;   maxEnergy = 0.
        for  excData  in photoexcitationData
            for  line  in excData.linesE
                if  minEnergy > line.omega     minEnergy = line.omega   end
                if  maxEnergy < line.omega     maxEnergy = line.omega   end
                    omega       = line.omega + eshift
                    fosc        = line.oscStrength
                    g0          = Basics.twice(line.initialLevel.J) + 1;   ge = Basics.twice(line.finalLevel.J) + 1
                    lambda_au   = convertUnits("energy: from atomic to Angstrom", omega) * A_au
                for  ivalue = 1:NoValues
                    lmd_over_dl    = lambda_over_dlambda(dependence, omega, kT, values[ivalue])
                    if  lmd_over_dl == 0.  continue  end
                    tau_Cou        = pi * alpha * ne * lambda_au * exptime_au * ge / g0 * fosc.Coulomb   * exp(-omega/kT)
                    tau_Bab        = pi * alpha * ne * lambda_au * exptime_au * ge / g0 * fosc.Babushkin * exp(-omega/kT)
                    term_Cou       = factor * lmd_over_dl * (1.0 - exp(-tau_Cou))
                    term_Bab       = factor * lmd_over_dl * (1.0 - exp(-tau_Bab))
                    kappas[ivalue] = kappas[ivalue] + Basics.EmProperty(term_Cou, term_Bab)
                end
            end
        end
        #
        if  printout
            Cascade.displayExpansionOpacities(stdout, name, property, (minEnergy,maxEnergy), kappas)     
            if  printSummary   
                Cascade.displayExpansionOpacities(iostream, name, property, (minEnergy,maxEnergy), kappas)      end
                
              return( nothing )
        else  return( kappas )
        end

    end
    

    """
    `Cascade.simulateIonDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
        ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
            propagation, it runs through all levels and propagates the probabilty until no level probability changes anymore. The final 
            level distribution is then used to derive the ion distribution or the level distribution, if appropriate. Nothing is returned.
    """
    function simulateIonDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        # Specify and display the initial (relative) occupation
        if length(simulation.property.initialOccupations) > 0
                Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations) 
        else    Cascade.specifyInitialOccupation!(levels, simulation.property.leadingConfigs) 
        end
        Cascade.displayRelativeOccupation(stdout, levels)
        #
        Cascade.propagateProbability!(levels, collectPhotonIntensities=false)  
        #
        ionDist = Cascade.displayIonDistribution(stdout, simulation.name, levels)     
        if  printSummary   Cascade.displayIonDistribution(iostream, simulation.name, levels)      end

        return( ionDist )
    end
    

    """
    `Cascade.simulateLevelDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
        ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
            propagation, it runs through all levels and propagates the probabilty until no level probability changes anymore. The final 
            level distribution is then used to derive the ion distribution or the level distribution, if appropriate. Nothing is returned.
    """
    function simulateLevelDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        Cascade.propagateProbability!(levels)   ##x , simulation.cascadeData)
        #
        if  typeof(simulation.property) == Cascade.IonDistribution   
            Cascade.displayIonDistribution(stdout, simulation.name, levels)     
            if  printSummary   Cascade.displayIonDistribution(iostream, simulation.name, levels)      end
        end
        if  typeof(simulation.property) == Cascade.FinalLevelDistribution    
            Cascade.displayLevelDistribution(stdout, simulation.name, levels)   
            if  printSummary   Cascade.displayLevelDistribution(iostream, simulation.name, levels)    end
        end

        return( nothing )
    end
    

    """
    `Cascade.simulateMeanRelaxationTime(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
        ... determine the mean relaxation time until 70%, 80%, 90%  of the initially occupied levels decay down to
            the ground configurations. An relaxTimes::Array{Float64,1} is returned that contains the
            mean relaxation times for 70%, 80%, 90%, ...
    """
    function simulateMeanRelaxationTime(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        relaxTimes = zeros(3);     relaxPercentage = [0.7, 0.8, 0.9];
        time = 0.0;     dt = simulation.property.timeStep;   nx = 0
        # Specify and display the initial (relative) occupation
        if length(simulation.property.initialOccupations) > 0
                Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations) 
        else    Cascade.specifyInitialOccupation!(levels, simulation.property.leadingConfigs) 
        end
        Cascade.displayRelativeOccupation(stdout, levels)
        # Determine the smallest and largest rate for the given cascade tree
        minRate = 1.0e100;    maxRate = 0.0
        for  level in levels
            for daugther in level.daugthers
                if      daugther.process == Basics.Auger()
                    aLine  = daugther.lineSet.linesA[daugther.index]
                    @show aLine.totalRate
                    if  minRate > aLine.totalRate > 0.  minRate = aLine.totalRate   end
                    if  maxRate < aLine.totalRate       maxRate = aLine.totalRate   end
                elseif  daugther.process == Basics.Radiative()
                    rLine  = daugther.lineSet.linesR[daugther.index]
                    @show rLine.photonRate.Coulomb
                    if  minRate > rLine.photonRate.Coulomb > 0.  minRate = rLine.photonRate.Coulomb   end
                    if  maxRate < rLine.photonRate.Coulomb       maxRate = rLine.photonRate.Coulomb   end
                else    error("stop a")
                end
            end
        end
        println(">> Simulate cascade with time step $(dt) for minRate = $minRate and  maxRate = $maxRate; " *
                "minRate*timeStep = $(minRate*dt) ")
        wocc = Cascade.extractOccupation(levels, simulation.property.groundConfigs)
        @show wocc
        #
        # Now propagate the occupation in time
        goon = true;    pr70 = true;    pr80 = true;    pr90 = true
        while  goon
            time = time + dt;   nx = nx + 1
            Cascade.propagateOccupationInTime!(levels, dt)
            wocc = Cascade.extractOccupation(levels, simulation.property.groundConfigs)
            if  rem(nx,10000) == 0     println(">>> $nx) time = $time [a.u.]   ground conf. occupation = $wocc")   end
            for  i = 1:length(relaxPercentage) 
                if wocc < relaxPercentage[i]    relaxTimes[i] = time    end
            end
            if  pr70  &&  wocc > 0.70   pr70 = false;   println("70%  at  time = $(relaxTimes[1])")   end
            if  pr80  &&  wocc > 0.80   pr80 = false;   println("80%  at  time = $(relaxTimes[2])")   end
            if  pr90  &&  wocc > 0.90   pr90 = false;   println("90%  at  time = $(relaxTimes[3])")   end
            if  wocc > 0.91     goon = false    end
        end

        return( relaxPercentage, relaxTimes )
    end
    

    """
    `Cascade.simulatePhotoAbsorptionSpectrum(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
        ... cycle through all levels and (incident photon) energies to derive the overall photo-absorption spectrum.
            All absorption cross sections are displayed in a neat table but nothing is returned otherwise.
    """
    function simulatePhotoAbsorptionSpectrum(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        crossSections = Cascade.AbsorptionCrossSection[]
        # Display the initial (relative) occupation and the photon energies for which photo-ionization cross sections are available
        Cascade.displayRelativeOccupation(stdout, levels, simulation.settings)
        excitationData = Cascade.extractPhotoExcitationData(simulation.computationData)
        ionizationData = Cascade.extractPhotoIonizationData(simulation.computationData)
        photonEnergies = Float64[];   for data  in  ionizationData   push!(photonEnergies, data.photonEnergy)   end
        println("\n* Photo-ionization data are defined for the photon energies [in a.u.]: \n   $photonEnergies")
        crossSections = Cascade.generateAbsorptionCrossSections(excitationData, ionizationData, simulation.settings)
        #
        Cascade.displayPhotoAbsorptionSpectrum(stdout, crossSections, simulation.settings)
        if  printSummary   Cascade.displayPhotoAbsorptionSpectrum(iostream, crossSections, simulation.settings)    end

        return( nothing )
    end
    

    """
    `Cascade.simulatePhotonIntensities(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
        ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
            propagation, it runs through all levels and propagates the probabilty until no level probability changes anymore. The final 
            level distribution is then used to derive the ion distribution or the level distribution, if appropriate. Nothing is returned.
    """
    function simulatePhotonIntensities(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        # Specify and display the initial (relative) occupation
        if length(simulation.property.initialOccupations) > 0
                Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations) 
        else    Cascade.specifyInitialOccupation!(levels, simulation.property.leadingConfigs) 
        end
        Cascade.displayRelativeOccupation(stdout, levels)
        #
        prop = simulation.property
        energiesInts = Cascade.propagateProbability!(levels, collectPhotonIntensities=true)  
        energiesInts = Cascade.truncateEnergiesIntensities(energiesInts, prop.minPhotonEnergy, prop.maxPhotonEnergy)
        #
        Cascade.displayIntensities(stdout, simulation.property, energiesInts)
        if  printSummary   Cascade.displayIntensities(iostream, simulation.property, energiesInts)      end

        return( energiesInts )
    end
    

    """
    `Cascade.simulateRosselandOpacities(photoexcitationData::Array{Cascade.ExcitationData,1}, simulation::Cascade.Simulation)` 
        ... runs through all excitation lines, sums up their contributions and form a (list of) Rosseland opacities, based on
            the expansion opacities, for the given parameters. Nothing is returned.
    """
    function simulateRosselandOpacities(photoexcitationData::Array{Cascade.ExcitationData,1}, simulation::Cascade.Simulation)
        ulist, wlist = FastGaussQuadrature.gausslaguerre(8);            kappaList    = Float64[]
        
        opacityDependence = simulation.property.opacityDependence
        for  rho in simulation.property.ionDensities
            for  T in simulation.property.temperatures
                property=Cascade.ExpansionOpacities(Basics.BoltzmannLevelPopulation(), opacityDependence, rho, T, 
                                 simulation.property.expansionTime, simulation.property.transitionEnergyShift, ulist)
                                 
                kappas = Cascade.simulateExpansionOpacities(photoexcitationData, "expansion opacity for rho = $rho & T=$T", 
                                                            property, printout=false)    
                #
                # Form the u-integral of the Rosseland opacity
                rosseland = Basics.EmProperty(0.)
                for (i,u)  in  enumerate(ulist)
                    rosseland = rosseland + 15.0/ (4*pi^4) * u^4 * wlist[i] / ( (1.0 - exp(-u))^2 ) * kappas[i]
                end
                #
                sa = "> Rosseland opacity for rho = " * @sprintf("%.3e",rho) * "[g/cm^3]  &  T="    * @sprintf("%.3e",T) *
                     " [K]  is  kappa^Rosseland [cm^2/g] = " * @sprintf("%.5e",rosseland.Coulomb)   * " [Coulomb]  "  *
                                                               @sprintf("%.5e",rosseland.Babushkin) * " [Babushkin]"
                println(sa)
            end
        end
        
        return( nothing )
    end
    

    """
    `Cascade.simulateRrRateCoefficients(lines::Array{PhotoRecombination.Line,1}, simulation)` 
        ... Integrates over all selected cross sections in order to determine the RR plasma rate coefficients for a Maxwellian
            distribution.
    """
    function simulateRrRateCoefficients(lines::Array{PhotoRecombination.Line,1}, simulation)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        # Check consistency of data
        isym = LevelSymmetry(lines[1].initialLevel.J, lines[1].initialLevel.parity)
        for  line  in  lines
             if  LevelSymmetry(line.initialLevel.J, line.initialLevel.parity) != isym   error("error a")    end
        end
        
        # Convert units into cm^3 / s
        factor  = Defaults.convertUnits("length: from atomic to fm", 1.0)^3 * 1.0e-39 * 
                  Defaults.convertUnits("rate: from atomic to 1/s", 1.0) 
        
        temperatures = simulation.property.temperatures
        alphaRR      = EmProperty[]
        #
        for temp  in  temperatures
            temp_au  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
            wa       = EmProperty(0., 0.)
            #
            ## @warn "Cross sections not yet properly set."
            for  line  in  lines
                # Determine cross section of this line
                cs  = EmProperty(0., 0.)
                if   simulation.property.finalLevelSelection.active  &&  
                     !(line.finalLevel.index  in  simulation.property.finalLevelSelection.indices)    continue    end
                pcs = PhotoRecombination.computeCrossSectionForMultipoles(simulation.property.multipoles, line.channels)
                cs  = cs + pcs
                #==
                # Determine cross section if only one final level contributes; for test purposes
                if   line.finalLevel.index == 1
                      wb = PhotoRecombination.crossSectionKramers(line.electronEnergy, 26.0::Float64, (1,50))
                      ## wb = PhotoRecombination.crossSectionStobbe(line.electronEnergy, 26.0::Float64)
                      cs = EmProperty(wb, wb)
                else  cs = EmProperty(0., 0.)
                end  ==#
                wa = wa + 2 / sqrt(2pi) / temp_au^(3/2) * factor * line.electronEnergy * 
                          exp(-line.electronEnergy / temp_au ) * line.weight * cs
            end
            push!(alphaRR, wa)
        end
        
        # Display the RR plasma rate coefficients
        PhotoRecombination.displayRateCoefficients(stdout, isym, simulation.property.temperatures, alphaRR)

        return( alphaRR )
    end
    

    """
    `Cascade.sortByEnergy(levels::Array{Cascade.Level,1}; ascendingOrder::Bool=false)` 
        ... sorts all levels by energy and assigns the occupation as given by the simulation
    """
    function sortByEnergy(levels::Array{Cascade.Level,1}; ascendingOrder::Bool=false)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        newlevels = Cascade.Level[]
        #
        
        # Sort the levels by energy in reversed order
        energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
        if     ascendingOrder   enIndices = sortperm(energies, rev=true)
        else                    enIndices = sortperm(energies, rev=false)
        end
        #
        newlevels = Cascade.Level[]
        for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
        
        ##x # Assign the relative occupation of the levels in this list due to the given settings
        ##x for  pair  in  settings.initialOccupations
        ##x     newlevels[pair[1]].relativeOcc = pair[2]
        ##x end

        println("> Sort a total of $(length(newlevels)) levels, and to which all level numbers refer below. " *
                "Here all charged states are considered together in the overall cascade.")
        if printSummary     
            println(iostream, "> Sort a total of $(length(newlevels)) levels, and to which all level numbers refer below. " *
                              "Here all charged states are considered together in the overall cascade.")
        end

        return( newlevels )
    end
    

    """
    `Cascade.specifyInitialOccupation!(levels::Array{Cascade.Level,1}, initialOccupations::Array{Tuple{Int64,Float64},1})` 
        ... specifies the initial occupation of levels for the given relOccupation; it modifies the occupation 
            but returns nothing otherwise.
    """
    function specifyInitialOccupation!(levels::Array{Cascade.Level,1}, initialOccupations::Array{Tuple{Int64,Float64},1})
        #
        for  initialOcc in  initialOccupations
            idx = initialOcc[1];   occ = initialOcc[2]
            if   idx < 1   ||   idx > length(levels)       error("In appropriate choice of initial occupation; idx = $idx")    end
            levels[idx].relativeOcc = occ
        end

        return( nothing )
    end
    

    """
    `Cascade.specifyInitialOccupation!(levels::Array{Cascade.Level,1}, leadingConfigs::Array{Configuration,1})` 
        ... specifies the initial occupation of levels for the given leadingConfigs; it modifies the occupation 
            but returns nothing otherwise.
    """
    function specifyInitialOccupation!(levels::Array{Cascade.Level,1}, leadingConfigs::Array{Configuration,1})
        #
        nx = 0
        for  lev = 1:length(levels)
            if  Basics.extractLeadingConfiguration(levels[lev])  in  leadingConfigs   nx = nx + 1    end
        end
        if  nx == 0     error("Inappropriate selection of leading configurations for the given set of cascade levels.")     end
        #
        # Now distribute the occupation
        for  lev = 1:length(levels)
            if  Basics.extractLeadingConfiguration(levels[lev]) in leadingConfigs   levels[lev].relativeOcc = 1/nx    end
        end

        return( nothing )
    end
    

    """
    `Cascade.truncateEnergiesIntensities(energiesInts::Array{Tuple{Float64,Float64},1}, minPhotonEnergy::Float64, maxPhotonEnergy::Float64)` 
        ... reduces and truncates the energies & intensities  energiesInts; 'reduce' hereby refer to omit all intensity < 1.0e-8,
            while 'truncate' omits all energies outside of the interval [minPhotonEnergy, miaxPhotonEnergy]. An 
            newEnergiesInts::Array{Tuple{Float64,Float64},1} is returned.
    """
    function truncateEnergiesIntensities(energiesInts::Array{Tuple{Float64,Float64},1}, minPhotonEnergy::Float64, maxPhotonEnergy::Float64)
        # Firt, truncate contributions to given energy range
        @show minPhotonEnergy, maxPhotonEnergy
        w1EnergiesInts = Tuple{Float64,Float64}[];   we = Float64[]
        for  enInt in energiesInts
            if  minPhotonEnergy <= enInt[1] <= maxPhotonEnergy    push!(w1EnergiesInts, enInt);     push!(we, enInt[1])   end
        end
        ##x @show we
        # Add contributions with equal energies contributions
        w2EnergiesInts = Tuple{Float64,Float64}[];  wasConsidered = falses(length(we))
        for  (en, enInt) in  enumerate(w1EnergiesInts)
            if      wasConsidered[en]   continue
            else    
                wa = findall(isequal(enInt[1]), we);    tInt = 0. 
                ##x @show en, wa
                for  a in wa
                    tInt = tInt + w1EnergiesInts[a][2];     wasConsidered[a] = true
                end
                push!(w2EnergiesInts, (enInt[1], tInt))
            end
        end
        ##x @show wasConsidered
        # Finally, omit all small contributions
        newEnergiesInts = Tuple{Float64,Float64}[]
        for  enInt in w2EnergiesInts
            if enInt[2] > 1.0e-8    push!(newEnergiesInts, enInt)   end
        end
        
        return( newEnergiesInts )
    end
