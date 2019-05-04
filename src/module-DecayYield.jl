
"""
`module  JAC.DecayYield`  
    ... a submodel of JAC that contains all methods for computing fluorescence and Auger yields for some level(s).
"""
module DecayYield

    using Printf, JAC, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..TableStrings

    """
    `struct  DecayYield.Settings`  ... defines a type for the details and parameters of computing fluorescence and Auger yields.

        + approach                 ::String           ... Identifier to determine the cascade approach that is applied for the 
                                                          yield computation: {AverageSCA, SCA}
        + printBeforeComputation   ::Bool             ... True if a list of selected levels is to be printed before the 
                                                          actual computations start. 
        + selectLevels             ::Bool             ... True if individual levels are selected for the computation.
        + selectedLevels           ::Array{Level,1}   ... List of selected levels.
    """
    struct Settings 
        approach                   ::String 
        printBeforeComputation     ::Bool
        selectLevels               ::Bool
        selectedLevels             ::Array{Level,1}
    end 


    """
    `DecayYield.Settings()`  
        ... constructor for an `empty` instance of DecayYield.Settings for the computation of fluorescence and Auger yields.
    """
    function Settings()
        Settings("not defined", false, false, Level[])
    end


    # `Base.show(io::IO, settings::DecayYield.Settings)`  ... prepares a proper printout of the variable settings::DecayYield.Settings.
    function Base.show(io::IO, settings::DecayYield.Settings) 
        println(io, "approach:                 $(settings.approach)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
    end


    """
    `struct  DecayYield.Outcome`  
        ... defines a type to keep the outcome of a fluorescence and Auger yield computation as well other results.

        + level                     ::Level              ... Atomic level to which the outcome refers to.
        + NoPhotonLines             ::Int64              ... Number of photon decay lines
        + NoAugerLines              ::Int64              ... Number of Auger decay lines
        + rateR                     ::EmProperty         ... Total fluorescence rate
        + rateA                     ::Float64            ... Total Auger rate
        + omegaR                    ::EmProperty         ... fluorescence yield
        + omegaA                    ::EmProperty         ... Auger yield
    """
    struct Outcome 
        level                       ::Level 
        NoPhotonLines               ::Int64
        NoAugerLines                ::Int64 
        rateR                       ::EmProperty
        rateA                       ::Float64
        omegaR                      ::EmProperty
        omegaA                      ::EmProperty 
    end 


    """
    `DecayYield.Outcome()`  
        ... constructor for an `empty` instance of DecayYield.Outcome for the computation of fluorescence and Auger yields.
    """
    function Outcome()
        Outcome(Level(), 0, 0, EmProperty(0.), 0., EmProperty(0.), EmProperty(0.))
    end


    # `Base.show(io::IO, outcome::DecayYield.Outcome)`  ... prepares a proper printout of the variable outcome::DecayYield.Outcome.
    function Base.show(io::IO, outcome::DecayYield.Outcome) 
        println(io, "level:                $(outcome.level)  ")
        println(io, "NoPhotonLines:        $(outcome.NoPhotonLines)  ")
        println(io, "NoAugerLines:         $(outcome.NoAugerLines)  ")
        println(io, "NoPhotonLines:        $(outcome.NoPhotonLines)  ")
        println(io, "NoPhotonLines:        $(outcome.NoPhotonLines)  ")
    end


    """
    `DecayYield.computeOutcomes(configs::Array{Configuration,1},asfSettings::AsfSettings, 
                                multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::DecayYield.Settings; output=true)` 
        ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified by 
            the given settings. The results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(initialConfigs::Array{Configuration,1}, asfSettings::AsfSettings, 
                             multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::DecayYield.Settings; output=true)
        println("")
        printstyled("DecayYield.computeOutcomes(): The computation of the fluorescence & Auger yields starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------------- \n", color=:light_green)
        #
        outcomes = DecayYield.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    DecayYield.displayOutcomes(outcomes)    end
        #
        # Overwrite the default settings for the generation of the continuum orbitals, if needed
        ## setDefaults("method: continuum, asymptotic Coulomb")
        ## setDefaults("method: normalization, pure sine")
        #
        println("\nPerform a cascade computation for all radiative decay channels of the levels from the initial configurations:")
        if      settings.approach == "AverageSCA"    cApproach = JAC.Cascade.AverageSCA()
        elseif  settings.approach == "SCA"           cApproach = JAC.Cascade.SCA()
        else    error("Improper (cascade) approach for yield computations; approach = $(settings.approach)")    
        end
        wa = JAC.Cascade.Computation("photon lines", nm, grid, asfSettings, cApproach, [Radiative], 
                                     initialConfigs, Tuple{Int64,Float64}[],  0, 0, Shell[], Shell[], JAC.Cascade.Step[])
        wb = perform(wa; output=false);   linesR = wb.linesR
        println("\nPerform a cascade computation for all (single-electron) Auger decay channels of the levels from the initial configurations:")
        wa = JAC.Cascade.Computation("Auger lines", nm, grid, asfSettings, cApproach, [Auger], 
                                      initialConfigs, Tuple{Int64,Float64}[],  1, 0, Shell[], Shell[], JAC.Cascade.Step[])
        wb = perform(wa; output=false);   linesA = wb.linesA
        #
        # Calculate all amplitudes and requested properties
        newOutcomes = DecayYield.Outcome[]
        for  outcome in outcomes
            newOutcome = JAC.Cascade.computeDecayYieldOutcome(outcome, linesR, linesA, settings)
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        DecayYield.displayResults(stdout, newOutcomes, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    DecayYield.displayResults(iostream, newOutcomes, settings)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `DecayYield.determineOutcomes(multiplet::Multiplet, settings::DecayYield.Settings)`  
        ... to determine a list of Outcomes's for the computation of the alpha-variation parameters for the given 
            multiplet. It takes into account the particular selections and settings. An Array{DecayYield.Outcome,1} is 
            returned. Apart from the level specification, all physical properties are set to zero during the 
            initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::DecayYield.Settings) 
        if    settings.selectLevels   selectLevels   = true;   selectedLevels = copy(settings.selectedLevels)
        else                          selectLevels   = false
        end
    
        outcomes = DecayYield.Outcome[]
        for  i = 1:length(multiplet.levels)
            if  selectLevels  &&  !(haskey(selectedLevels, i))    continue   end
            push!( outcomes, DecayYield.Outcome(multiplet.levels[i], 0, 0, EmProperty(0.), 0., EmProperty(0.), EmProperty(0.)) )
        end
        return( outcomes )
    end


    """
    `DecayYield.displayOutcomes(outcomes::Array{DecayYield.Outcome,1})`  
        ... to display a list of levels that have been selected for the computations. A small neat table of all 
            selected levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{DecayYield.Outcome,1})
        println(" ")
        println("  Selected DecayYield levels:")
        println(" ")
        println("  ", TableStrings.hLine(43))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(43)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
            println( sa )
        end
        println("  ", TableStrings.hLine(43))
        #
        return( nothing )
    end


    """
    `DecayYield.displayResults(stream::IO, outcomes::Array{DecayYield.Outcome,1}, settings::DecayYield.Settings)`  
        ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{DecayYield.Outcome,1}, settings::DecayYield.Settings)
        println(stream, " ")
        println(stream, "  Fluorescence and Auger decay yields:")
        println(stream, " ")
        println(stream, "    + No_R, rate_R, omega_R  ... number of fluorescence lines, total fluorescence rate and yield. ")
        println(stream, "    + No_A, rate_A, omega_A  ... number of Auger lines, total Auger rate and Auger yield. ")
        if    settings.approach in ["AverageSCA", "SCA"]
        println(stream, "    + Approach:  $(settings.approach)  ... all fluorescence rates and yields only in Babushkin gauge. ")
        end
        #
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(149))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=3)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center( 9, "No_R"    ; na=0);                          sb = sb * TableStrings.hBlank(9)
        sa = sa * TableStrings.center(24, "Cou -- rate_R  -- Bab"; na=1);             sb = sb * TableStrings.hBlank(25)
        sa = sa * TableStrings.center(24, "Cou -- omega_R -- Bab"; na=3);             sb = sb * TableStrings.hBlank(27)
        sa = sa * TableStrings.center( 6, "No_A"    ; na=3);                          sb = sb * TableStrings.hBlank(9)
        sa = sa * TableStrings.center(10, "rate_A"; na=2);                            sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(24, "Cou -- omega_A -- Bab"; na=3);             sb = sb * TableStrings.hBlank(27)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(149)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            energy = outcome.level.energy
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))    * "     "
            sa = sa * string(outcome.NoPhotonLines, pad=4)                                      * "    "
            sa = sa * @sprintf("%.4e", outcome.rateR.Coulomb)                                   * " "
            sa = sa * @sprintf("%.4e", outcome.rateR.Babushkin)                                 * "    "
            sa = sa * @sprintf("%.4e", outcome.omegaR.Coulomb)                                  * " "
            sa = sa * @sprintf("%.4e", outcome.omegaR.Babushkin)                                * "     "
            sa = sa * string(outcome.NoAugerLines, pad=4)                                       * "    "
            sa = sa * @sprintf("%.4e", outcome.rateA)                                           * "    "
            sa = sa * @sprintf("%.4e", outcome.omegaA.Coulomb)                                  * " "
            sa = sa * @sprintf("%.4e", outcome.omegaA.Babushkin)                                * "  "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(149), "\n\n")
        #
        return( nothing )
    end


end # module
