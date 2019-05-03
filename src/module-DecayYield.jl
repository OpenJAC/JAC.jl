
"""
`module  JAC.DecayYield`  
    ... a submodel of JAC that contains all methods for computing fluorescence and Auger yields for some level(s).
"""
module DecayYield

    using Printf, JAC, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..TableStrings

    """
    `struct  DecayYield.Settings`  ... defines a type for the details and parameters of computing fluorescence and Auger yields.

        + calcX                    ::Bool             ... True if the ... need to be calculated, and false otherwise.
        + printBeforeComputation   ::Bool             ... True if a list of selected levels is to be printed before the 
                                                          actual computations start. 
        + selectLevels             ::Bool             ... True if individual levels are selected for the computation.
        + selectedLevels           ::Array{Level,1}   ... List of selected levels.
    """
    struct Settings 
        calcX                      ::Bool
        printBeforeComputation     ::Bool
        selectLevels               ::Bool
        selectedLevels             ::Array{Level,1}
    end 


    """
    `DecayYield.Settings()`  
        ... constructor for an `empty` instance of DecayYield.Settings for the computation of fluorescence and Auger yields.
    """
    function Settings()
        Settings(false, false, false, Level[])
    end


    # `Base.show(io::IO, settings::DecayYield.Settings)`  ... prepares a proper printout of the variable settings::DecayYield.Settings.
    function Base.show(io::IO, settings::DecayYield.Settings) 
        println(io, "calcX:                    $(settings.calcX)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
    end


    
    """
    `struct  DecayYield.Outcome`  
        ... defines a type to keep the outcome of a fluorescence and Auger yield computation as well other results.

        + level                     ::Level              ... Atomic level to which the outcome refers to.
        + K                         ::Float64            ... K enhancement parameter
    """
    struct Outcome 
        level                       ::Level 
        K                           ::Float64
    end 


    """
    `DecayYield.Outcome()`  
        ... constructor for an `empty` instance of DecayYield.Outcome for the computation of fluorescence and Auger yields.
    """
    function Outcome()
        Outcome(Level(), 0.)
    end


    # `Base.show(io::IO, outcome::DecayYield.Outcome)`  ... prepares a proper printout of the variable outcome::DecayYield.Outcome.
    function Base.show(io::IO, outcome::DecayYield.Outcome) 
        println(io, "level:                   $(outcome.level)  ")
        println(io, "K:                       $(outcome.K)  ")
    end


    """
    `DecayYield.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::DecayYield.Settings; output=true)` 
        ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified by 
            the given settings. The results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::DecayYield.Settings; output=true)
        println("")
        printstyled("DecayYield.computeOutcomes(): The computation of the fluorescence & Auger yields starts now ... \n", color=:light_green)
        printstyled("----------------------------------------------------------------------------------------------- \n", color=:light_green)
        #
        outcomes = DecayYield.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    DecayYield.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = DecayYield.Outcome[]
        for  outcome in outcomes
            ## newOutcome = DecayYield.computeAmplitudesProperties(outcome, nm, grid, settings) 
            newOutcome = DecayYield.Outcome(outcome.level, 2.0)
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        DecayYield.displayResults(stdout, newOutcomes)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    DecayYield.displayResults(iostream, newOutcomes)   end
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
            push!( outcomes, DecayYield.Outcome(multiplet.levels[i], 0.) )
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
    `DecayYield.displayResults(stream::IO, outcomes::Array{DecayYield.Outcome,1})`  
        ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{DecayYield.Outcome,1})
        println(stream, " ")
        println(stream, "  Fluorescence and Auger decay yields:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(64))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=4)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(14, "K";     na=4)              
        sb = sb * TableStrings.center(14, "    " ; na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(64)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            energy = 1.0
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))    * "    "
            sa = sa * @sprintf("%.8e", outcome.K)                                               * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(64), "\n\n")
        #
        return( nothing )
    end


end # module
