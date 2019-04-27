
"""
`module  JAC.GreenFunction`  ... a submodel of JAC that contains all methods for computing fluorescence and Auger yields for some level; 
                              it is using JAC, JAC.ManyElectron, JAC.Radial.
"""
module GreenFunction

    using Printf, JAC, ..Basics,  ..Basics,  ..Defaults, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0


    """
    `struct  GreenFunction.Settings`  ... defines a type for the details and parameters of computing fluorescence and Auger yields.

        + calcX                    ::Bool             ... True if the ... need to be calculated, and false otherwise.
        + printBeforeComputation   ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
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
    `JAC.GreenFunction.Settings()`  ... constructor for an `empty` instance of GreenFunction.Settings for the computation of fluorescence and Auger yields.
    """
    function Settings()
        Settings(false, false, false, Level[])
    end


    # `Base.show(io::IO, settings::GreenFunction.Settings)`  ... prepares a proper printout of the variable settings::GreenFunction.Settings.
    function Base.show(io::IO, settings::GreenFunction.Settings) 
        println(io, "calcX:                    $(settings.calcX)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
    end


    
    """
    `struct  GreenFunction.Outcome`  ... defines a type to keep the outcome of a fluorescence and Auger yield computation
                                      as well other results.

        + level                     ::Level              ... Atomic level to which the outcome refers to.
        + K                         ::Float64            ... K enhancement parameter
    """
    struct Outcome 
        level                       ::Level 
        K                           ::Float64
    end 


    """
    `JAC.GreenFunction.Outcome()`  ... constructor for an `empty` instance of GreenFunction.Outcome for the computation of fluorescence and Auger yields.
    """
    function Outcome()
        Outcome(Level(), 0.)
    end


    # `Base.show(io::IO, outcome::GreenFunction.Outcome)`  ... prepares a proper printout of the variable outcome::GreenFunction.Outcome.
    function Base.show(io::IO, outcome::GreenFunction.Outcome) 
        println(io, "level:                   $(outcome.level)  ")
        println(io, "K:                       $(outcome.K)  ")
    end


    """
    `JAC.GreenFunction.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::GreenFunction.Settings; output=true)` 
         ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified by the given settings. 
         The results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::GreenFunction.Settings; output=true)
        println("")
        printstyled("JAC.GreenFunction.computeOutcomes(): The computation of approximate Green functions starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------------------- \n", color=:light_green)
        #
        outcomes = JAC.GreenFunction.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    JAC.GreenFunction.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = GreenFunction.Outcome[]
        for  outcome in outcomes
            ## newOutcome = JAC.GreenFunction.computeAmplitudesProperties(outcome, nm, grid, settings) 
            newOutcome = JAC.GreenFunction.Outcome(outcome.level, 2.0)
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        JAC.GreenFunction.displayResults(stdout, newOutcomes)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    JAC.GreenFunction.displayResults(iostream, newOutcomes)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `JAC.GreenFunction.determineOutcomes(multiplet::Multiplet, settings::GreenFunction.Settings)`  ... to determine a list of Outcomes's for 
         the computation of the alpha-variation parameters for the given multiplet. It takes into account the particular selections and 
         settings. An Array{GreenFunction.Outcome,1} is returned. Apart from the level specification, all physical properties are set to zero 
         during the initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::GreenFunction.Settings) 
        if    settings.selectLevels   selectLevels   = true;   selectedLevels = copy(settings.selectedLevels)
        else                          selectLevels   = false
        end
    
        outcomes = GreenFunction.Outcome[]
        for  i = 1:length(multiplet.levels)
            if  selectLevels  &&  !(haskey(selectedLevels, i))    continue   end
            push!( outcomes, GreenFunction.Outcome(multiplet.levels[i], 0.) )
        end
        return( outcomes )
    end


    """
    `JAC.GreenFunction.displayOutcomes(outcomes::Array{GreenFunction.Outcome,1})`  ... to display a list of levels that have been selected 
         for the computations. A small neat table of all selected levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{GreenFunction.Outcome,1})
        println(" ")
        println("  Selected GreenFunction levels:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(43))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(43)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.level.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
            println( sa )
        end
        println("  ", JAC.TableStrings.hLine(43))
        #
        return( nothing )
    end


    """
    `JAC.GreenFunction.displayResults(stream::IO, outcomes::Array{GreenFunction.Outcome,1})`  ... to display the energies, M_ms and F-parameters, etc. for the 
         selected levels. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{GreenFunction.Outcome,1})
        println(stream, " ")
        println(stream, "  Green's function summary:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(64))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4)              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(14, "K";     na=4)              
        sb = sb * JAC.TableStrings.center(14, "    " ; na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(64)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.level.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            energy = 1.0
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))              * "    "
            sa = sa * @sprintf("%.8e", outcome.K)                                               * "    "
            println(stream, sa )
        end
        println(stream, "  ", JAC.TableStrings.hLine(64), "\n\n")
        #
        return( nothing )
    end


end # module
