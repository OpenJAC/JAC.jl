
"""
`module  JAC.ReducedDensityMatrix`  
    ... a submodel of JAC that contains all methods for computing reduced density matrices, natural orbitals, radial density
        distributions and related information for some level(s).
"""
module ReducedDensityMatrix

    using Printf, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals,
                  ..SpinAngular, ..TableStrings

    """
    `struct  ReducedDensityMatrix.Settings  <:  AbstractPropertySettings`  
        ... defines a type for providing and parameters for reduced density matrices and other parameters.

        + calc1pRDM                ::Bool             ... True if 1-particle RDM need to be calculated, and false otherwise.
        + calc2pRDM                ::Bool             ... True if 2-particle RDM need to be calculated, and false otherwise.
        + calcCumulant             ::Bool             ... True if the cumulant need to be calculated, and false otherwise.
        + calcNatural              ::Bool             ... True if natural orbitals need to be calculated, and false otherwise.
        + calcDensity              ::Bool             ... True if the radial density need to be calculated, and false otherwise.
        + calcIpq                  ::Bool             ... True if orbital interaction I_pq need to be calculated, and false otherwise.
        + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
        + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
    """
    struct Settings  <:  AbstractPropertySettings 
        calc1pRDM                  ::Bool 
        calc2pRDM                  ::Bool 
        calcCumulant               ::Bool
        calcNatural                ::Bool 
        calcDensity                ::Bool 
        calcIpq                    ::Bool 
        printBefore                ::Bool 
        levelSelection             ::LevelSelection
    end 


    """
    `ReducedDensityMatrix.Settings(; calc1pRDM::Bool=true,` calc2pRDM::Bool=false, calcCumulant::Bool=false, calcNatural::Bool=false, 
                                     calcDensity::Bool=false, calcIpq::Bool=false, printBefore::Bool=true, 
                                     levelSelection::LevelSelection=LevelSelection()) 
        ... keyword constructor to overwrite selected value of reduced-density matrix computations.
    """
    function Settings(; calc1pRDM::Bool=true, calc2pRDM::Bool=false, calcCumulant::Bool=false, calcNatural::Bool=false, 
                        calcDensity::Bool=false, calcIpq::Bool=false, printBefore::Bool=true, 
                        levelSelection::LevelSelection=LevelSelection())
        Settings(calc1pRDM, calc2pRDM, calcCumulant, calcNatural, calcDensity, calcIpq, printBefore, levelSelection)
    end

    # `Base.show(io::IO, settings::ReducedDensityMatrix.Settings)`  
    #       ... prepares a proper printout of the variable settings::ReducedDensityMatrix.Settings.
    function Base.show(io::IO, settings::ReducedDensityMatrix.Settings) 
        println(io, "calc1pRDM:                $(settings.calc1pRDM)  ")
        println(io, "calc2pRDM:                $(settings.calc2pRDM)  ")
        println(io, "calcCumulant:             $(settings.calcCumulant)  ")
        println(io, "calcNatural:              $(settings.calcNatural)  ")
        println(io, "calcDensity:              $(settings.calcDensity)  ")
        println(io, "calcIpq:                  $(settings.calcIpq)  ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "levelSelection:           $(settings.levelSelection)  ")
    end


    """
    `struct  ReducedDensityMatrix.Outcome`  
        ... defines a type to keep the outcome of a reduced density-matrix and natural-orbital computation.

        + level                     ::Level                 ... Atomic level to which the outcome refers to.
        + reduced1p                 ::Array{Float64,2}      ... 1-particle RDM (ip,iq)
        + reduced2p                 ::Array{Float64,4}      ... 1-particle RDM (ip,iq,ir,is)
        + orbitalInteraction        ::Array{Float64,2}      ... orbital interaction I_pq = (ip,iq)
        + naturalOccupation         ::Array{Float64,1}      ... Occupation numbers of natural orbitals.
        + naturalSubshells          ::Array{Subshell,1}     ... List of natural orbitals (subshells).
        + naturalOrbitals           ::Dict{Subshell, Orbital} . Dictionary of (one-electron) natural orbitals.
        + electronDensity           ::Radial.Density        ... Radial density distribution.
    """
    struct Outcome 
        level                       ::Level 
        reduced1p                   ::Array{Float64,2}
        reduced2p                   ::Array{Float64,4}
        orbitalInteraction          ::Array{Float64,2}
        naturalOccupation           ::Array{Float64,1}
        naturalSubshells            ::Array{Subshell,1}
        naturalOrbitals             ::Dict{Subshell, Orbital}
        electronDensity             ::Radial.Density
    end 


    """
    `ReducedDensityMatrix.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of isotope-shift properties.
    """
    function Outcome()
        Outcome(Level(), zeros(1,1), zeros(1,1,1,1), zeros(1,1), Float64[], Subshell[], Dict{Subshell, Orbital}(), Radial.Density() )
    end


    # `Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome)`  ... prepares a proper printout of the variable outcome::ReducedDensityMatrix.Outcome.
    function Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome) 
        println(io, "level:                   $(outcome.level)  ")
        println(io, "reduced1p:               $(outcome.reduced1p)  ")
        println(io, "reduced2p:               $(outcome.reduced2p)  ")
        println(io, "orbitalInteraction:      $(outcome.orbitalInteraction)  ")
        println(io, "naturalOccupation:       $(outcome.naturalOccupation)  ")
        println(io, "naturalSubshells:        $(outcome.naturalSubshells)  ")
        println(io, "naturalOrbitals:         $(outcome.naturalOrbitals)  ")
        println(io, "electronDensity:         $(outcome.electronDensity)  ")
    end

    
    """
    `ReducedDensityMatrix.computeProperties(outcome::ReducedDensityMatrix.Outcome, nm::Nuclear.Model, grid::Radial.Grid, 
                                            settings::ReducedDensityMatrix.Settings) 
        ... to compute all properties for a given level; an outcome::ReducedDensityMatrix.Outcome is returned for 
            which the properties are now evaluated explicitly.
    """
    function  computeProperties(outcome::ReducedDensityMatrix.Outcome, nm::Nuclear.Model, grid::Radial.Grid,
                                settings::ReducedDensityMatrix.Settings)
        reduced1p = outcome.reduced1p;     naturalOccupation = outcome.naturalOccupation;       
        reduced2p = outcome.reduced2p;     naturalSubshells  = outcome.naturalSubshells;     naturalOrbitals = outcome.naturalOrbitals;
        orbitalInteraction = outcome.orbitalInteraction;       electronDensity = outcome.electronDensity;    
        
        if  settings.calc1pRDM      @warn("Option calc1pRDM not yet implemented.")
        end

        if  settings.calc2pRDM      @warn("Option calc2pRDM not yet implemented.")
        end

        if  settings.calcCumulant   println("Option calcCumulant not yet implemented.")
        end

        if  settings.calcNatural    println("Option calcNatural not yet implemented.")
        end

        if  settings.calcDensity    println("Option calcDensity not yet implemented.")
        end

        if  settings.calcIpq        println("Option calcIpq not yet implemented.")
        end

        newOutcome = ReducedDensityMatrix.Outcome( outcome.level, reduced1p, reduced2p, orbitalInteraction, 
                                                   naturalOccupation, naturalSubshells, naturalOrbitals, electronDensity )
        return( newOutcome )
    end


    """
    `ReducedDensityMatrix.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                          settings::ReducedDensityMatrix.Settings; output=true)` 
        ... to compute (as selected) the 1-particle and 2-particle reduced density matrices, natural orbitals or other
            requested information for the levels of interest.
    """
    function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                             settings::ReducedDensityMatrix.Settings; output=true)
        println("")
        printstyled("ReducedDensityMatrix.computeOutcomes(): The computation of the reduced density matrices starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------------------ \n", color=:light_green)
        #
        outcomes = ReducedDensityMatrix.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBefore    ReducedDensityMatrix.displayOutcomes(outcomes, settings)    end
        # Calculate all requested properties
        newOutcomes = ReducedDensityMatrix.Outcome[]
        for  outcome in outcomes
            newOutcome = ReducedDensityMatrix.computeProperties(outcome, nm, grid, settings) 
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        ReducedDensityMatrix.displayResults(stdout, newOutcomes, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    ReducedDensityMatrix.displayResults(iostream, newOutcomes, settings)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `ReducedDensityMatrix.determineOutcomes(multiplet::Multiplet, settings::ReducedDensityMatrix.Settings)`  
        ... to determine a list of Outcomes's for the computation of reduced density matrices and other information 
            for the given multiplet. It takes into account the particular selections and settings. 
            An Array{ReducedDensityMatrix.Outcome,1} is returned. Apart from the level specification, all physical 
            properties are set to zero during the initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::ReducedDensityMatrix.Settings) 
        outcomes = ReducedDensityMatrix.Outcome[]
        @show settings.levelSelection
        for  level  in  multiplet.levels
            if  Basics.selectLevel(level, settings.levelSelection)
                push!( outcomes, ReducedDensityMatrix.Outcome(level, zeros(1,1), zeros(1,1,1,1), zeros(1,1), 
                                                              Float64[], Subshell[], Dict{Subshell, Orbital}(), Radial.Density()) )
            end
        end
        return( outcomes )
    end


    """
    `ReducedDensityMatrix.displayOutcomes(outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)`  
        ... to display a list of levels that have been selected for the computations. A small neat table of all 
            selected levels and their energies is printed but nothing is returned otherwise. Moreover, the
            selected properties and information is printed as well, though not yet calculated.
    """
    function  displayOutcomes(outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)
        nx = 43
        println(" ")
        println("  Selected levels for reduced density matrix & natural orbital calculations:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
            println( sa )
        end
        println("  ", TableStrings.hLine(nx))
        println("\n  The following properties and information will be computed:\n")
        if  settings.calc1pRDM      println("  ++ compute 1-particle RDM.")                       end
        if  settings.calc2pRDM      println("  ++ compute 2-particle RDM.")                       end
        if  settings.calcCumulant   println("  ++ compute cumulant.")                             end
        if  settings.calcNatural    println("  ++ compute natural orbitals & occupation.")        end
        if  settings.calcDensity    println("  ++ compute radial density.")                       end
        if  settings.calcIpq        println("  ++ compute orbital interaction I(p,q).")           end
        #
        return( nothing )
    end


    """
    `ReducedDensityMatrix.displayResults(stream::IO, outcomes::Array{ReducedDensityMatrix.Outcome,1}, 
                                         settings::ReducedDensityMatrix.Settings)`  
        ... to display the reduced density matrices, natural orbitals, etc. for the selected levels. 
            A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)
        nx = 133
        println(stream, " ")
        println(stream, "  ReducedDensityMatrix parameters and information:  !! Not yet worked out !! ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  ";   sc = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5);  sc = sc * TableStrings.hBlank(45)
        sa = sa * TableStrings.center(11, "K_nms"; na=4)              
        sb = sb * TableStrings.center(11, "[GHz u]" ; na=4)
        sc = sc * TableStrings.center(11, "[a.u.]" ; na=4)
        println(stream, sa);    println(stream, sb);    ## println(stream, sc);    
        println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy))            * "    "
            println(stream, sa )
            ## sb = TableStrings.hBlank(47)
        end
        println(stream, "  ", TableStrings.hLine(nx), "\n\n")

        return( nothing )
    end

end # module
