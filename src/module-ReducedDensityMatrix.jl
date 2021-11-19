
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

        + calcNatural              ::Bool             ... True if natural orbitals need to be calculated, and false otherwise.
        + calcDensity              ::Bool             ... True if the radial density need to be calculated, and false otherwise.
        + calcIpq                  ::Bool             ... True if orbital interaction I_pq need to be calculated, and false otherwise.
        + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
        + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
    """
    struct Settings  <:  AbstractPropertySettings 
        calcNatural                ::Bool 
        calcDensity                ::Bool 
        calcIpq                    ::Bool 
        printBefore                ::Bool 
        levelSelection             ::LevelSelection
    end 


    """
    `ReducedDensityMatrix.Settings(; calcNatural::Bool=true, calcDensity::Bool=true, calcIpq::Bool=false, 
                                     printBefore::Bool=true, levelSelection::LevelSelection=LevelSelection()) 
        ... keyword constructor to overwrite selected value of reduced-density matrix computations.
    """
    function Settings(; calcNatural::Bool=true, calcDensity::Bool=true, calcIpq::Bool=false, 
                                     printBefore::Bool=true, levelSelection::LevelSelection=LevelSelection())
        Settings(calcNatural, calcDensity, calcIpq, printBefore, levelSelection)
    end

    # `Base.show(io::IO, settings::ReducedDensityMatrix.Settings)`  
    #       ... prepares a proper printout of the variable settings::ReducedDensityMatrix.Settings.
    function Base.show(io::IO, settings::ReducedDensityMatrix.Settings) 
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
        + naturalSubshells          ::Array{Subshell,1}     ... List of natural orbitals (subshells).
        + naturalOccupation         ::Array{Float64,1}      ... Occupation numbers of natural orbitals.
        + naturalOrbitalExpansion   ::Dict{Subshell, Array{Float64,1}} 
            ... Dictionary of the expansion of (one-electron) natural orbitals in terms of the standard orbitals
                from the given basis.
        + naturalOrbitals           ::Dict{Subshell, Orbital} . Dictionary of (one-electron) natural orbitals.
        + orbitalInteraction        ::Array{Float64,2}      ... orbital interaction I_pq = (ip,iq)
        + electronDensity           ::Radial.Density        ... Radial density distribution.
    """
    struct Outcome 
        level                       ::Level 
        naturalSubshells            ::Array{Subshell,1}
        naturalOccupation           ::Array{Float64,1}
        naturalOrbitalExpansion     ::Dict{Subshell, Array{Float64,1}}
        naturalOrbitals             ::Dict{Subshell, Orbital}
        orbitalInteraction          ::Array{Float64,2}
        electronDensity             ::Radial.Density
    end 


    """
    `ReducedDensityMatrix.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of isotope-shift properties.
    """
    function Outcome()
        Outcome(Level(), Subshell[], Float64[], Dict{Subshell, Array{Float64,1}}, Dict{Subshell, Orbital}(), zeros(1,1), Radial.Density() )
    end


    # `Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome)`  ... prepares a proper printout of the variable outcome::ReducedDensityMatrix.Outcome.
    function Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome) 
        println(io, "level:                   $(outcome.level)  ")
        println(io, "naturalSubshells:        $(outcome.naturalSubshells)  ")
        println(io, "naturalOccupation:       $(outcome.naturalOccupation)  ")
        println(io, "naturalOrbitalExpansion: $(outcome.naturalOrbitalExpansion)  ")
        println(io, "naturalOrbitals:         $(outcome.naturalOrbitals)  ")
        println(io, "orbitalInteraction:      $(outcome.orbitalInteraction)  ")
        println(io, "electronDensity:         $(outcome.electronDensity)  ")
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
        # Calculate all requested properties for each outcome
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
    `ReducedDensityMatrix.computeNaturalOrbitalExpansion(naturalSubshells::Array{Subshell,1}, level::Level)`  
        ... to perform the expansion of the natural orbitals for the given level in terms of the standard orbitals
            as defined by the basis. A tuple (naturalOcc::Array{Float64,1}, naturalExp::Dict{Subshell, Array{Float64,1}})
            is returned which provides the natural occpuation numbers and expansion coefficients with regard
            to the given list of natural orbitals. This method makes use of angular coefficients for a zero-rank operator 
            to calculate and diagonalize the expansion matrix.
    """
    function  computeNaturalOrbitalExpansion(naturalSubshells::Array{Subshell,1}, level::Level)
        naturalOcc = Float64[], naturalExp = Dict{Subshell, Array{Float64,1}}()

        return( naturalOcc, naturalExp )
    end


    """
    `ReducedDensityMatrix.computeNaturalOrbitals(naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)`  
        ... to compute the natural orbitals as superposition of the standard orbitals; for each subshell in naturalExp,
            an orbital is computed, and naturalOrbitals::Dict{Subshell, Orbital} returned
    """
    function  computeNaturalOrbitals(naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)
        naturalOrbitals = Dict{Subshell, Orbital}()

        return( naturalOrbitals )
    end

    
    """
    `ReducedDensityMatrix.computeProperties(outcome::ReducedDensityMatrix.Outcome, nm::Nuclear.Model, grid::Radial.Grid, 
                                            settings::ReducedDensityMatrix.Settings) 
        ... to compute all properties for a given level; an outcome::ReducedDensityMatrix.Outcome is returned for 
            which the properties are now evaluated explicitly.
    """
    function  computeProperties(outcome::ReducedDensityMatrix.Outcome, nm::Nuclear.Model, grid::Radial.Grid,
                                settings::ReducedDensityMatrix.Settings)
        naturalOccupation  = outcome.naturalOccupation;        naturalOccupationExp = outcome.naturalOccupationExpansion;  
        naturalSubshells   = outcome.naturalSubshells;         naturalOrbitals      = outcome.naturalOrbitals;
        orbitalInteraction = outcome.orbitalInteraction;       electronDensity      = outcome.electronDensity;    
        
        if  settings.calcNatural    println("Option calcNatural not yet implemented.")
        end

        if  settings.calcDensity    println("Option calcDensity not yet implemented.")
        end

        if  settings.calcIpq        println("Option calcIpq not yet implemented.")
        end

        newOutcome = ReducedDensityMatrix.Outcome( outcome.level, naturalOccupationExp, naturalOccupation, 
                                                   naturalSubshells, naturalOrbitals, orbitalInteraction, electronDensity )
        return( newOutcome )
    end


    """
    `ReducedDensityMatrix.computeRadialDistribution(naturalExp::Dict{Subshell, Array{Float64,1}}, grid::Radial.Grid, level::Level)`  
        ... to compute the radial distribution function from the given standard orbitals and natural expansion
            coefficients. An radial disribution D::Array{Float64,1} is returned that specifies the radial distribution
            w.r.t. the given grid.
    """
    function  computeRadialDistribution(naturalExp::Dict{Subshell, Array{Float64,1}}, grid::Radial.Grid, level::Level)
        D = Float64[]

        return( D )
    end


    """
    `ReducedDensityMatrix.compute1pRDM(p::Subshell, q::Subshell, naturalSubshells::Array{Subshell,1}, 
                                       naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)`  
        ... to compute 1p RDM for the given pair (p,q) of subshells; a rdm::Array{Float64,2} is returned.
    """
    function  compute1pRDM(p::Subshell, q::Subshell, naturalSubshells::Array{Subshell,1}, 
                           naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)
        ns = length(naturalSubshells);  rdm = zeros(ns,ns)

        return( rdm )
    end


    """
    `ReducedDensityMatrix.compute1pRDM(p::Subshell, q::Subshell, r::Subshell, s::Subshell, naturalSubshells::Array{Subshell,1}, 
                                       naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)`  
        ... to compute 2p RDM for the given pair (p,q;r,s) of subshells; a rdm::Array{Float64,4} is returned.
    """
    function  compute2pRDM(p::Subshell, q::Subshell, r::Subshell, s::Subshell, naturalSubshells::Array{Subshell,1}, 
                           naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)
        ns = length(naturalSubshells);  rdm = zeros(ns,ns,ns,ns)

        return( rdm )
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
                push!( outcomes, ReducedDensityMatrix.Outcome(level, Subshell[], Float64[], Dict{Subshell, Array{Float64,1}}(), 
                                                              Dict{Subshell, Orbital}(), zeros(1,1), Radial.Density()) )
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
        println("  Results for the natural occupation number, natural orbitals and kp RDM are printed in turn for the following levels:")
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
        #
        println("\n\n  The following properties are calculated: \n")
        println("  + Natural occupations numbers")
        println("  + Natural orbital expansion in terms of the standard orbitals.")
        if  settings.calcNatural   println("  + Representation of the natural orbitals.")  end
        if  settings.calcDensity   println("  + Radial electron density distribution.")    end
        if  settings.calcIpq       println("  + Orbital interaction I_pq.")                end
        
        return( nothing )
    end


    """
    `ReducedDensityMatrix.displayResults(outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)`  
        ... to display a list of levels that have been selected for the computations. A small neat table of all 
            selected levels and their energies is printed but nothing is returned otherwise. Moreover, the
            selected properties and information is printed as well, though not yet calculated.
    """
    function  displayResults(outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)
        nx = 43
        println(" ")
        println("  Results for the natural occupation number, natural orbitals and kp RDM are printed in turn for the following levels:")
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
        #
        # Now print all selected result for each outcome
        for  outcome in outcomes
            nx = 43
            println(" ")
            println("  Natural occupation numbers:")
            println(" ")
            println("  ", TableStrings.hLine(nx))
            # Use table.strings ... to split a vector in strings of 10 subshells/numbers ... and which are returned
            # as a tuple.
            # Split  and generate the header of the table
            stuple = ( "aa", "bb")
            for  (i, st)  in  enumerate(stuple)
                if       i == 1  sa = "1****" * st
                else             sa = "     " * st      end
                println(sa)
            end
            println("  ", TableStrings.hLine(nx))
            # Split and generate the occupation numbers of the table
            stuple = ( "cc", "dd")
            for  (i, st)  in  enumerate(stuple)
                if       i == 1  sa = "1****" * st
                else             sa = "     " * st      end
                println(sa)
            end
            println("  ", TableStrings.hLine(nx))
            #
            nx = 43
            println(" ")
            println("  Natural orbital expansion:")
            println(" ")
            println("  ", TableStrings.hLine(nx))
            # Use table.strings ... to split a vector in strings of 10 subshells/numbers ... and which are returned
            # as a tuple.
            # Split  and generate the header of the table
            stuple = ( "aa", "bb")
            for  (i, st)  in  enumerate(stuple)
                if       i == 1  sa = "1****" * st
                else             sa = "     " * st      end
                println(sa)
            end
            println("  ", TableStrings.hLine(nx))
            #
            for subsh in naturalSubshells
                # Split and generate the occupation numbers of the table
                stuple = ( "cc", "dd")
                for  (i, st)  in  enumerate(stuple)
                    if       i == 1  sa = "1****" * st
                    else             sa = "     " * st      end
                end
            end
            println("  ", TableStrings.hLine(nx))
            #
        end  # outcomes
        #
        #
        if  settings.calcNatural
            println("\n  Natural orbitals are calculated for the following subshells and are available by outcome.naturalOrbitals:\n")
            sa = "  ";      for  subsh in outcome.naturalSubshells    sa = sa * "   $subsh"     end 
            println(sa)
        end
        #
        #
        if  settings.calcDensity
            println("\n  Radial electron distribution:\n")
            for  i = 1:20:grid.NoPoints
                println("   " * ("      "*string(i))[end-6:end] * ")    " *
                        @sprintf("%.4e", grid.r[i]) * @sprintf("%.4e", outcome.electronDensity[i])  )
            end
        end
        #
        #
        if  settings.calcIpq
            nx = 43
            println("\n  Orbital interactions I_pq:\n")
            println("  ", TableStrings.hLine(nx))
            # Use table.strings ... to split a vector in strings of 10 subshells/numbers ... and which are returned
            # as a tuple.
            # Split  and generate the header of the table
            stuple = ( "aa", "bb")
            for  (i, st)  in  enumerate(stuple)
                if       i == 1  sa = "NO / NO    " * st
                else             sa = "           " * st      end
                println(sa)
            end
            println("  ", TableStrings.hLine(nx))
            # Split and generate the occupation numbers of the table
            for  subsh in outcome.naturalSubshells
                stuple = ( "cc", "dd")
                for  (i, st)  in  enumerate(stuple)
                    if       i == 1  sa = "1****" * st
                    else             sa = "     " * st      end
                end
            end
            println("  ", TableStrings.hLine(nx))
        end
        #
        return( nothing )
    end

end # module
