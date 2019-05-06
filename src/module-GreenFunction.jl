
"""
`module  JAC.GreenFunction`  
    ... a submodel of JAC that contains all methods for computing approximate many-electron Green functions.
"""
module GreenFunction

    using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings

    """
    `abstract type GreenFunction.AbstractGreenApproach` 
        ... defines an abstract and a number of singleton types for approximating a many-electron Green
            function for calculating second-order processes.

      + struct SingleCSFwithoutCI        
        ... to approximate the many-electron multiplets (gMultiplet) for every chosen level symmetry by dealing with each CSF 
            independently, and without any configuration interaction. This is a fast but very rough approximation.
            
      + struct CoreSpaceCI                
        ... to approximate the many-electron multiplets (gMultiplet) by including configuration between the bound-state
            orbitals into account.
    """
    abstract type  AbstractGreenApproach                           end
    struct         SingleCSFwithoutCI  <:  AbstractGreenApproach   end
    struct         CoreSpaceCI         <:  AbstractGreenApproach   end

    """
    `struct  GreenFunction.Settings`  
        ... defines a type for defining the details and parameters of the approximated Green functions.

        + approach                 ::AbstractGreenApproach    ... Approach used to approximate the representation.
        + excitationScheme         ::Basics.AbstractConfigurationScheme    
                                                      ... Applied excitation scheme with regard to the given bound 
                                                          configuration.
        + nMax                     ::Int64            ... maximum principal quantum numbers of (single-electron) 
                                                          excitations to be included into the representation.
        + lValues                  ::Array{Int64,1}   ... List of (non-relativistic) orbital angular momenta for which
                                                          (single-electron) excitations are to be included.
        + levelSymmetries          ::Array{LevelSymmetry,1}   ... Symmetries J^P to be included into the Green 
                                                                  function representation.
        + printBeforeComputation   ::Bool             ... True if a short overview is to be printed before. 
        + selectLevels             ::Bool             ... True if individual levels are selected for the computation.
        + selectedLevels           ::Array{Int64,1}   ... List of selected levels.
    """
    struct Settings 
        approach                   ::AbstractGreenApproach
        excitationScheme           ::Basics.AbstractConfigurationScheme
        nMax                       ::Int64
        lValues                    ::Array{Int64,1}
        levelSymmetries            ::Array{LevelSymmetry,1}
        printBeforeComputation     ::Bool 
        selectLevels               ::Bool
        selectedLevels             ::Array{Int64,1}
    end 


    """
    `GreenFunction.Settings()`  ... constructor for an `empty` instance of GreenFunction.
    """
    function Settings()
        Settings( GreenFunction.SingleCSFwithoutCI(), Basics.NoExcitationScheme(), 0, Int64[], LevelSymmetry[], false, false, Int64[])
    end


    # `Base.show(io::IO, settings::GreenFunction.Settings)`  ... prepares a proper printout of the variable settings::GreenFunction.Settings.
    function Base.show(io::IO, settings::GreenFunction.Settings) 
        println(io, "approach:                 $(settings.approach)  ")
        println(io, "excitationScheme:         $(settings.excitationScheme)  ")
        println(io, "nMax:                     $(settings.nMax)  ")
        println(io, "lValues:                  $(settings.lValues)  ")
        println(io, "levelSymmetries:          $(settings.levelSymmetries)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
    end


    """
    `struct  GreenFunction.Channel`  
        ... defines a type for a single symmetry channel of an (approximate) Green function.

        + symmetry          ::LevelSymmetry    ... Level symmetry of this part of the representation.
        + gMultiplet        ::Multiplet        ... Multiplet of (scattering) levels of this symmetry.
    """
    struct Channel 
        symmetry            ::LevelSymmetry
        gMultiplet          ::Multiplet
    end   


    """
    `GreenFunction.Channel()`  ... constructor for an `empty` instance of GreenFunction.Channel.
    """
    function Channel()
        Channel( LevelSymmetry(0, Basics.plus), ManyElectron.Multiplet)
    end


    # `Base.show(io::IO, channel::GreenFunction.Channel)`  ... prepares a proper printout of the variable channel::GreenFunction.Channel.
    function Base.show(io::IO, channel::GreenFunction.Channel) 
        println(io, "symmetry:                $(channel.symmetry)  ")
        println(io, "gMultiplet:              $(channel.gMultiplet)  ")
    end


    """
    `struct  GreenFunction.Representation`  
        ... defines a type to keep a representation of an (approximate) Green function that is associated
            with a given set of bound configurations.

        + approach          ::AbstractGreenApproach          ... Approach used to approximate the representation.
        + excitationScheme  ::Basics.AbstractConfigurationScheme    ... Applied excitation scheme with regard to the
                                                                 given bound configuration.
        + boundConfigs      ::Array{Configuration,1}         ... List of bound configurations to which the 
                                                                 Green function refers to.
        + NoElectrons       ::Int64                          ... Number of electrons
        + channels          ::Array{GreenFunction.Channel,1} ... List of channels with different level symmetry.
    """
    struct Representation 
        approach            ::AbstractGreenApproach
        excitationScheme    ::Basics.AbstractConfigurationScheme 
        boundConfigs        ::Array{Configuration,1}
        NoElectrons         ::Int64 
        channels            ::Array{GreenFunction.Channel,1}
    end   


    """
    `GreenFunction.Representation()`  ... constructor for an `empty` instance of GreenFunction.Representation.
    """
    function Representation()
        Representation( GreenFunction.SingleCSFwithoutCI(), Basics.NoExcitationScheme(), Configuration[], 0, GreenFunction.Channel[])
    end


    # `Base.show(io::IO, representation::GreenFunction.Representation)`  ... prepares a proper printout of the 
    #       variable representation::GreenFunction.Representation.
    function Base.show(io::IO, representation::GreenFunction.Representation) 
        println(io, "approach:                $(representation.approach)  ")
        println(io, "excitationScheme:        $(representation.excitationScheme)  ")
        println(io, "boundConfigs:            $(representation.boundConfigs)  ")
        println(io, "NoElectrons:             $(representation.NoElectrons)  ")
        println(io, "channels:                $(representation.channels)  ")
    end


    """
    `GreenFunction.computeRepresentation(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::GreenFunction.Settings; output=true)` 
        ... to compute (as selected) an approximate Green function representation for the levels of the given bound configurations.
    """
    function computeRepresentation(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::GreenFunction.Settings; output=true)
        println("")
        printstyled("GreenFunction.computeRepresentation(): The computation of approximate Green functions starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------- \n", color=:light_green)
        #
        # Generate all configurations from the bound configurations due to the given scheme ... and print them, if required
        # Print further information about the requested Greens function representation ... to stdout and summary file
        # Determine all (selected) levels and generate an averaged potential for these levels
        # Generate a full single-electron spectrum for this potential
        # Cycle over all selected level symmetries
        # * Generate all CSFs for this symmetry
        # * Define a basis for this symmetry
        # * Calculate a multiplet representation due to the given approach
        # Set up the final representation
        # Summarize and display the results about this representation in a neat format
        # Return the representation if required
        error("stop Greens here")
        outcomes = GreenFunction.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    GreenFunction.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = GreenFunction.Outcome[]
        for  outcome in outcomes
            ## newOutcome = GreenFunction.computeAmplitudesProperties(outcome, nm, grid, settings) 
            newOutcome = GreenFunction.Outcome(outcome.level, 2.0)
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        GreenFunction.displayResults(stdout, newOutcomes)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    GreenFunction.displayResults(iostream, newOutcomes)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end

    
    
    
    
    
    
    

    """
    `GreenFunction.determineOutcomes(multiplet::Multiplet, settings::GreenFunction.Settings)`  
        ... to determine a list of Outcomes's for the computation of the alpha-variation parameters for the 
            given multiplet. It takes into account the particular selections and settings. An Array{GreenFunction.Outcome,1} 
            is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
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
    `GreenFunction.displayRepresentation(representation::GreenFunction.Representation)`  
        ... to display
    """
    function  displayRepresentation(representation::GreenFunction.Representation)
        println(" ")
        println("  Selected GreenFunction levels:")
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
    `GreenFunction.displayResults(stream::IO, representation::GreenFunction.Representation)`  
        ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but nothing is 
            returned otherwise.
    """
    function  displayResults(stream::IO, representation::GreenFunction.Representation)
        println(stream, " ")
        println(stream, "  Green's function summary:")
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
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))              * "    "
            sa = sa * @sprintf("%.8e", outcome.K)                                               * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(64), "\n\n")
        #
        return( nothing )
    end

end # module
