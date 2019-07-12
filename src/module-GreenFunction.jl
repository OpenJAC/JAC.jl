
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
            independently, and without any configuration interaction. This is a fast but also very rough approximation.
            
      + struct CoreSpaceCI                
        ... to approximate the many-electron multiplets (gMultiplet) by taking the electron-electron interaction between the 
            bound-state orbitals into account.
            
      + struct DampedSpaceCI                
        ... to approximate the many-electron multiplets (gMultiplet) by taking the electron-electron interaction for all, the 
            bound and free-electron, orbitals into account but by including a damping factor e^{tau*r}, tau > 0 into the 
            electron densities rho_ab (r) --> rho_ab (r) * e^{tau*r}
    """
    abstract type  AbstractGreenApproach                           end
    struct         SingleCSFwithoutCI  <:  AbstractGreenApproach   end
    struct         CoreSpaceCI         <:  AbstractGreenApproach   end
    struct         DampedSpaceCI       <:  AbstractGreenApproach   end

    """
    `struct  GreenFunction.Settings`  
        ... defines a type for defining the details and parameters of the approximated Green functions.

        + approach                 ::AbstractGreenApproach            ... Approach used to approximate the representation.
        + excitationScheme         ::Basics.AbstractExcitationScheme  ... Applied excitation scheme with regard to the 
                                                                          given bound configuration.
        + nMax                     ::Int64            ... maximum principal quantum numbers of (single-electron) 
                                                          excitations that are to be included into the representation.
        + lValues                  ::Array{Int64,1}   ... List of (non-relativistic) orbital angular momenta for which
                                                          (single-electron) excitations are to be included.
        + levelSymmetries          ::Array{LevelSymmetry,1}           ... Total level symmetries J^P to be included into 
                                                                          the  Green function representation.
        + printBeforeComputation   ::Bool             ... True if a short overview is to be printed before. 
        + selectLevels             ::Bool             ... True if individual levels are selected for the computation.
        + selectedLevels           ::Array{Int64,1}   ... List of selected levels.
    """
    struct Settings 
        approach                   ::AbstractGreenApproach
        excitationScheme           ::Basics.AbstractExcitationScheme
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

        + approach          ::AbstractGreenApproach           ... Approach used to approximate the representation.
        + excitationScheme  ::Basics.AbstractExcitationScheme ... Applied excitation scheme with regard to the
                                                                  given bound configuration.
        + boundConfigs      ::Array{Configuration,1}          ... List of bound configurations to which the Green 
                                                                  function refers to.
        + NoElectrons       ::Int64                           ... Number of electrons.
        + channels          ::Array{GreenFunction.Channel,1}  ... List of channels with different level symmetry.
    """
    struct Representation 
        approach            ::AbstractGreenApproach
        excitationScheme    ::Basics.AbstractExcitationScheme 
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
    `GreenFunction.computeRepresentation(configs::Array{Configuration,1}, multiplet::Multiplet, nm::Nuclear.Model, 
                                         grid::Radial.Grid, settings::GreenFunction.Settings; output=true)` 
        ... to compute an approximate Green function representation for the levels of the given bound configurations and for the
            given excitation scheme.
    """
    function computeRepresentation(configs::Array{Configuration,1}, multiplet::Multiplet, nm::Nuclear.Model, 
                                   grid::Radial.Grid, settings::GreenFunction.Settings; output=true)
        println("")
        printstyled("GreenFunction.computeRepresentation(): The computation of approximate Green functions starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------- \n", color=:light_green)
        
        gChannels = GreenFunction.Channel[]
        
        # Generate all (non-relativistic) configurations from the bound configurations due to the given excitation scheme 
        confList = Basics.generateConfigurationsForExcitationScheme(configs, settings.excitationScheme, settings.nMax, settings.lValues)
        # Print (if required) information about the generated configuration list
        if  settings.printBeforeComputation    GreenFunction.displayConfigurationList(confList)    end
        # Determine all (selected) levels and generate an averaged potential for these levels; use a ... potential for the selected levels
        pot = GreenFunction.generateMeanPotential(multiplet, nm, grid, settings)
        # Generate a full single-electron spectrum for this potential
        spectrum = Bspline.generateSpectrum()
        # Cycle over all selected level symmetries to generate in turn all channels of the requires Green function representation
        for  levelSymmetry  in  settings.levelSymmetries
            # Generate all CSFs for this symmetry
            csfs  = Basics.generateCsfList(levelSymmetry, confList)
            # Compute a multiplet representation due to the given approach
            gMultiplet = GreenFunction.computeMultipletForApproach(settings.approach, csfs, spectrum)
            channel    = GreenFunction.Channel( levelSymmetry, gMultiplet )
        end        
        
        gfRep = GreenFunction.Representation( settings.approach, settings.excitationScheme, configs, confList[1].NoElectrons , gChannels) 
        #
        # Print all results to screen
        GreenFunction.displayRepresentation(stdout, gfRep)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    GreenFunction.displayRepresentation(iostream, gfRep)   end
        #
        if    output    return( gfRep )
        else            return( nothing )
        end
    end

    
    
    """
    `GreenFunction.displayConfigurationList(configs::Array{Configuration,1})`  
        ... displays the gnerated list of configurations in a compact form; nothing is returned in this case.
    """
    function  displayConfigurationList(configs::Array{Configuration,1})
        println(" ")
        println("  Generated list of (non-relativistic) configurations that contribute to the Green functions representation:")
        println(" ")
        println("  ", TableStrings.hLine(60))
        for  i = 1:length(configs)
            sa = string("(", i, ")" );    sa = TableStrings.flushright(9, sa, na=3);   sa = sa * string(configs[i])
            println(sa)
        end 
        println("  ", TableStrings.hLine(60))
        
        return( nothing )
    end


    """
    `GreenFunction.displayRepresentation(representation::GreenFunction.Representation)`  
        ... to display
    """
    function  displayRepresentation(representation::GreenFunction.Representation)
        println(" ")
        println("  Selected GreenFunction levels:  ... not yet implemented")
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
    `GreenFunction.generateMeanPotential(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::GreenFunction.Settings)`  
        ... generates a mean radial potential for the selected levels of the multiplet that is used to compute the 
            single-electron spectrum for the representation of the Green function; a potential::Radial.Potential
            is returned.
    """
    function  generateMeanPotential(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::GreenFunction.Settings)
        nuclearPotential  = Nuclear.nuclearPotential(nm, grid)
        ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
        ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
        ## wp3 = compute("radial potential: Kohn-Sham", grid, wLevel)
        ## wp4 = compute("radial potential: Dirac-Fock-Slater", grid, wLevel)
        
        if  settings.selectLevels  
            wpx = zeros[grid.NoPoints];   nx = 0
            for  i = 1:length(multiplet.levels)
                if haskey(settings.selectedLevels, i)   
                    wp = compute("radial potential: Dirac-Fock-Slater", grid, multiplet.level[i]);   np = length(wp.Zr)
                    wpx[1:np] = wpx[1:np] + wp.Zr[1:np];    nx = nx + 1
                end
            end
            wpx = wpx / nx
            wp = Radial.Potential(wp.name, wpx, grid)
        else   
            wp = compute("radial potential: Dirac-Fock-Slater", grid, multiplet.levels[1].basis)   
        end
        
        
        pot = Basics.add(nuclearPotential, wp)
        println("Mean potential for the generation of the (quasi-complete) single-electron spectrum:  $(pot)")
        
        return( pot )
    end

end # module
