
"""
`module  JAC.GreenFunction`  
    ... a submodel of JAC that contains all methods for computing approximate many-electron Green functions.
"""
module GreenFunction

    using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings


    """
    `GreenFunction.computeRepresentation(configs::Array{Configuration,1}, multiplet::Multiplet, nm::Nuclear.Model, 
                                         grid::Radial.Grid, settings::GreenFunction.Settings; output=true)` 
        ... to compute an approximate Green function representation for the levels of the given bound configurations and for the
            given excitation scheme.
    """
    function computeRepresentation(configs::Array{Configuration,1}, multiplet::Multiplet, nm::Nuclear.Model, 
                                   grid::Radial.Grid, settings::GreenFunction.Settings; output=true)
        #
        if    output    return( gfRep )
        else            return( nothing )
        end
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
