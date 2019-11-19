
"""
`module  JAC.BasicsHP`  
    ... a submodel of JAC that contains methods for setting-up and performing atomic structure computations,
        i.e. self-consistent-fields, configuration-interaction, etc.
"""
module BasicsHP

    using Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..Hfs, ..ManyElectron, ..Nuclear, ..PhotoEmission,
                  ..Radial, ..TableStrings
    
    export isless

    """
    `Basics.interpolateOnGridGrasp92((F::Array{Float64,1}, oldGrid::Radial.Grid), (G::Void, newGrid::Radial.Grid) )`  
        ... to interpolate the (radial) function F from the oldgrid to newgrid, by a call to the Grasp92 Fortran procedures; 
            an function G::Array{Float64,1} is returned. **Not yet implemented !**    
    """
    function Basics.interpolateOnGridGrasp92(from::Tuple{Array{Float64,1},Radial.Grid}, to::Tuple{DataType,Radial.Grid})
        # First prepare all fields and settings to call quad_grasp92()
        error("Not yet implemented !")
    end


    """
    `Basics.interpolateOnGridTrapezRule((F::Array{Float64,1}, oldGrid::Radial.Grid), (G::Void, newGrid::Radial.Grid)` 
        ... to integrate the (radial) function F from the oldgrid to newgrid, by using a trapez rule from the Numerical Recipies. 
            A function G::Array{Float64,1} is returned. **Not yet implemented !**.
    """
    function Basics.interpolateOnGridTrapezRule(from::Tuple{Array{Float64,1},Radial.Grid}, to::Tuple{DataType,Radial.Grid})
        error("Not yet implemented !")
    end


    """
    `Basics.integrate()`  
        ... integrates a one- or two-dimensional function by different numerical methods, either on a given grid or for a general function.

    + `("function: on radial grid, Newton-Cotes", F::Array{Float64,1}, grid::Radial.Grid)`  
        ... to integrate the (radial) function F over the given (radial) grid, int_0^infinity dr F(r) by using a 5-point Newton-Cotes 
            integration formula; a value::Float64 is returned.   

    + `("function: on radial grid, Simpson rule", F::Array{Float64,1}, grid::Radial.Grid)`  
        ... to integrate the (radial) function F over the given (radial) grid, int_0^infinity dr F(r) by using Simpson's rule; 
            a value::Float64 is returned. 

    + `("function: on radial grid, trapez rule", F::Array{Float64,1}, grid::Radial.Grid)`  
        ... to integrate the (radial) function F over the given (radial) grid, int_0^infinity dr F(r), by using a simple trapez rule; 
            a value::Float64 is returned. 
    """
    function Basics.integrate(sa::String, F::Array{Float64,1}, grid::Radial.Grid)
        if       sa == "function: on radial grid, Newton-Cotes"     wa = integrateOnGridNewtonCotes(F, grid)   
        elseif   sa == "function: on radial grid, Simpson rule"     wa = integrateOnGridSimpsonRule(F, grid)
        elseif   sa == "function: on radial grid, trapez rule"      wa = integrateOnGridTrapezRule( F, grid)   
        else     error("Unsupported keystring = $sa.")
        end
        
    end


    """
    `Basics.integrateOnGridNewtonCotes()`  ... integrates by using a 5-point Newton-Cotes integration formula; see Basic.integrate().
    """
    function Basics.integrateOnGridNewtonCotes(F::Array{Float64,1}, grid::Radial.Grid)
        coefficients = Array{Float64}([2 * 7, 32, 12, 32, 7]) * 2 / 45 * grid.h
        n = size(coefficients, 1)
    
        result = 0.;    i0     = 1
    
        while  abs(F[i0]) == 0
            i0 += 1
        end
        gamma  = log(F[i0]/F[i0+1] * grid.rp[i0+1]/grid.rp[i0]) / log(grid.r[i0]/grid.r[i0+1])
        result = 1/(gamma + 1) * grid.r[i0] * F[i0] / grid.rp[i0]
    
        result += F[i0] * coefficients[n]
        for i = i0+1:size(F, 1)
            result += coefficients[(i - i0) % (n - 1) + 1] * F[i]
        end
    
        if (size(F, 1) - i0) % (n - 1) + 1 == 1
            result -= coefficients[n] * F[size(F, 1)]
        end
    
        return( result )
    end


    """
    `Basics.integrateOnGridSimpsonRule()`  ... integrates by using Simpson's rule; see Basics.integrate().
    """
    function Basics.integrateOnGridSimpsonRule(F::Array{Float64,1}, grid::Radial.Grid)
        coefficients = Array{Float64}([2 * 1, 4, 1]) / 3 * grid.h
        n = size(coefficients, 1)
    
        result = 0.
        result += F[1] * coefficients[n]
        for i = 2:size(F, 1)
            result += coefficients[(i - 1) % (n - 1) + 1] * F[i]
        end
    
        if (size(F, 1) - 1) % (n - 1) + 1 == 1
            result -= coefficients[n] * F[size(F, 1)]
        end
    
        return( result )
    end


    """
    `Basics.integrateOnGridTrapezRule()`  ... integrates by using a simple trapez rule; see Basics.integrate().
    """
    function Basics.integrateOnGridTrapezRule(F::Array{Float64,1}, grid::Radial.Grid)
        coefficients = Array{Float64}([2 * 1, 1]) / 2 * grid.h
        n = size(coefficients, 1)
    
        result = 0.
    
        result += F[1] * coefficients[n]
        for i = 2:size(F, 1)
            result += coefficients[(i - 1) % (n - 1) + 1] * F[i]
        end
    
        if (size(F, 1) - 1) % (n - 1) + 1 == 1
            result -= coefficients[n] * F[size(F, 1)]
        end
    
        return( result )
    end


    function Base.isless(x::ManyElectron.Level, y::ManyElectron.Level)
        return x.energy < y.energy
    end


    function Base.isless(x::Hfs.IJF_Level, y::Hfs.IJF_Level)
        return x.energy < y.energy
    end


    function Base.isless(x::Basics.AngularJ64, y::Basics.AngularJ64)
        return x.num < y.num
    end

        
    """
    `Basics.isSimilar()`  ... returns true if two instances are similar to each other, and false otherwise.

    + `(keya::LevelKey, keyb::LevelKey, relAcc::Float64)`  
        ... returns true if two level keys refer to the same level, i.e. level with the same symmetry and
            if the relative energy abs( (E_a - E_b)/E_a ) < relAcc. It returns false otherwise.
    """
    function Basics.isSimilar(keya::LevelKey, keyb::LevelKey, relAcc::Float64)
        ##x println("** isSimilar abs() = $(abs( (keya.energy - keyb.energy)/keya.energy ))")
        if  keya.sym == keyb.sym   &&   abs( (keya.energy - keyb.energy)/keya.energy ) < relAcc    return(true)
        else                                                                                       return(false)
        end
    end

        
    """
    `Basics.isStandardSubshellList(basis::Basis)`  
        ... returns true if the subshell list basis.subshells is in standard order, and false otherwise.
            It requests especially that both subshells of the same shell (n,l) occur in the sequence j = l-1/2, j = l+1/2
            and that l increases before n increases.
    """
    function Basics.isStandardSubshellList(basis::Basis)
        function  decimal(n,l,jnum)
            return( 10000n + 100l + jnum)
        end
        subshells = basis.subshells
        na = subshells[1].n;    la = Basics.subshell_l(subshells[1]);    ja = Basics.subshell_j(subshells[1])
        
        for  i = 2:length(subshells)
            nb = subshells[i].n;    lb = Basics.subshell_l(subshells[i]);    jb = Basics.subshell_j(subshells[1])
            if  decimal(nb,lb,jb.num) <= decimal(na,la,ja.num)    return( false )   end
        end
        
        return( true )
    end



        
    """
    Basics.merge(bases::Array{Basis,1})`  
        ... to merge two (or more) atomic bases into a single basis::Basis. This method assumes the same number of electrons in all basis and 
            that the subshell lists are the same or can be made 'consistent' to each other. Two bases have a consistent subshell list if all 
            subshells, that appear in any of the two lists appear always in the same sequence (if they are not missing at all). In the merged 
            basis, the radial orbitals are taken from the basis (in the bases-array}, from where they are found first.
    """
    function Basics.merge(bases::Array{Basis,1})

        if  length(bases) > 2 
            println("Number of bases = $(length(bases)) ...")   
            bs       = Basics.merge([ bases[1], bases[2] ])
            basesNew = [bs]
            for i = 3:length(bases)    push!( basesNew, bases[i])    end
            Basics.merge(basesNew)
        elseif  length(bases) == 1    return( basis[1] )
        elseif  length(bases) == 2
            # This is the essential step here; first create a consistent subshell list
            basisA = bases[1];    basisB = bases[2]
            wa = deepcopy( basisA.subshells );   wb = [ 100i  for i = 1:length(wa) ]
            # Add the subshells of basisB
            for  i = 1:length( basisB.subshells )
                if    basisB.subshells[i] in wa
                else  push!( wa, basisB.subshells[i] );    push!( wb, wb[i-1]+1 )   
                end
            end
            # Now sort wb and create a new subshell list from wa due to this sorting
            wc = sortperm( wb )
            newSubshells = Subshell[]
            for  i in wc   push!( newSubshells, wa[i] )    end
            println("Sorted subshell list: $(newSubshells) ")
            error("Not yet fully implemented, see source code ")
            # Transform all CSF from basisA and basisB to newSubshells ... and apply CsfRExcludeDouble()
            # return new basis
            #
        else    error("stop a")
        end
    end


    """
    `Basics.merge(multiplets::Array{Multiplet,1})`  
        ... to merge two (or more) atomic multiplets into a single multiplet::Multiplet. This method assumes (and checks) that all 
            levels have level.hasStateRep = true and that all levels refer to the same basis.
    """
    function Basics.merge(multiplets::Array{Multiplet,1})

    if      length(multiplets) == 1    return( multiplets[1] )
        elseif  length(multiplets) == 2
            levels = copy(multiplets[1].levels);   basis = multiplets[1].levels[1].basis
            for  lev  in  multiplets[2].levels
                if  !(lev.hasStateRep)   ||   lev.basis != basis    error("Levels of multiplets do not refer to the same basis.")   end
                push!(levels, lev)
            end
            mp = Multiplet( multiplets[1].name * "+" * multiplets[2].name, levels )
            return( mp )
        elseif  length(multiplets) > 2
            mp = Basics.merge([ multiplets[1], multiplets[2] ])
            for  k = 3:length(multiplets)
                mp = Basics.merge([ mp, multiplets[k] ])
            end
            return( mp )
        else    error("stop a")
        end
    end

    
    """
    `Basics.modifyLevelEnergies(multiplet::Multiplet)`  
        ... to shift the total energies of the atomic levels in multiplet interactively; a (new) multiplet::Multiplet is 
            returned in which the total energies are just modified accordingly. **Not yet implemented !**
    """
    function Basics.modifyLevelEnergies(multiplet::Multiplet)
        error("Not yet implemented !")

        return( nothing )
    end

    
    # using Plots
    # pyplot()

    using RecipesBase

    """
    `Basics.plot()`  ... plots various quantities, often in a new window.

    + `("radial potentials", potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N::Int64 = 0)`  
        ... to plot one or more radial potentials, and where N::Int64 describes the number of grid points to be considered.
            call:  using Plots; pyplot()    ... to access this method by plot(...)
    """
    function Basics.plot(sa::String, potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N::Int64 = 0)   
        error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...)")
    end


    @recipe function f(sa::String, potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N = 0)
        !(sa == "radial potentials")   &&   error("Unsupported keystring = $sa")
        wa = Float64[];   wc = [NaN for i=1:N];   labels = String[];   np = length(potentials)
        for  pot in potentials
            wb = wc
            nx = min(length(pot.Zr), N);   wb[1:nx] = pot.Zr[1:nx]
            append!(wa, wb)
            push!(labels, pot.name)
        end
        x = grid.r[1:N];     y = reshape(wa, (N, np))
        label --> reshape(labels, (1,np))
        x, y
    end



    """
    + `("radial orbitals: large", orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`  
        ... to plot the large component of one or more radial orbitals, and where N::Int64 describes the number of grid 
            points to be considered.
    + `("radial orbitals: small", orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`  
        ... to plot the small component of one or more radial orbitals.
    + `("radial orbitals: both",  orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`  
        ... to plot the large and small component of one or more radial orbitals.

        call:  using Plots; pyplot()    ... to access this method by plot(...)
    """
    function Basics.plot(sa::String, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
        error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...)")
    end


    @recipe function f(sa::String, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
        wa = Float64[];   wc = [NaN for i=1:N];   labels = String[];   np = length(orbitals)

        if       sa == "radial orbitals: large"
            for  orb in orbitals
                wb = wc
                nx = min(length(orb.P), N);   wb[1:nx] = orb.P[1:nx]
                append!(wa, wb)
                push!(labels, "$(orb.subshell):large")
            end
            ##x print("labels = $labels ")
            x = grid.r[1:N];     y = reshape(wa, (N, np))

        elseif   sa == "radial orbitals: small"
            for  orb in orbitals
                wb = wc
                nx = min(length(orb.Q), N);   wb[1:nx] = orb.Q[1:nx]
                append!(wa, wb)
                push!(labels, "$(orb.subshell):small")
            end
            ##x print("labels = $labels ")
            x = grid.r[1:N];     y = reshape(wa, (N, np))

        elseif   sa == "radial orbitals: both"
            for  orb in orbitals
                wb = wc;    nx = min(length(orb.P), N);   wb[1:nx] = orb.P[1:nx];    append!(wa, wb)
                wb = wc;    nx = min(length(orb.Q), N);   wb[1:nx] = orb.Q[1:nx];    append!(wa, wb)
                push!(labels, "$(orb.subshell):large");   push!(labels, "$(orb.subshell):small")
            end
            ##x print("labels = $labels ")
            x = grid.r[1:N];     y = reshape(wa, (N, 2np))

        else   error("Unsupported keystring = $sa") 
        end
        label --> permutedims(labels)
        x, y
    end



    """
    + `("spectrum: transition rates over energy", lines::Array{PhotoEmission.Line,1})`  
        ... to plot the transition rates of all lines as function of their transition energies. The plot is shown in a new 
            window but nothing is returned otherwise. **Not yet implemented !**

    + `("spectrum: oscillator strength over energy, emission", lines::Array{PhotoEmission.Line,1})` or
        `("spectrum: oscillator strength over energy, absorption", lines::Array{PhotoEmission.Line,1})` 
        ... to plot the absorption oscillator strength of all lines as function of their transition energies. Again, a new 
            window is opened but nothing returned by this method. **Not yet implemented !**
    """
    function Basics.plot(sa::String, lines::Array{PhotoEmission.Line,1})
        error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...) ... not yet implemented !")
    end



    """
    + `("spectrum: transition rates over energy, Gaussian", lines::Array{PhotoEmission.Line,1}; widths=value::Float64)` or
        `("spectrum: transition rates over energy, Lorentzian", lines::Array{PhotoEmission.Line,1}; widths=value::Float64)`
        ... to plot the transition rates of all lines as function of their transition energies but with a Gaussian or Lorentzian 
            distribution. Again, a new window is opened but nothing returned by this method. It still need to be decided how 
            the widths (and, perhaps, other parameters) will be communicated to the method. **Not yet implemented !**
    """
    function Basics.plot(sa::String, lines::Array{PhotoEmission.Line,1}, widths::Float64)
        error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...) ... not yet implemented !")
    end


end # module