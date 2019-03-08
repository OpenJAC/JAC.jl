
## export  merge

"""
`JAC.merge()`  ... to merge data from different instances into a single one.

  + `("atomic bases", bases::Array{Basis,1})`  ... to merge two (or more) atomic bases into a single basis::Basis. This method assumes the same 
                      number of electrons in all basis and that the subshell lists are the same or can be made 'consistent' to each other. 
                      Two bases have a consistent subshell list if all subshells, that appear in any of the two lists appear always in the same 
                      sequence (if they are not missing at all). In the merged basis, the radial orbitals are taken from the basis (in the
                      bases-array}, from where they are found first.
"""
function merge(sa::String, bases::Array{Basis,1})
    sa != "atomic bases"   &&   error("Unsupported keystring = $sa")

    if  length(bases) > 2 
        println("Number of bases = $(length(bases)) ...")   
        bs       = JAC.merge("atomic bases", [ bases[1], bases[2] ])
        basesNew = [bs]
        for i = 3:length(bases)    push!( basesNew, bases[i])    end
        JAC.merge("atomic bases", basesNew )
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
  + `("multiplets", multiplets::Array{Multiplet,1})`  ... to merge two (or more) atomic multiplets into a single multiplet::Multiplet.
                    This method assumes (and checks) that all levels have level.hasStateRep = true and that all levels
                    refer to the same basis.
"""
function merge(sa::String, multiplets::Array{Multiplet,1})
    sa != "multiplets"   &&   error("Unsupported keystring = $sa")
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
        mp = JAC.merge("multiplets", [ multiplets[1], multiplets[2] ])
        for  k = 3:length(multiplets)
            mp = JAC.merge("multiplets", [ mp, multiplets[k] ])
        end
        return( mp )
    else    error("stop a")
    end
end

