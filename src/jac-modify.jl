
export  modify

"""
`JAC.modify()`  ... modifies some data, either interactively or due to some given rule.

  + `("level energies: interactive", multiplet::Multiplet)`  ... to shift the total energies of the atomic levels in multiplet interactively; 
                       a (new) multiplet::Multiplet is returned in which the total energies are just modified accordingly.  
                       **Not yet implemented !**
"""
function modify(sa::String, multiplet::Multiplet)
    !(sa == "level energies: interactive")   &&   error("Unsupported keystring = $sa")
    error("Not yet implemented !")

    return( nothing )
end

