
"""
`Basics.estimate()`  
    ... estimates various quantities, typically, due to semi-empirical formulas and parametrizations.

  + `("ionization potential: inner-shell", shell::Shell, Z::Int64)`  
    ... to estimate the ionization potential = mean binding energy of an electron in the given shell; an energy value::Float64 
        is returned. These ionization potentials are taken from the tabulation ....
"""
function Basics.estimate(sa::String, shell::Shell, Z::Int64)
    if      sa == "ionization potential: inner-shell"
        shellIndex = Dict( Subshell("1s_1/2") => 1, Subshell("2s_1/2") => 2, Subshell("2p_1/2") => 3, Subshell("2p_3/2") => 4, 
                           Subshell("3s_1/2") => 5, Subshell("3p_1/2") => 6, Subshell("3p_3/2") => 7, 
                           Subshell("3d_3/2") => 8, Subshell("3d_3/2") => 9  )
        idx = shellIndex[shell]
        return( JAC.store("binding energies: Williams (2000)", Z)[idx] )
    else
        error("Unsupported keystring = $sa")
    end
end


