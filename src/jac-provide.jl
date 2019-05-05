
"""
`Basics.provide()`  ... provides various data and transformations. 

  + `("subshell states: antisymmetric, seniority", sh::Subshell, occ::Int64)`   
    ... to provide all antisymmetric SubshellStates within the relativistic seniority scheme for the given 
        shell and occupation; an Array{SubshellStateR,1} is returned.
"""
function Basics.provide(sa::String, sh::Subshell, occ::Int64)

    !(sa == "subshell states: antisymmetric, seniority")   &&   error("Unsupported keystring = $sa")

    # Use the j-value = wa[4] of the subshell to select the allowed SubshellStateR's
    wb = Basics.SubshellStateR[]

    if       Basics.subshell_2j(sh) == 1
        if       occ == 0  || occ == 2      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );       return(wb)
        elseif   occ == 1                   push!( wb, Basics.SubshellStateR( sh, occ, 1, 1) );       return(wb)
        end
    elseif   Basics.subshell_2j(sh) == 3
        if       occ == 0  || occ == 4      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );       return(wb)
        elseif   occ == 1  || occ == 3      push!( wb, Basics.SubshellStateR( sh, occ, 1, 3) );       return(wb)
        elseif   occ == 2                   push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 4) );       return(wb)
        end
    elseif   Basics.subshell_2j(sh) == 5
        if       occ == 0  || occ == 6      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );       return(wb)
        elseif   occ == 1  || occ == 5      push!( wb, Basics.SubshellStateR( sh, occ, 1, 5) );       return(wb)
        elseif   occ == 2  || occ == 4      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 4) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 8) );       return(wb)
        elseif   occ == 3                   push!( wb, Basics.SubshellStateR( sh, occ, 1, 5) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 3) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 9) );       return(wb)
        end
    elseif   Basics.subshell_2j(sh) == 7
        if       occ == 0  || occ == 8      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );       return(wb)
        elseif   occ == 1  || occ == 7      push!( wb, Basics.SubshellStateR( sh, occ, 1, 7) );       return(wb)
        elseif   occ == 2  || occ == 6      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 4) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 8) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2,12) );       return(wb)
        elseif   occ == 3  || occ == 5      push!( wb, Basics.SubshellStateR( sh, occ, 1, 7) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 3) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 5) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 9) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3,11) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3,15) );       return(wb)
        elseif   occ == 4                   push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 4) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 8) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2,12) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 4, 4) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 4, 8) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 4,10) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 4,16) );       return(wb)
        end
    elseif   Basics.subshell_2j(sh) in [ 9, 11, 13, 15]
        if       occ == 0                   push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );                 return(wb)
        elseif   occ == 1                   push!( wb, Basics.SubshellStateR( sh, occ, 1, subshell_2j(sh)) );   return(wb)
        else     error("stop b")
        end

    # For an extension to other subshell-occupations, see the Racah-manual; Table ...
    else  error("stop a")
    end
end


"""
  + `("binding energy", Z::Float64, sh::Subshell)`   
    ... to provide the binding energy of a subshell electron, taken from a semi-empirical tabulations by Williams et al., 
        https://userweb.jlab.org/~gwyn/ebindene.html. A energy::Float64 in  Hartree is returned.
"""
function Basics.provide(sa::String, Z::Int64, sh::Subshell)
    if     sa == "binding energy"
        wa = JAC.store_Williams2000(Z)
        if      sh == Subshell("1s_1/2")    wb = wa[1]
        elseif  sh == Subshell("2s_1/2")    wb = wa[2]
        elseif  sh == Subshell("2p_1/2")    wb = wa[3]
        elseif  sh == Subshell("2p_3/2")    wb = wa[4]
        elseif  sh == Subshell("3s_1/2")    wb = wa[5]
        elseif  sh == Subshell("3p_1/2")    wb = wa[6]
        elseif  sh == Subshell("3p_3/2")    wb = wa[7]
        elseif  sh == Subshell("3d_3/2")    wb = wa[8]
        elseif  sh == Subshell("3d_5/2")    wb = wa[9]
        elseif  sh == Subshell("4s_1/2")    wb = wa[10]
        elseif  sh == Subshell("4p_1/2")    wb = wa[11]
        elseif  sh == Subshell("4p_3/2")    wb = wa[12]
        else    error("No binding energy available for Z = $Z and subshell $sh ")
        end
    else   error("Unsupported keystring")
    end
    #
    if     wb == -1.   error("No binding energy available for Z = $Z and subshell $sh ")
    else   wb = Defaults.convertUnits("energy: from eV to atomic", wb)
    end
    # 
    return( wb )
end



"""
  + `("binding energy", Z::Float64, conf::Configuration)`  
    ... to provide an approximate binding energy of a given electron configuration. This estimate adds the binding 
        energies of all subshell, taken frogm a semi-empirical tabulations by Williams et al., 
        https://userweb.jlab.org/~gwyn/ebindene.html. No relaxation effects are included if several hole states
        occur with regard to the neutral atom. An energy::Float64 in  Hartree is returned.
"""
function Basics.provide(sa::String, Z::Int64, conf::Configuration)
    if     sa == "binding energy"
        wa = JAC.PeriodicTable.bindingEnergies_Williams2000(Z);    wb = 0.
        for (sh,v) in  conf.shells
            if      sh == Shell("1s")    wb = wb + v * wa[1];    if wa[1]  == -1.   error("stop aa")   end
            elseif  sh == Shell("2s")    wb = wb + v * wa[2];    if wa[2]  == -1.   error("stop ab")   end
            elseif  sh == Shell("2p")    wb = wb + v * wa[3];    if wa[3]  == -1.   error("stop ac")   end
            elseif  sh == Shell("3s")    wb = wb + v * wa[5];    if wa[5]  == -1.   error("stop ae")   end
            elseif  sh == Shell("3p")    wb = wb + v * wa[6];    if wa[6]  == -1.   error("stop af")   end
            elseif  sh == Shell("3d")    wb = wb + v * wa[8];    if wa[8]  == -1.   error("stop ah")   end
            elseif  sh == Shell("4s")    wb = wb + v * wa[10];   if wa[10] == -1.   error("stop aj")   end
            elseif  sh == Shell("4p")    wb = wb + v * wa[11];   if wa[11] == -1.   error("stop ak")   end
            else    error("No binding energy available for Z = $Z and subshell $sh ")
            end
        end
    else   error("Unsupported keystring")
    end
    #
    wb = Defaults.convertUnits("energy: from eV to atomic", wb)
    # 
    return( wb )
end

