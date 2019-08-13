
"""
`module  JAC.ManyElectronAZ`  
    ... a submodel of JAC that contains methods that support tasks related to spectroscopic computation.
"""
module ManyElectronAZ

    using Printf, ..Basics, ..Defaults, ..ManyElectron, ..Radial
    
    export provideSubshellStates


    """
    `ManyElectron.provideSubshellStates(sh::Subshell, occ::Int64)`   
        ... to provide all antisymmetric SubshellStates within the relativistic seniority scheme for the given 
            shell and occupation; an Array{SubshellStateR,1} is returned.
    """
    function ManyElectron.provideSubshellStates(sh::Subshell, occ::Int64)

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

end # module
