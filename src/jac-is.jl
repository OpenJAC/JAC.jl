



export isless

function isless(x :: JAC.ManyElectron.Level, y :: JAC.ManyElectron.Level)
    return x.energy < y.energy
end


function isless(x :: JAC.Hfs.IJF_Level, y :: JAC.Hfs.IJF_Level)
    return x.energy < y.energy
end


function isless(x :: JAC.Basics.AngularJ64, y :: JAC.Basics.AngularJ64)
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



