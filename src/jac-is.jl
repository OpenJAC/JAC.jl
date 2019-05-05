



export isless

function isless(x :: JAC.ManyElectron.Level, y :: JAC.ManyElectron.Level)
    return x.energy < y.energy
end


function isless(x :: JAC.Hfs.IJF_Level, y :: JAC.Hfs.IJF_Level)
    return x.energy < y.energy
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
