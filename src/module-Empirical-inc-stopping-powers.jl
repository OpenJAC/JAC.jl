

"""
`Empirical.stoppingPowerElectronBohr1913(energy::Float64, nRatio::Float64, omegaPlasma::Float64; printout::Bool=false)`  
    ... to evaluate the stopping power - d E / dx of an electron with kinetic energy in a plasma with
        nRatio = n_electrons / n_ions and for the plasma frequency omegaPlasma. A sPower::Float64 is returned.
        If printout = true, the value is printed in user-defined values along with short other explanations.
"""
function stoppingPowerElectronBohr1913(energy::Float64, nRatio::Float64, omegaPlasma::Float64; printout::Bool=false)
    error("Not yet implemented.")
    wa = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * wa
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( wa )
end


"""
`Empirical.stoppingPowerElectronBethe1931(energy::Float64, nRatio::Float64, omegaPlasma::Float64; printout::Bool=false)`  
    ... to evaluate the stopping power - d E / dx of an electron with kinetic energy in a plasma with
        nRatio = n_electrons / n_ions and for the plasma frequency omegaPlasma. A sPower::Float64 is returned.
        If printout = true, the value is printed in user-defined values along with short other explanations.
"""
function stoppingPowerElectronBethe1931(energy::Float64, nRatio::Float64, omegaPlasma::Float64; printout::Bool=false)
    error("Not yet implemented.")
    wa = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * wa
        println("> Stopping power due to Bethe (1931):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( wa )
end
