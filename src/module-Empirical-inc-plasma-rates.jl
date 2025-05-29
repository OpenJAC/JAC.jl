
#################################################################################################################################
### Photoemission (PE) ##########################################################################################################

"""
`Empirical.photoemissionSpontaneousRate(iConf::Configuration, fConf::Configuration; printout::Bool=false)`  
    ... to estimate the photoemission rate for a transition from iConf -> fConf. A rate::Float64 is returned. 
"""
function photoemissionSpontaneousRate(iConf::Configuration, fConf::Configuration; printout::Bool=false)
    error("Not yet implemented.")
    rate = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( rate )
end


#################################################################################################################################
### Photoexcitation (PX) ########################################################################################################


"""
`Empirical.photoexcitationCrossSection(iConf::Configuration, fConf::Configuration; 
                                       approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), printout::Bool=false)`  
    ... to estimate the photoexcitation cross section for a transition from iConf -> fConf by means of
        Kramer (1923) approximation. A cs::Float64 is returned. 
"""
function photoexcitationCrossSection(iConf::Configuration, fConf::Configuration; 
                                     approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), printout::Bool=false)
    error("Not yet implemented.")
    cs = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( cs )
end


"""
`Empirical.photoexcitationPlasmaRate(iConf::Configuration, fConf::Configuration; 
                                     flux::Distribution.AbstractPlasmaDistribution=SpectralFluxPlanck(),
                                     profile::Distribution.AbstractLineProfile=GaussianProfile(), 
                                     printout::Bool=false)`  
    ... to estimate the photoexcitation plasma for a transition from iConf -> fConf by applying a given
        photon-flux distribution and line profile. A rate::Float64 is returned. 
"""
function photoexcitationPlasmaRate(iConf::Configuration, fConf::Configuration; 
                                   flux::Distribution.AbstractPlasmaDistribution=SpectralFluxPlanck(),
                                   profile::Distribution.AbstractLineProfile=GaussianProfile(), 
                                   printout::Bool=false)
    error("Not yet implemented.")
    rate = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( rate )
end


#################################################################################################################################
### Photoionization (PI) ########################################################################################################


"""
`Empirical.photoionizationCrossSection(omega::Float64, iConf::Configuration, fConf::Configuration; 
                                       approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), printout::Bool=false)`  
    ... to estimate the photoionization cross section for a transition from iConf --> fConf by means of
        Kramer (1923) approximation. A cs::Float64 is returned. 
"""
function photoexcitationCrossSection(omega::Float64, iConf::Configuration, fConf::Configuration; 
                                     approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), printout::Bool=false)
    error("Not yet implemented.")
    cs = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( cs )
end


"""
`Empirical.photoionizationPlasmaRate(iConf::Configuration, fConf::Configuration; 
                                     flux::Distribution.AbstractPlasmaDistribution=SpectralFluxPlanck(),
                                     approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), 
                                     printout::Bool=false)`  
    ... to estimate the photoexcitation plasma for a transition from iConf -> fConf by applying a given
        photon-flux distribution and line profile. A rate::Float64 is returned. 
"""
function photoionizationPlasmaRate(iConf::Configuration, fConf::Configuration; 
                                   flux::Distribution.AbstractPlasmaDistribution=SpectralFluxPlanck(),
                                   approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), 
                                   printout::Bool=false)
    error("Not yet implemented.")
    rate = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( rate )
end


#################################################################################################################################
### Photorecombination (PR) ################################################################################################


"""
`Empirical.photorecombinationPlasmaRateAlpha(temp::Float64, iConf::Configuration; 
                                     approx::Empirical.AbstractEmpiricalApproximation=Axelrod1980(), 
                                     printout::Bool=false)`  
    ... to estimate the photorecombination plasma rate coefficient alpha for an ion in the initial level iConf 
        by applying a given approximation. An alpha::Float64 is returned. 
"""
function photorecombinationPlasmaRateAlpha(temp::Float64, iConf::Configuration; 
                                           approx::Empirical.AbstractEmpiricalApproximation=Axelrod1980(), 
                                           printout::Bool=false)
    error("Not yet implemented.")
    alpha = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( alpha )
end



#################################################################################################################################
### Three-Body Recombination (TBR) ##############################################################################################


"""
`Empirical.threeBRCrossSection(iConf::Configuration, fConf::Configuration; 
                               approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), printout::Bool=false)`  
    ... to estimate the TBR cross section for a transition from iConf --> fConf by means of
        Kramer (1923) approximation. A cs::Float64 is returned. 
"""
function threeBRCrossSection(iConf::Configuration, fConf::Configuration; 
                             approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), printout::Bool=false)
    error("Not yet implemented.")
    cs = 0.
    return( cs )
end


"""
`Empirical.threeBRPlasmaRateAlpha(temp::Float64, iConf::Configuration; 
                                  approx::Empirical.AbstractEmpiricalApproximation=Axelrod1980(), 
                                  printout::Bool=false)`  
    ... to estimate the 3BR plasma rate coefficient alpha for an ion in the initial level iConf 
        by applying a given approximation. An alpha::Float64 is returned. 
"""
function threeBRPlasmaRateAlpha(temp::Float64, iConf::Configuration; 
                                approx::Empirical.AbstractEmpiricalApproximation=Axelrod1980(), 
                                printout::Bool=false)
    error("Not yet implemented.")
    alpha = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( alpha )
end


#################################################################################################################################
### Autoionization (AI) #########################################################################################################


"""
`Empirical.autoionizationPlasmaRate(iConf::Configuration, fConf::Configuration; 
                                    flux::Distribution.AbstractPlasmaDistribution=SpectralFluxPlanck(),
                                    approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), 
                                    printout::Bool=false)`  
    ... to estimate the autoionization plasma rate for a transition from iConf -> fConf by applying a given
        photon-flux distribution. A rate::Float64 is returned. 
"""
function autoionizationPlasmaRate(iConf::Configuration, fConf::Configuration; 
                                  flux::Distribution.AbstractPlasmaDistribution=SpectralFluxPlanck(),
                                  approx::Empirical.AbstractEmpiricalApproximation=Kramers1923(), 
                                  printout::Bool=false)
    error("Not yet implemented.")
    rate = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( rate )
end


#################################################################################################################################
### Electron-impact excitation (EIE) ############################################################################################

"""
`Empirical.eieCollisionStrength(energy::Float64, iConf::Configuration, fConf::Configuration; 
                               approx::Empirical.AbstractEmpiricalApproximation=Regemorter1962(), 
                               printout::Bool=false)`  
    ... to estimate the EIE collision strength for a transition from iConf -> fConf by applying a given
        approximation. A cs::Float64 is returned. 
"""
function eieCollisionStrength(energy::Float64, iConf::Configuration, fConf::Configuration; 
                              approx::Empirical.AbstractEmpiricalApproximation=Regemorter1962(), 
                              printout::Bool=false)
    error("Not yet implemented.")
    cs = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( cs )
end


"""
`Empirical.eieCrossSection(energy::Float64, iConf::Configuration, fConf::Configuration; 
                           approx::Empirical.AbstractEmpiricalApproximation=Regemorter1962(), 
                           printout::Bool=false)`  
    ... to estimate the EIE cross section for a transition from iConf -> fConf by applying a given
        approximation. A cs::Float64 is returned. 
"""
function eieCrossSection(energy::Float64, iConf::Configuration, fConf::Configuration; 
                         approx::Empirical.AbstractEmpiricalApproximation=Regemorter1962(), 
                         printout::Bool=false)
    error("Not yet implemented.")
    cs = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( cs )
end


"""
`Empirical.eiePlasmaRateAlpha(temp::Float64, iConf::Configuration; 
                              approx::Empirical.AbstractEmpiricalApproximation=Axelrod1980(), 
                              printout::Bool=false)`  
    ... to estimate the EIE plasma rate coefficient alpha for an ion in the initial level iConf 
        by applying a given approximation. An alpha::Float64 is returned. 
"""
function eiePlasmaRateAlpha(temp::Float64, iConf::Configuration; 
                            approx::Empirical.AbstractEmpiricalApproximation=Axelrod1980(), 
                            printout::Bool=false)
    error("Not yet implemented.")
    alpha = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( alpha )
end


#################################################################################################################################
### Electron-impact ionization (EII) ############################################################################################

"""
`Empirical.eiiCollisionStrength(energy::Float64, iConf::Configuration, fConf::Configuration; 
                               approx::Empirical.AbstractEmpiricalApproximation=Regemorter1962(), 
                               printout::Bool=false)`  
    ... to estimate the EII collision strength for a transition from iConf -> fConf by applying a given
        approximation. A cs::Float64 is returned. 
"""
function eiiCollisionStrength(energy::Float64, iConf::Configuration, fConf::Configuration; 
                              approx::Empirical.AbstractEmpiricalApproximation=Regemorter1962(), 
                              printout::Bool=false)
    error("Not yet implemented.")
    cs = 0.
    
    if  printout
        unit = "aaa"
        wb   = 3.0 * cs
        println("> Stopping power due to Bohr (1913):  -dE / dx ($energy, $mRatio, $omegaPlasma´ = $wb  $unit" *
                "\n\n> ...")
    end
    
    return( cs )
end
