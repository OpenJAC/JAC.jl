
"""
`module  JAC.Distribution`  
	... a submodel of JAC that contains all methods for computing different distributions, such as the Maxwell-Boltzmann
        (normalized) velocity and energy distributions, Planck's spectral black-body photon-energy distributions
        and various others.
"""
module Distribution


using  Printf, ..Basics, ..Defaults


"""
`abstract type Distribution.AbstractLineProfile` 
    ... defines an abstract and a number of singleton types to deal with different line profiles, for instance,
        in a plasma.

    + GaussianProfile       ... assumes a Gaussian line profile L^Gaussian (omega)
    + LorentzianProfile     ... assumes a Lorentzian line profile L^Lorentzian (omega)
"""
abstract type  AbstractLineProfile                                   end
struct     GaussianProfile               <:  AbstractLineProfile     end
struct     LorentzianProfile             <:  AbstractLineProfile     end

export  AbstractLineProfile, GaussianProfile, LorentzianProfile

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Distribution.AbstractPlasmaDistribution` 
    ... defines an abstract and a number of singleton types to deal with different electron and photon-flux
        distributions in a plasma.

    + ElectronMaxwell       ... applies the (normalized) non-relativistic Maxwell-Boltzmann electron distribution.
    + ElectronJuettner      ... applies the (normalized) non-relativistic Maxwell-Juettner velocity distribution.
    + SpectralFluxPlanck    ... applies Planck's (normalized, black-body) spectral-photon distribution.
"""
abstract type  AbstractPlasmaDistribution                            end
struct     ElectronMaxwell        <:  AbstractPlasmaDistribution     end
struct     ElectronJuettner       <:  AbstractPlasmaDistribution     end
struct     SpectralFluxPlanck     <:  AbstractPlasmaDistribution     end

export  AbstractPlasmaDistribution, ElectronMaxwell, ElectronJuettner, SpectralFluxPlanck

#################################################################################################################################
#################################################################################################################################


"""
`Distribution.electronEnergyMaxwell(energy::Float64, temp::Float64)`  
    ... to provide the probability density f(en) to find a electron with energy en in a Maxwell-Boltzmann
        distributed plasma. A value::Float64 is returned. 
"""
function  electronEnergyMaxwell(energy::Float64, temp::Float64)
    error("Not yet implemented.")
    f = 0.
    return( f )
end


"""
`Distribution.electronVelocityMaxwell(velocity::Float64, temp::Float64)`  
    ... to provide the probability density f(velocity) to find a electron with velocity v in a Maxwell-Boltzmann
        distributed plasma. A value::Float64 is returned. 
"""
function  electronVelocityMaxwell(velocity::Float64, temp::Float64)
    error("Not yet implemented.")
    f = 0.
    return( f )
end


"""
`Distribution.electronVelocityJuettner(velocity::Float64, temp::Float64)`  
    ... to provide the probability density f(velocity) to find a electron with velocity v in a Maxwell-Juettner
        distributed plasma. A value::Float64 is returned. 
"""
function  electronVelocityJuettner(velocity::Float64, temp::Float64)
    error("Not yet implemented.")
    f = 0.
    return( f )
end


end # module
