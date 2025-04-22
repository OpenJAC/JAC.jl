
using Printf, ..Basics, ..Defaults, ..Nuclear, ..ManyElectron, ..Radial, ..TableStrings

#################################################################################################################################
#################################################################################################################################

## Photoionization plasma rates and rate coefficients

"""
`struct  Plasma.RateSchemePhotoionization  <:  Plasma.AbstractPlasmaScheme`  
    ... a struct to perform the computation of empirical configuration-averaged photoionization rates and rate coefficients.

    + calcRates             ::Bool                 ... True, if photoionization plasma rates are to be calculated.
    + calcAlphas            ::Bool                 ... True, if photoionization plasma rate coefficents are to be calculated.
    + calcCrossSections     ::Bool                 ... True, if photoionization cross sections are to be calculated.
    + printExplanations     ::Bool                 ... True, if additional explanations about atomic data need to be printed.
    + photonDistribution    ::Bool                 ... ::Basics.AbstractPhotonDistribution
    + calcTotalRates        ::Bool                 ... True, if total rates are calculated for each configuation-averaged level
    + calcPairwiseRates     ::Bool                 ... True, if rates are calculated for each pair of configuation-averaged levels
"""
struct   RateSchemePhotoionization  <:  Plasma.AbstractPlasmaScheme
    calcRates               ::Bool
    calcAlphas              ::Bool
    calcCrossSections       ::Bool 
    printExplanations       ::Bool 
    photonDistribution      ::Bool    
end  


"""
`Plasma.RateSchemePhotoionization()`  ... constructor for an 'default' instance of a Plasma.RateSchemePhotoionization.
"""
function RateSchemePhotoionization()
    RateSchemePhotoionization( false, false, false, false, false )
end


# `Base.show(io::IO, scheme::RateSchemePhotoionization)`  ... prepares a proper printout of the scheme::RateSchemePhotoionization.
function Base.show(io::IO, scheme::RateSchemePhotoionization)
    println(io, "calcRates:             $(scheme.calcRates)  ")
    println(io, "calcAlphas:            $(scheme.calcAlphas)  ")
    println(io, "calcCrossSections:     $(scheme.calcCrossSections)  ")
    println(io, "printExplanations:     $(scheme.printExplanations)  ")
end


"""
`Plasma.computeRatesEtc(scheme::RateSchemePhotoionization, refConfigs::Array{Configuration,1}, nm::Nuclear.Model, 
                        grid::Radial.Grid, asfSettings::AsfSettings, settings::Plasma.Settings; output=true)` 
    ... to compute ... nothing::Nothing is returned.
"""
function computeRatesEtc(scheme::RateSchemePhotoionization, refConfigs::Array{Configuration,1}, nm::Nuclear.Model, 
                         grid::Radial.Grid, asfSettings::AsfSettings, settings::Plasma.Settings; output=true)
    println("")
    printstyled("Plasma.computeRatesEtc(): The computation of photoionization plasma-rates starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    # Compute R(T; i --> f),  R(T; i --> sum f) = R(T; i) ... + Ausgabe kommentieren mit printExplanations
    #
    return( nothing )
end
