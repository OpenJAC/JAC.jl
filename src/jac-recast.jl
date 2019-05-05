
"""
`Basics.recast()`  
    ... recasts some data from one number/representation into another one; cf. the supported keystrings and return values.

  + `("rate: radiative, to decay width", value::Float64)`  
    ... to recast a given radiative rate (Einstein A in atomic units) into a decay withs, taking the selected energy unit 
        into account. A Float64 is returned.

  + `("rate: radiative, to Einstein A", value::Float64)`  
    ... to recast a given spontaneous radiative rate (= Einstein A-coefficient), taking the selected unit into account. 
        A Float64 is returned.

  + `("rate: radiative, to Einstein B", value::Float64)`  
    ... to recast a given radiative rate (Einstein A in atomic units) into an Einstein B-coefficient, taking the selected 
        unit into account. A Float64 is returned.

  + `("rate: radiative, to g_f", value::Float64)`  
    ... to recast a given radiative rate (Einstein A in atomic units) into an oscillator strength g_f; 
        a Float64 is returned.
"""
function Basics.recast(sa::String, line::Union{Einstein.Line,PhotoEmission.Line}, wa::Float64)
    ##x global  CONVERT_ENERGY_AU_TO_EV, CONVERT_ENERGY_AU_TO_KAYSERS, CONVERT_ENERGY_AU_TO_PER_SEC, CONVERT_TIME_AU_TO_SEC,
    ##x         CONVERT_CROSS_SECTION_AU_TO_BARN,  CONVERT_RATE_AU_TO_PER_SEC 
   
    if       sa == "rate: radiative, to decay width"
        width = Defaults.convertUnits("energy: from atomic", wa)
        return( width )

    elseif   sa == "rate: radiative, to Einstein A"
        einsteinA = Defaults.convertUnits("rate: from atomic", wa)
        return( einsteinA )

    elseif   sa == "rate: radiative, to Einstein B"
        einsteinB = (JAC.AngularMomentum.twoJ(line.initialLevel.J) + 1) / (JAC.AngularMomentum.twoJ(line.finalLevel.J) + 1) * pi^2 *
                    Defaults.getDefaults("speed of light: c")^3 / line.omega^3  * wa
        einsteinB = Defaults.convertUnits("Einstein B: from atomic", einsteinB)
        return( einsteinB )

    elseif   sa == "rate: radiative, to g_f"
        gf = (JAC.AngularMomentum.twoJ(line.initialLevel.J) + 1) / (JAC.AngularMomentum.twoJ(line.finalLevel.J) + 1) / 2. * 
              Defaults.getDefaults("speed of light: c")^3 / line.omega^2 * wa     
        return( gf )

    else     error("Unsupported keystring = $sa")
    end
end

