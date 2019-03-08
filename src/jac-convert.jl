
# export  convert

"""
`JAC.convert()`  ... converts some data from one format/unit into another one; cf. the supported keystrings and return values.

  + `("cross section: from atomic to predefined unit", value::Float64)`  or  `("energy: from atomic", value::Float64)`  ... to convert 
                                                an energy value from atomic to the predefined cross section unit; a Float64 is returned.

  + `("cross section: from atomic to barn", value::Float64)`  or  `("cross section: from atomic to Mbarn", value::Float64)`  or
    `("cross section: from atomic to Hz", value::Float64)`  or  `("energy: from atomic to Angstrom", value::Float64)` ... to convert an 
                                                energy value from atomic to the speficied cross section unit; a Float64 is returned.

  + `("cross section: from predefined to atomic unit", value::Float64)`  or  `("cross section: to atomic", value::Float64)`... to convert a
                                                cross section value from the predefined to the atomic cross section unit; a Float64 is returned.

  + `("einstein B: from atomic", value::Float64)`  ... to convert a Einstein B coefficient from atomic to the speficied energy units; 
                                                       a Float64 is returned.or

  + `("energy: from atomic to eV", value::Float64)`  or  `("energy: from atomic to Kayser", value::Float64)`  or
    `("energy: from atomic to Hz", value::Float64)`  or  `("energy: from atomic to Angstrom", value::Float64)` ... to convert an energy value 
                                                from atomic to the speficied energy unit; a Float64 is returned.

  + `("energy: from predefined to atomic unit", value::Float64)`  or  `("energy: to atomic", value::Float64)`... to convert an energy value 
                                                from the predefined to the atomic energy unit; a Float64 is returned.
  + `("energy: from eV to atomic unit", value::Float64)` ... to convert an energy value from eV to the atomic energy unit; a Float64 is returned.

  + `("kinetic energy to wave number: atomic units", value::Float64)`  ... to convert a kinetic energy value (in a.u.) into a wave number
                                                k (a.u.); a Float64 is returned.

  + `("kinetic energy to wavelength: atomic units", value::Float64)`  ... to convert a kinetic energy value (in a.u.) into a wavelength (a.u.); 
                                                a Float64 is returned.

  + `("length: from fm to atomic", value::Float64)`  ... to convert a length value (in fm) into a.u.;  a Float64 is returned.
  + `("length: from atomic to fm", value::Float64)`  ... to convert a length value (in Bohr's a.u.) into fm;  a Float64 is returned.

  + `("rate: from atomic to predefined unit", value::Float64)`  or  ("rate: from atomic", value::Float64)  ... to convert a rate value 
                                              from atomic to the predefined rate unit; a Float64 is returned.

  + `("rate: from atomic to 1/s", value::Float64)`  ... to convert an rate value from atomic to the speficied rate unit; a Float64 is returned.

  + `("rate: from predefined to atomic unit", value::Float64)`  or  `("rate: to atomic", value::Float64)'... to convert a
                                                rate value from the predefined to the atomic rate unit; a Float64 is returned.

  + `("time: from atomic to predefined unit", value::Float64)`  or  ("time: from atomic", value::Float64)  ... to convert an time value 
                                              from atomic to the predefined time unit; a Float64 is returned.

  + `("time: from atomic to sec", value::Float64)`  or  `("time: from atomic to fs", value::Float64)'  or
    `("time: from atomic to as", value::Float64)`  ... to convert a time value from atomic to the speficied time unit; a Float64 is returned.

  + `("time: from predefined to atomic unit", value::Float64)`  or  `("time: to atomic", value::Float64)'  ... to convert a
                                              time value from the predefined to the atomic time unit; a Float64 is returned.

  + `("wave number to total electron energy: atomic units", value::Float64)`  ... to convert a wavenumber (a.u.) into the total electron
                                              energy, including the rest energy; a Float64 is returned.

  + `("wave number to kinetic energy: atomic units", value::Float64)`  ... to convert a wavenumber (a.u.) into the kinetic energy; 
                                              a Float64 is returned.
"""
function convert(sa::String, wa::Float64)
    global  CONVERT_ENERGY_AU_TO_EV, CONVERT_ENERGY_AU_TO_KAYSERS, CONVERT_ENERGY_AU_TO_PER_SEC, CONVERT_TIME_AU_TO_SEC,
            CONVERT_CROSS_SECTION_AU_TO_BARN,  CONVERT_RATE_AU_TO_PER_SEC,  CONVERT_LENGTH_AU_TO_FEMTOMETER 
   
    if       sa in ["cross section: from atomic to predefined unit", "cross section: from atomic"]
        if      JAC.give("unit: cross section") == "a.u."   return( wa )
        elseif  JAC.give("unit: cross section") == "barn"   return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN )
        elseif  JAC.give("unit: cross section") == "Mbarn"  return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN * 10.0e-6 )
        else    error("stop a")
        end
      
    elseif   sa in ["cross section: from atomic to barn"]   return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN )
    elseif   sa in ["cross section: from atomic to Mbarn"]  return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN * 10.0e-6 )

    elseif   sa in ["cross section: from predefined to atomic unit", "cross section: to atomic"]
        if      JAC.give("unit: cross section") == "a.u."   return( wa )
        elseif  JAC.give("unit: cross section") == "barn"   return( wa / CONVERT_CROSS_SECTION_AU_TO_BARN )
        elseif  JAC.give("unit: cross section") == "Mbarn"  return( wa / CONVERT_CROSS_SECTION_AU_TO_BARN * 10.0e6 )
        else    error("stop b")
        end

    elseif   sa in ["Einstein B: from atomic"]              return( wa )  # Einstein B not yet properly converted.

    elseif   sa in ["energy: from atomic to predefined unit", "energy: from atomic"]
        if      JAC.give("unit: energy") == "eV"            return( wa * CONVERT_ENERGY_AU_TO_EV )
        elseif  JAC.give("unit: energy") == "Kayser"        return( wa * CONVERT_ENERGY_AU_TO_KAYSERS )
        elseif  JAC.give("unit: energy") == "Hartree"       return( wa )
        elseif  JAC.give("unit: energy") == "Hz"            return( wa * CONVERT_ENERGY_AU_TO_PER_SEC )
        elseif  JAC.give("unit: energy") == "A"             return( 1.0e8 / ( CONVERT_ENERGY_AU_TO_KAYSERS * wa ) )
        else    error("stop c")
        end
      
    elseif   sa in ["energy: from atomic to eV"]            return( wa * CONVERT_ENERGY_AU_TO_EV )
    elseif   sa in ["energy: from atomic to Kayser"]        return( wa * CONVERT_ENERGY_AU_TO_KAYSERS )
    elseif   sa in ["energy: from atomic to Hz"]            return( wa * CONVERT_ENERGY_AU_TO_PER_SEC )
    elseif   sa in ["energy: from atomic to Angstrom"]      return( 1.0e8 / ( CONVERT_ENERGY_AU_TO_KAYSERS * wa ) )

    elseif   sa in ["energy: from predefined to atomic unit", "energy: to atomic"]
        if      JAC.give("unit: energy") == "eV"            return( wa / CONVERT_ENERGY_AU_TO_EV )
        elseif  JAC.give("unit: energy") == "Kayser"        return( wa / CONVERT_ENERGY_AU_TO_KAYSERS )
        elseif  JAC.give("unit: energy") == "Hartree"       return( wa )
        elseif  JAC.give("unit: energy") == "Hz"            return( wa / CONVERT_ENERGY_AU_TO_PER_SEC )
        elseif  JAC.give("unit: energy") == "A"             return( CONVERT_ENERGY_AU_TO_KAYSERS / (1.0e8 * wa) )
        else    error("stop d")
        end

    elseif   sa in ["energy: from eV to atomic"]            return( wa / CONVERT_ENERGY_AU_TO_EV )

    elseif    sa in ["rate: from atomic to predefined unit", "rate: from atomic"]
        if       JAC.give("unit: rate") == "1/s"            return( wa * CONVERT_RATE_AU_TO_PER_SEC )
        elseif   JAC.give("unit: rate") == "a.u."           return( wa )
        else     error("stop e")
        end
      
    elseif   sa in ["rate: from atomic to 1/s"]             return( wa * CONVERT_RATE_AU_TO_PER_SEC )

    elseif   sa in ["rate: from predefined to atomic unit", "rate: to atomic"]
        if      JAC.give("rate: time") == "1/s"             return( wa / CONVERT_RATE_AU_TO_PER_SEC )
        elseif  JAC.give("rate: time") == "a.u."            return( wa  )
        else    error("stop f")
        end

    elseif  sa in ["time: from atomic to predefined unit", "time: from atomic"]
        if      JAC.give("unit: time") == "sec"             return( wa * CONVERT_TIME_AU_TO_SEC )
        elseif  JAC.give("unit: time") == "fs"              return( wa * CONVERT_TIME_AU_TO_SEC * 10.0e15 )
        elseif  JAC.give("unit: time") == "as"              return( wa * CONVERT_TIME_AU_TO_SEC * 10.0e18 )
        elseif  JAC.give("unit: time") == "a.u."            return( wa )
        else    error("stop g")
        end
      
    elseif   sa in ["time: from atomic to sec"]             return( wa * CONVERT_TIME_AU_TO_SEC )
    elseif   sa in ["time: from atomic to fs"]              return( wa * CONVERT_TIME_AU_TO_SEC * 10.0e15 )
    elseif   sa in ["time: from atomic to as"]              return( wa * CONVERT_TIME_AU_TO_SEC * 10.0e18 )

    elseif  sa in ["time: from predefined to atomic unit", "time: to atomic"]
        if      JAC.give("unit: time") == "sec"             return( wa / CONVERT_TIME_AU_TO_SEC )
        elseif  JAC.give("unit: time") == "fs"              return( wa / (CONVERT_TIME_AU_TO_SEC * 10.0e15) )
        elseif  JAC.give("unit: time") == "as"              return( wa / (CONVERT_TIME_AU_TO_SEC * 10.0e18) )
        elseif  JAC.give("unit: time") == "a.u."            return( wa )
        else    error("stop h")
        end

    elseif  sa in ["length: from fm to atomic"]        return (wa / CONVERT_LENGTH_AU_TO_FEMTOMETER )
    elseif  sa in ["length: from atomic to fm"]        return (wa * CONVERT_LENGTH_AU_TO_FEMTOMETER )

    elseif  sa in ["kinetic energy to wave number: atomic units"]
        c = INVERSE_FINE_STRUCTURE_CONSTANT
        wb = sqrt( wa*wa/(c*c) + 2wa );                     return( wb )

    elseif  sa in ["kinetic energy to wavelength: atomic units"]
        c = INVERSE_FINE_STRUCTURE_CONSTANT
        wb = sqrt( wa*wa/(c*c) + 2wa );                     return( 2pi / wb )

    elseif  sa in ["wave number to total electron energy: atomic units"]
        c = INVERSE_FINE_STRUCTURE_CONSTANT
        wb = sqrt( wa*wa/(c*c) + 2wa );                     return( sqrt(wb*wb*c*c + c^4) )

    elseif  sa in ["wave number to kinetic energy: atomic units"]
        c = INVERSE_FINE_STRUCTURE_CONSTANT
        wb = sqrt( wa*wa/(c*c) + 2wa );                     return( c + sqrt(wb*wb*c*c + c^2) )

    else
        error("Unsupported keystring = $sa")
    end
end

