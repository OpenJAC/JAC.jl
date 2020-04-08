
"""
`module  JAC.Defaults`  
    ... a submodel of JAC that contains all constants and global parameters of the program.
"""
module Defaults
    
    using Dates,  JAC, ..Basics,  ..Radial,  ..Math
    # 2014 CODATA recommended values, obtained from http://physics.nist.gov/cuu/Constants/
    
    export  convertUnits, getDefaults,  setDefaults

    # Dimensionless constants
    const FINE_STRUCTURE_CONSTANT         = 7.297_352_566_4e-3 
    const INVERSE_FINE_STRUCTURE_CONSTANT = 137.035_999_139

    # Constants in SI units
    const BOHR_RADIUS_SI                  = 0.529_177_210_67e-10
    const BOLTZMANN_CONSTANT_SI           = 1.380_648_52e-23
    const ELECTRON_MASS_SI                = 9.109_383_56e-31
    const ELEMENTARY_CHARGE_SI            = 1.602_176_620_8e-19
    const HARTREE_ENERGY_SI               = 4.359_744_650e-18
    const PLANCK_CONSTANT_SI              = 6.626_070_040e-34
    const PLANCK_CONSTANT_OVER_2_PI_SI    = 1.054_571_800e-34
    const SPEED_OF_LIGHT_IN_VACUUM_SI     = 299_792_458

    # Constants with eV instead of J
    const BOLTZMANN_CONSTANT_EV           = 8.617_330_3e-5
    const HARTREE_ENERGY_EV               = 27.211_386_02
    const PLANCK_CONSTANT_EV              = 4.135_667_662e-15
    const PLANCK_CONSTANT_OVER_2_PI_EV    = 6.582_119_514e-16

    # Constants with other atomic units
    const RYDBERG_IN_KAYSERS              = 109.73731534e3
    const ELECTRON_MASS_IN_G              = 9.1093897e-28
    const HBAR_IN_ERGS                    = 1.05457266e-27
    const ELECTRON_CHARGE_IN_ESU          = 4.80320680e-10
    
    # Constants in unified atomic mass units (u)
    const ELECTRON_MASS_U                 = 5.485_799_090_70e-4
    const NEUTRON_MASS_U                  = 1.008_664_915_88
    const PROTON_MASS_U                   = 1.007_276_466_879

    # Relationships between energy equivalents
    const ELECTRON_VOLT_ATOMIC_MASS_UNIT_RELATIONSHIP = 1.073_544_110_5e-9

    # Predefined conversion factors
    const CONVERT_ENERGY_AU_TO_EV           = HARTREE_ENERGY_EV
    const CONVERT_ENERGY_AU_TO_KAYSERS      = 2.0 * RYDBERG_IN_KAYSERS
    const CONVERT_ENERGY_AU_TO_PER_SEC      = 6.57968974479e15
    const CONVERT_TIME_AU_TO_SEC            = 2.418_884_254e-17
    const CONVERT_CROSS_SECTION_AU_TO_BARN  = BOHR_RADIUS_SI^2 * 1.0e28
    const CONVERT_RATE_AU_TO_PER_SEC        = (ELECTRON_MASS_IN_G/HBAR_IN_ERGS) * ((ELECTRON_CHARGE_IN_ESU^2/HBAR_IN_ERGS)^2)
    const CONVERT_LENGTH_AU_TO_FEMTOMETER   = BOHR_RADIUS_SI * 1.0e15

    # Predefined coefficients for numerical integration
    const FINITE_DIFFERENCE_NPOINTS         = 6
    const fivePointCoefficients             = Array{Float64}([2 * 7, 32, 12, 32, 7]) * 2 / 45  # * grid.h     
    const newtonCotesCoefficients           = fivePointCoefficients

    n = FINITE_DIFFERENCE_NPOINTS     
    ##x const weights = Array{Float64}(2*n + 1, 2*n + 1)     
    const weights = Array{Float64}(undef, 2*n + 1, 2*n + 1)
    for i = 1:2*n + 1       
        weights[i,:] = Math.finiteDifferenceWeights(-n + i - 1, 2*n + 1, order=1)[2,:]     
    end

    #
    # Global settings that can be (re-) defined by the user.
    GBL_FRAMEWORK               = "relativistic"
    GBL_CONT_SOLUTION           = AsymptoticCoulomb         ###  ContBessel    ContSine     AsymptoticCoulomb    NonrelativisticCoulomb    BsplineGalerkin
    GBL_CONT_NORMALIZATION      = PureSine                  ###  PureSine    CoulombSine    OngRussek
    GBL_QED_HYDROGENIC_LAMBDAC  = [1.0,  1.0,  1.0,  1.0,  1.0]
    GBL_QED_NUCLEAR_CHARGE      = 0.1
    ##x GBL_QED_MODEL               = QedPetersburg             ###  QedPetersburg  QedSydney    
    ##x GBL_CONT_PHASE              = WrtSin                ###  WrtSin
    GBL_WARNINGS                = String[]

    GBL_ENERGY_UNIT             = "eV"
    GBL_CROSS_SECTION_UNIT      = "barn"
    GBL_RATE_UNIT               = "1/s"
    GBL_TIME_UNIT               = "sec"

    GBL_PRINT_SUMMARY           = false
    GBL_PRINT_TEST              = false
    GBL_PRINT_DEBUG             = false
    
    GBL_SUMMARY_IOSTREAM        = nothing

    GBL_STANDARD_GRID           = Radial.Grid(true, printout=false)        
    ##x GBL_STANDARD_GRID       = Radial.Grid("grid: exponential", printout=false)        

    
    """
    `Defaults.convertUnits()`  
        ... converts some data from one format/unit into another one; cf. the supported keystrings and return values.

    + `("cross section: from atomic to predefined unit", value::Float64)`  or  `("cross section: from atomic", value::Float64)`  
        ... to convert an energy value from atomic to the predefined cross section unit; a Float64 is returned.

    + `("cross section: from atomic to barn", value::Float64)`  or  `("cross section: from atomic to Mbarn", value::Float64)`  or
      `("cross section: from atomic to Hz", value::Float64)`  or  `("energy: from atomic to Angstrom", value::Float64)` 
        ... to convert an energy value from atomic to the speficied cross section unit; a Float64 is returned.

    + `("cross section: from predefined to atomic unit", value::Float64)`  or  `("cross section: to atomic", value::Float64)`
        ... to convert a cross section value from the predefined to the atomic cross section unit; a Float64 is returned.

    + `("einstein B: from atomic", value::Float64)`  
        ... to convert a Einstein B coefficient from atomic to the speficied energy units; a Float64 is returned.or

    + `("energy: from atomic to eV", value::Float64)`  or  `("energy: from atomic to Kayser", value::Float64)`  or
      `("energy: from atomic to Hz", value::Float64)`  or  `("energy: from atomic to Angstrom", value::Float64)` 
        ... to convert an energy value from atomic to the speficied energy unit; a Float64 is returned.

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
    function convertUnits(sa::String, wa::Float64)
        global  CONVERT_ENERGY_AU_TO_EV, CONVERT_ENERGY_AU_TO_KAYSERS, CONVERT_ENERGY_AU_TO_PER_SEC, CONVERT_TIME_AU_TO_SEC,
                CONVERT_CROSS_SECTION_AU_TO_BARN,  CONVERT_RATE_AU_TO_PER_SEC,  CONVERT_LENGTH_AU_TO_FEMTOMETER 
    
        if       sa in ["cross section: from atomic to predefined unit", "cross section: from atomic"]
            if      Defaults.getDefaults("unit: cross section") == "a.u."   return( wa )
            elseif  Defaults.getDefaults("unit: cross section") == "barn"   return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN )
            elseif  Defaults.getDefaults("unit: cross section") == "Mbarn"  return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN * 10.0e-6 )
            else    error("stop a")
            end
        
        elseif   sa in ["cross section: from atomic to barn"]               return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN )
        elseif   sa in ["cross section: from atomic to Mbarn"]              return( wa * CONVERT_CROSS_SECTION_AU_TO_BARN * 10.0e-6 )

        elseif   sa in ["cross section: from predefined to atomic unit", "cross section: to atomic"]
            if      Defaults.getDefaults("unit: cross section") == "a.u."   return( wa )
            elseif  Defaults.getDefaults("unit: cross section") == "barn"   return( wa / CONVERT_CROSS_SECTION_AU_TO_BARN )
            elseif  Defaults.getDefaults("unit: cross section") == "Mbarn"  return( wa / CONVERT_CROSS_SECTION_AU_TO_BARN * 10.0e6 )
            else    error("stop b")
            end

        elseif   sa in ["Einstein B: from atomic"]                          return( wa )  # Einstein B not yet properly converted.

        elseif   sa in ["energy: from atomic to predefined unit", "energy: from atomic"]
            if      Defaults.getDefaults("unit: energy") == "eV"            return( wa * CONVERT_ENERGY_AU_TO_EV )
            elseif  Defaults.getDefaults("unit: energy") == "Kayser"        return( wa * CONVERT_ENERGY_AU_TO_KAYSERS )
            elseif  Defaults.getDefaults("unit: energy") == "Hartree"       return( wa )
            elseif  Defaults.getDefaults("unit: energy") == "Hz"            return( wa * CONVERT_ENERGY_AU_TO_PER_SEC )
            elseif  Defaults.getDefaults("unit: energy") == "A"             return( 1.0e8 / ( CONVERT_ENERGY_AU_TO_KAYSERS * wa ) )
            else    error("stop c")
            end
        
        elseif   sa in ["energy: from atomic to eV"]                        return( wa * CONVERT_ENERGY_AU_TO_EV )
        elseif   sa in ["energy: from atomic to Kayser"]                    return( wa * CONVERT_ENERGY_AU_TO_KAYSERS )
        elseif   sa in ["energy: from atomic to Hz"]                        return( wa * CONVERT_ENERGY_AU_TO_PER_SEC )
        elseif   sa in ["energy: from atomic to Angstrom"]                  return( 1.0e8 / ( CONVERT_ENERGY_AU_TO_KAYSERS * wa ) )

        elseif   sa in ["energy: from predefined to atomic unit", "energy: to atomic"]
            if      Defaults.getDefaults("unit: energy") == "eV"            return( wa / CONVERT_ENERGY_AU_TO_EV )
            elseif  Defaults.getDefaults("unit: energy") == "Kayser"        return( wa / CONVERT_ENERGY_AU_TO_KAYSERS )
            elseif  Defaults.getDefaults("unit: energy") == "Hartree"       return( wa )
            elseif  Defaults.getDefaults("unit: energy") == "Hz"            return( wa / CONVERT_ENERGY_AU_TO_PER_SEC )
            elseif  Defaults.getDefaults("unit: energy") == "A"             return( CONVERT_ENERGY_AU_TO_KAYSERS / (1.0e8 * wa) )
            else    error("stop d")
            end

        elseif   sa in ["energy: from eV to atomic"]                        return( wa / CONVERT_ENERGY_AU_TO_EV )

        elseif    sa in ["rate: from atomic to predefined unit", "rate: from atomic"]
            if       Defaults.getDefaults("unit: rate") == "1/s"            return( wa * CONVERT_RATE_AU_TO_PER_SEC )
            elseif   Defaults.getDefaults("unit: rate") == "a.u."           return( wa )
            else     error("stop e")
            end
        
        elseif   sa in ["rate: from atomic to 1/s"]                         return( wa * CONVERT_RATE_AU_TO_PER_SEC )

        elseif   sa in ["rate: from predefined to atomic unit", "rate: to atomic"]
            if      Defaults.getDefaults("rate: time") == "1/s"             return( wa / CONVERT_RATE_AU_TO_PER_SEC )
            elseif  Defaults.getDefaults("rate: time") == "a.u."            return( wa  )
            else    error("stop f")
            end

        elseif  sa in ["time: from atomic to predefined unit", "time: from atomic"]
            if      Defaults.getDefaults("unit: time") == "sec"             return( wa * CONVERT_TIME_AU_TO_SEC )
            elseif  Defaults.getDefaults("unit: time") == "fs"              return( wa * CONVERT_TIME_AU_TO_SEC * 10.0e15 )
            elseif  Defaults.getDefaults("unit: time") == "as"              return( wa * CONVERT_TIME_AU_TO_SEC * 10.0e18 )
            elseif  Defaults.getDefaults("unit: time") == "a.u."            return( wa )
            else    error("stop g")
            end
        
        elseif   sa in ["time: from atomic to sec"]                         return( wa * CONVERT_TIME_AU_TO_SEC )
        elseif   sa in ["time: from atomic to fs"]                          return( wa * CONVERT_TIME_AU_TO_SEC * 10.0e15 )
        elseif   sa in ["time: from atomic to as"]                          return( wa * CONVERT_TIME_AU_TO_SEC * 10.0e18 )

        elseif  sa in ["time: from predefined to atomic unit", "time: to atomic"]
            if      Defaults.getDefaults("unit: time") == "sec"             return( wa / CONVERT_TIME_AU_TO_SEC )
            elseif  Defaults.getDefaults("unit: time") == "fs"              return( wa / (CONVERT_TIME_AU_TO_SEC * 10.0e15) )
            elseif  Defaults.getDefaults("unit: time") == "as"              return( wa / (CONVERT_TIME_AU_TO_SEC * 10.0e18) )
            elseif  Defaults.getDefaults("unit: time") == "a.u."            return( wa )
            else    error("stop h")
            end

        elseif  sa in ["length: from fm to atomic"]                         return (wa / CONVERT_LENGTH_AU_TO_FEMTOMETER )
        elseif  sa in ["length: from atomic to fm"]                         return (wa * CONVERT_LENGTH_AU_TO_FEMTOMETER )

        elseif  sa in ["kinetic energy to wave number: atomic units"]
            c = INVERSE_FINE_STRUCTURE_CONSTANT
            wb = sqrt( wa*wa/(c*c) + 2wa );                                 return( wb )

        elseif  sa in ["kinetic energy to wavelength: atomic units"]
            c = INVERSE_FINE_STRUCTURE_CONSTANT
            wb = sqrt( wa*wa/(c*c) + 2wa );                                 return( 2pi / wb )

        elseif  sa in ["wave number to total electron energy: atomic units"]
            c = INVERSE_FINE_STRUCTURE_CONSTANT
            wb = sqrt( wa*wa/(c*c) + 2wa );                                 return( sqrt(wb*wb*c*c + c^4) )

        elseif  sa in ["wave number to kinetic energy: atomic units"]
            c = INVERSE_FINE_STRUCTURE_CONSTANT
            wb = sqrt( wa*wa/(c*c) + 2wa );                                 return( c + sqrt(wb*wb*c*c + c^2) )

        else
            error("Unsupported keystring = $sa")
        end
    end



    """
    `Defaults.setDefaults()`  
        ... (re-) defines some 'standard' settings which are common to all the computations with the JAC module, and which can 
            be 'overwritten' by the user. --- An improper setting of some variable may lead to an error message, if recognized
            immediately. The following defaults apply if not specified otherwise by the user: the framework is 'relativistic', 
            energies are given in eV and cross sections in barn. Note that, internally, atomic units are used throughout for 
            all the computations within the program. nothing is returned if not indicated otherwise.

    + `("framework: relativistic")`  or  `("framework: non-relativistic")`   
        ... to define a relativistic or non-relativistic framework for all subsequent computations. 

    + `("method: continuum, spherical Bessel")`  or  `("method: continuum, pure sine")`  or  
      `("method: continuum, asymptotic Coulomb")`  or  `("method: continuum, nonrelativistic Coulomb")`  or  
      `("method: continuum, Galerkin")`  
        ... to define a a method for the generation of the continuum orbitals as (pure) spherical Bessel, pure sine,
            asymptotic Coulomb, nonrelativistic Coulomb orbital or by means of the B-spline-Galerkin method.

    + `("method: normalization, pure sine")`  or  `("method: normalization, pure Coulomb")`   
        ... to define a method for the normalization of the continuum orbitals as asymptotically (pure) sine or Coulomb 
            functions.

    + `("QED model: Petersburg")`  or  `("QED model: Sydney")`   
        ... to define a model for the computation of the QED corrections following the work by Shabaev et al. (2011; Petersburg) 
            or Flambaum and Ginges (2004; Syney).

    + `("unit: energy", "eV")`  or  `("unit: energy", "Kayser")`  or  `("unit: energy", "Hartree")`  or  
      `("unit: energy", "Hz")`  or  `("unit: energy", "Hz")`  
        ... to (pre-) define the energy units for all further printouts and communications with the JAC module.

    + `("unit: cross section", "a.u.")`  or  `("unit: cross section", "barn")`  or  `("unit: cross section", "Mbarn")`  
        ... to (pre-) define the unit for the printout of cross sections.

    + `("unit: rate", "a.u.")`  or  `("unit: rate", "1/s")`  ... to (pre-) define the unit for the printout of rates.

    + `("unit: time", "a.u.")`  or  `("unit: time", "sec")`  or  `("unit: time", "fs")`  or  `("unit: time", "as")`  
        ... to (pre-) define the unit for the printout and communications of times with the JAC module.
    """
    function setDefaults(sa::String)
        global GBL_FRAMEWORK, GBL_CONT_SOLUTION, GBL_CONT_NORMALIZATION  ##x , GBL_QED_MODEL

        if        sa == "framework: relativistic"                            GBL_FRAMEWORK           = "relativistic"    
        elseif    sa == "framework: non-relativistic"                        GBL_FRAMEWORK           = "non-relativistic"  
        elseif    sa == "method: continuum, spherical Bessel"                GBL_CONT_SOLUTION       = ContBessel  
        elseif    sa == "method: continuum, pure sine"                       GBL_CONT_SOLUTION       = ContSine  
        elseif    sa == "method: continuum, asymptotic Coulomb"              GBL_CONT_SOLUTION       = AsymptoticCoulomb  
        elseif    sa == "method: continuum, nonrelativistic Coulomb"         GBL_CONT_SOLUTION       = NonrelativisticCoulomb  
        elseif    sa == "method: continuum, Galerkin"                        GBL_CONT_SOLUTION       = BsplineGalerkin   
        elseif    sa == "method: normalization, pure sine"                   GBL_CONT_NORMALIZATION  = PureSine   
        elseif    sa == "method: normalization, pure Coulomb"                GBL_CONT_NORMALIZATION  = CoulombSine  
        ##x elseif    sa == "QED model: Petersburg"                              GBL_QED_MODEL           = QedPetersburg
        ##x elseif    sa == "QED model: Sydney"                                  GBL_QED_MODEL           = QedSydney
        else      error("Unsupported keystring:: $sa")
        end
        nothing
    end


    function setDefaults(sa::String, sb::String)
        ## global GBL_ENERGY_UNIT, GBL_CROSS_SECTION_UNIT, GBL_RATE_UNIT, GBL_TIME_UNIT, GBL_SUMMARY_IOSTREAM, GBL_PRINT_SUMMARY, 
        ##    GBL_TEST_IOSTREAM, GBL_PRINT_TEST 

        if        sa == "unit: energy"
            units = ["eV", "Kayser", "Hartree", "Hz", "A"]
            !(sb in units)    &&    error("Currently supported energy units: $(units)")
            GBL_ENERGY_UNIT = sb
        elseif    sa == "unit: cross section"
            units = ["a.u.", "barn", "Mbarn"]
            !(sb in units)    &&    error("Currently supported cross section units: $(units)")
            GBL_CROSS_SECTION_UNIT = sb
        elseif    sa == "unit: rate"
            units = ["a.u.", "1/s"]
            !(sb in units)    &&    error("Currently supported rate units: $(units)")
            GBL_RATE_UNIT = sb
        elseif    sa == "unit: time"
            units = ["a.u.", "sec", "fs", "as"]
            !(sb in units)    &&    error("Currently supported time units: $(units)")
            GBL_TIME_UNIT = sb
        elseif    sa == "print summary: open"
            global GBL_PRINT_SUMMARY    = true
            global GBL_SUMMARY_IOSTREAM = open(sb, "w") 
            println(GBL_SUMMARY_IOSTREAM, "Summary file opened at $( string(now())[1:16] ): \n" *
                                        "========================================   \n")
        elseif    sa == "print summary: append"
            GBL_PRINT_SUMMARY    = true
            GBL_SUMMARY_IOSTREAM = open(sb, "a") 
            println(GBL_SUMMARY_IOSTREAM, "Summary file re-opened (to append) at $( string(now())[1:16]): \n" *
                                        "======================================================= \n")
        elseif    sa == "print summary: close"
            GBL_PRINT_SUMMARY    = false
            close(GBL_SUMMARY_IOSTREAM)
        elseif    sa == "print test: open"
            GBL_PRINT_TEST    = true
            JAC.JAC_TEST_IOSTREAM = open(sb, "w") 
            println(JAC.JAC_TEST_IOSTREAM, "Test report file opened at $( string(now())[1:16] ): \n" *
                                    "============================================  \n")
        elseif    sa == "print test: close"
            GBL_PRINT_TEST    = false
            close(GBL_TEST_IOSTREAM)
        else      error("Unsupported keystring:: $sa")
        end

        nothing
    end


    """
    + `("relativistic subshell list", subshells::Array{Subshell,1}; printout::Bool=true)`  
        ... to (pre-) define internally the standard relativistic subshell list on which the standard order of orbitals is based.
    """
    function setDefaults(sa::String, subshells::Array{Subshell,1}; printout::Bool=true)
        if        sa == "relativistic subshell list"
            if printout    println("(Re-) Define a new standard subshell list.")    end
            global GBL_STANDARD_SUBSHELL_LIST = deepcopy(subshells)
        else
            error("Unsupported keystring:: $sa")
        end

        nothing
    end


    """
    + `("standard grid", grid::Radial.Grid; printout::Bool=true)`  
        ... to (pre-) define internally the standard radial grid which is used to represent most orbitals.
    """
    function setDefaults(sa::String, grid::Radial.Grid; printout::Bool=true)
        global GBL_STANDARD_GRID

        if        sa == "standard grid"
            if  printout    println("(Re-) Define the standard grid with $(grid.NoPoints) grid points.")    end
            GBL_STANDARD_GRID = grid
        else
        error("Unsupported keystring:: $sa")
        end

        nothing
    end


    """
    + `("QED: damped-hydrogenic", Znuc::Float64, wa::Array{Float64,1})`  
        ... to (re-) define the lambda-C damped overlap integrals of the lowest kappa-orbitals 
            [ wa_1s_1/2, wa_2p_1/2, wa_2p_3/2, wa_3d_3/2, wa_3d_5/2 ] for the (new) nuclear charge Znuc; 
            nothing is returned.
    """
    function setDefaults(sa::String, Znuc::Float64, wa::Array{Float64,1})
        global GBL_QED_HYDROGENIC_LAMBDAC,  GBL_QED_NUCLEAR_CHARGE

        if        sa == "QED: damped-hydrogenic"
            if  length(wa) != 5  error("stop a ")   end
            println("Re-define the damped overlap integrals < a | e^{-r/lambda_C} | a > of the lowest kappa-orbitals for nuclear charge Z = " *
                    "$(Znuc).")
            GBL_QED_NUCLEAR_CHARGE     = Znuc
            GBL_QED_HYDROGENIC_LAMBDAC = wa
        else
        error("Unsupported keystring:: $sa")
        end

        nothing
    end

        
        
    """
    `Defaults.getDefaults()`  
        ... gives/supplies different information about the (present) framework of the computation or about some 
            given data; cf. Defaults.setDefaults(). 

    + `("alpha")`  or  `("fine-structure constant alpha")`   
        ... to get the (current) value::Float64 of the fine-structure constant alpha.

    + `("electron mass: kg")`  or  `("electron mass: amu")`  
        ... to get the (current) value::Float64 of the electron mass in the specified unit.

    + `("framework")`  ... to give the (current) setting::String  of the overall framework.

    + `("electron rest energy")`  or  `("mc^2")`  ... to get the electron rest energy.

    + `("electron g-factor")`  ... to give the electron g-factor g_s = 2.00232.

    + `("unit: energy")`  or  `("unit: cross section")`  or  `("unit: rate")`  
        ... to get the corresponding (user-defined) unit::String for the current computations.

    + `("standard grid")`  
        ... to get the (current standard) grid::Array{Float64,1} to which all radial orbital functions usually refer.

    + `("speed of light: c")`  ... to get the speed of light in atomic units.

    + `("summary flag/stream")`  
        ... to get the logical flag and stream for printing a summary file; a tupel (flag, iostream) is returned.
    """
    function getDefaults(sa::String)
        global GBL_FRAMEWORK, GBL_CONTINUUM, GBL_ENERGY_UNIT, GBL_CROSS_SECTION_UNIT, GBL_RATE_UNIT
        global ELECTRON_MASS_SI, ELECTRON_MASS_U, FINE_STRUCTURE_CONSTANT, GBL_STANDARD_GRID, GBL_SUMMARY_IOSTREAM, GBL_PRINT_SUMMARY, 
               GBL_TEST_IOSTREAM, GBL_TEST_SUMMARY

        if        sa in ["alpha", "fine-structure constant alpha"]          return (FINE_STRUCTURE_CONSTANT)   
        elseif    sa == "electron mass: kg"                                 return (ELECTRON_MASS_SI)   
        elseif    sa == "electron mass: amu"                                return (ELECTRON_MASS_U)   
        elseif    sa == "electron g-factor"                                 return (2.00232)   
        elseif    sa == "framework"                                         return (GBL_FRAMEWORK)
        elseif    sa == "method: continuum"                                 return (GBL_CONTINUUM)
        elseif    sa in ["electron rest energy", "mc^2"]                    return (1/(FINE_STRUCTURE_CONSTANT^2))   
        elseif    sa == "unit: energy"                                      return (GBL_ENERGY_UNIT)
        elseif    sa == "unit: cross section"                               return (GBL_CROSS_SECTION_UNIT)
        elseif    sa == "unit: rate"                                        return (GBL_RATE_UNIT)
        elseif    sa == "unit: time"                                        return (GBL_TIME_UNIT)
        elseif    sa == "standard grid"                                     return (GBL_STANDARD_GRID)
        elseif    sa == "speed of light: c"                                 return (1.0/FINE_STRUCTURE_CONSTANT)
        elseif    sa == "summary flag/stream"                               return ( (GBL_PRINT_SUMMARY, GBL_SUMMARY_IOSTREAM) )
        elseif    sa == "test flag/stream"                                  return ( (GBL_PRINT_TEST,    JAC.JAC_TEST_IOSTREAM) )
        else      error("Unsupported keystring:: $sa")
        end
    end


    """
    + `("ordered shell list: non-relativistic", n_max::Int64)`  
        ... to give an ordered list of non-relativistic shells::Array{Shell,1} up to the (maximum) principal number n_max.
            
    + `("ordered subshell list: relativistic", n_max::Int64)`   
        ... to give an ordered list of relativistic subshells::Array{Subshell,1} up to the (maximum) principal number n_max.
    """
    function getDefaults(sa::String, n_max::Int64)
    
        !(1 <= n_max < 10)    &&    error("Unsupported value of n_max = $n_max")

        if        sa == "ordered shell list: non-relativistic"
            wa = Shell[]
            for  n = 1:n_max
            for l = 0:n_max-1    push!(wa, Shell(n,l) )    end
            end
        elseif    sa == "ordered subshell list: relativistic"  
            wa = Subshell[]
            for  n = 1:n_max
                for l = 0:n_max-1
                    if     l == 0   push!(wa, Subshell(n, -1) )
                    else   j2 = 2l - 1;    kappa = Int64(  (j2+1)/2 );   push!(wa, Subshell(n, kappa) )
                        j2 = 2l + 1;    kappa = Int64( -(j2+1)/2 );   push!(wa, Subshell(n, kappa) )
                    end
                end
            end
        else    error("Unsupported keystring:: $sa")
        end
        return( wa )
    end


        
        
    """
    `Defaults.warn()`  
        ... deals with warnings that occur during a run and session; it handles the global array GBL_WARNINGS.

    + `(AddWarning, warning::String)`  ... to add warning to the global array GBL_WARNINGS.

    + `(PrintWarnings)`  ... to print all warnings that are currently kept in the global array GBL_WARNINGS.

    + `(ResetWarnings)`  ... to reset the global array GBL_WARNINGS.
    """
    function warn(wa::Warnings)
        global GBL_WARNINGS

        if        wa == PrintWarnings
            iostream = open("jac-warn.report", "w") 
            println(iostream, " ")
            println(iostream, "\n\n",
                            "Warnings from the present sessions or run, printed at $( string(now())[1:16] ): \n",
                            "======================================================================= \n")
            for sa in  GBL_WARNINGS
                println(iostream, "++ " * sa)
            end
            close(iostream)
            #
        elseif    wa == ResetWarnings
            ## @warn("Reset global array GBL_WARNINGS.")
            printstyled("Constants.warn():  Reset global array GBL_WARNINGS.", color=:light_magenta)
            GBL_WARNINGS = String[]
        else      error("Unsupported Warnings:: $wa")
        end

        nothing
    end


    function warn(wa::Warnings, sa::String)
        global GBL_WARNINGS

        if        wa == AddWarning
            push!(GBL_WARNINGS, sa)
        else      error("Unsupported Warnings:: $wa")
        end

        nothing
    end


end # module

