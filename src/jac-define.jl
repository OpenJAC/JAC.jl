

export  define

"""
`JAC.define()`  ... (re-) defines some 'standard' settings which are common to all the computations with the JAC module, and which can 
                    be 'overwritten' by the user. --- An improper setting of some variable may lead to an error message, if recognized
                    immediately. The following defaults apply if not specified otherwise by the user: the framework is 'relativistic', 
                    energies are given in eV and cross sections in barn. Note that, internally, atomic units are used throughout for 
                    all the computations within the program. nothing is returned if not indicated otherwise.

  + `("framework: relativistic")`  or  `("framework: non-relativistic")`   ... to define a relativistic or non-relativistic framework 
                   for all subsequent computations. 

  + `("method: continuum, spherical Bessel")`  or  `("method: continuum, pure sine")`  or  `("method: continuum, asymptotic Coulomb")`  
                   or  `("method: continuum, nonrelativistic Coulomb")`  or  `("method: continuum, Galerkin")`  
                   ... to define a a method for the generation of the continuum orbitals as (pure) spherical Bessel, pure sine,
                   asymptotic Coulomb, nonrelativistic Coulomb orbital or by means of the B-spline-Galerkin method.

  + `("method: normalization, pure sine")`  or  `("method: normalization, pure Coulomb")`   ... to define a method for the 
                   normalization of the continuum orbitals as asymptotically (pure) sine or Coulomb functions.

  + `("QED model: Petersburg")`  or  `("QED model: Sydney")`   ... to define a model for the computation of the QED corrections following
                   the work by Shabaev et al. (2011; Petersburg) or Flambaum and Ginges (2004; Syney).

  + `("unit: energy", "eV")`  or  `("unit: energy", "Kayser")`  or  `("unit: energy", "Hartree")`  or  `("unit: energy", "Hz")`  or  
    `("unit: energy", "Hz")`  ... to (pre-) define the energy units for all further printouts and communications with the JAC module.

  + `("unit: cross section", "a.u.")`  or  `("unit: cross section", "barn")`  or  `("unit: cross section", "Mbarn")`  ... to (pre-) 
                   define the unit for the printout of cross sections.

  + `("unit: rate", "a.u.")`  or  `("unit: rate", "1/s")`  ... to (pre-) define the unit for the printout of rates.

  + `("unit: time", "a.u.")`  or  `("unit: time", "sec")`  or  `("unit: time", "fs")`  or  `("unit: time", "as")`  ... to (pre-) define 
                   the unit for the printout and communications of times with the JAC module.
"""
function define(sa::String)
    global JAC_FRAMEWORK, JAC_CONT_SOLUTION, JAC_CONT_NORMALIZATION, JAC_QED_MODEL

    if        sa == "framework: relativistic"                            JAC_FRAMEWORK           = "relativistic"    
    elseif    sa == "framework: non-relativistic"                        JAC_FRAMEWORK           = "non-relativistic"  
    elseif    sa == "method: continuum, spherical Bessel"                JAC_CONT_SOLUTION       = ContBessel  
    elseif    sa == "method: continuum, pure sine"                       JAC_CONT_SOLUTION       = ContSine  
    elseif    sa == "method: continuum, asymptotic Coulomb"              JAC_CONT_SOLUTION       = AsymptoticCoulomb  
    elseif    sa == "method: continuum, nonrelativistic Coulomb"         JAC_CONT_SOLUTION       = NonrelativisticCoulomb  
    elseif    sa == "method: continuum, Galerkin"                        JAC_CONT_SOLUTION       = BsplineGalerkin   
    elseif    sa == "method: normalization, pure sine"                   JAC_CONT_NORMALIZATION  = PureSine   
    elseif    sa == "method: normalization, pure Coulomb"                JAC_CONT_NORMALIZATION  = CoulombSine  
    elseif    sa == "QED model: Petersburg"                              JAC_QED_MODEL           = QedPetersburg
    elseif    sa == "QED model: Sydney"                                  JAC_QED_MODEL           = QedSydney
    else      error("Unsupported keystring:: $sa")
    end
    nothing
end


function define(sa::String, sb::String)
    global JAC_ENERGY_UNIT, JAC_CROSS_SECTION_UNIT, JAC_RATE_UNIT, JAC_TIME_UNIT, JAC_SUMMARY_IOSTREAM, JAC_PRINT_SUMMARY, 
           JAC_TEST_IOSTREAM, JAC_PRINT_TEST 

    if        sa == "unit: energy"
        units = ["eV", "Kayser", "Hartree", "Hz", "A"]
        !(sb in units)    &&    error("Currently supported energy units: $(units)")
        JAC_ENERGY_UNIT = sb
    elseif    sa == "unit: cross section"
        units = ["a.u.", "barn", "Mbarn"]
        !(sb in units)    &&    error("Currently supported cross section units: $(units)")
        JAC_CROSS_SECTION_UNIT = sb
    elseif    sa == "unit: rate"
        units = ["a.u.", "1/s"]
        !(sb in units)    &&    error("Currently supported rate units: $(units)")
        JAC_RATE_UNIT = sb
    elseif    sa == "unit: time"
        units = ["a.u.", "sec", "fs", "as"]
        !(sb in units)    &&    error("Currently supported time units: $(units)")
        JAC_TIME_UNIT = sb
    elseif    sa == "print summary: open"
        JAC_PRINT_SUMMARY    = true
        JAC_SUMMARY_IOSTREAM = open(sb, "w") 
        println(JAC_SUMMARY_IOSTREAM, "Summary file opened at $( string(now())[1:16] ): \n" *
                                      "========================================   \n")
    elseif    sa == "print summary: append"
        JAC_PRINT_SUMMARY    = true
        JAC_SUMMARY_IOSTREAM = open(sb, "a") 
        println(JAC_SUMMARY_IOSTREAM, "Summary file re-opened (to append) at $( string(now())[1:16]): \n" *
                                      "======================================================= \n")
    elseif    sa == "print summary: close"
        JAC_PRINT_SUMMARY    = false
        close(JAC_SUMMARY_IOSTREAM)
    elseif    sa == "print test: open"
        JAC_PRINT_TEST    = true
        JAC_TEST_IOSTREAM = open(sb, "w") 
        println(JAC_TEST_IOSTREAM, "Test report file opened at $( string(now())[1:16] ): \n" *
                                   "============================================  \n")
    elseif    sa == "print test: close"
        JAC_PRINT_TEST    = false
        close(JAC_TEST_IOSTREAM)
    else      error("Unsupported keystring:: $sa")
    end

    nothing
end


"""
  + `("relativistic subshell list", subshells::Array{Subshell,1}; printout::Bool=true)`  
                    ... to (pre-) define internally the standard relativistic subshell list on which the standard order of orbitals is based.
"""
function define(sa::String, subshells::Array{Subshell,1}; printout::Bool=true)
    if        sa == "relativistic subshell list"
        if printout    println("(Re-) Define a new standard subshell list.")    end
        global JAC_STANDARD_SUBSHELL_LIST = deepcopy(subshells)
    else
        error("Unsupported keystring:: $sa")
    end

    nothing
end


"""
  + `("standard grid", grid::Radial.Grid)`  ... to (pre-) define internally the standard radial grid which is used to represent most orbitals.
"""
function define(sa::String, grid::Radial.Grid)
    global JAC_STANDARD_GRID

    if        sa == "standard grid"
        println("Re-define the standard grid with $(grid.NoPoints) grid points.")
        JAC_STANDARD_GRID = grid
    else
       error("Unsupported keystring:: $sa")
    end

    nothing
end


"""
  + `("QED: damped-hydrogenic", Znuc::Float64, wa::Array{Float64,1})`  ... to (re-) define the lambda-C damped overlap integrals of 
                                the lowest kappa-orbitals [ wa_1s_1/2, wa_2p_1/2, wa_2p_3/2, wa_3d_3/2, wa_3d_5/2 ] for the (new) 
                                nuclear charge Znuc; nothing is returned.
"""
function define(sa::String, Znuc::Float64, wa::Array{Float64,1})
    global JAC_QED_HYDROGENIC_LAMBDAC,  JAC_QED_NUCLEAR_CHARGE

    if        sa == "QED: damped-hydrogenic"
        if  length(wa) != 5  error("stop a ")   end
        println("Re-define the damped overlap integrals < a | e^{-r/lambda_C} | a > of the lowest kappa-orbitals for nuclear charge Z = " *
                "$(Znuc).")
        JAC_QED_NUCLEAR_CHARGE     = Znuc
        JAC_QED_HYDROGENIC_LAMBDAC = wa
    else
       error("Unsupported keystring:: $sa")
    end

    nothing
end

