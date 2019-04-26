
export  give

"""
`JAC.give()`  ... gives/supplies different information about the (present) framework of the computation or about some given data; 
                  cf. JAC.define(). 

  + `("alpha")`  or  `("fine-structure constant alpha")`   ... to give the (current) value::Float64 of the fine-structure constant alpha.

  + `("electron mass: kg")`  or  `("electron mass: amu")`  ... to give the (current) value::Float64 of the electron mass in the specified unit.

  + `("framework")`  ... to give the (current) setting::String  of the overall framework.

  + `("electron rest energy")`  or  `("mc^2")`  ... to give the electron rest energy.

  + `("electron g-factor")`  ... to give the electron g-factor g_s = 2.00232.

  + `("unit: energy")`  or  `("unit: cross section")`  or  `("unit: rate")`  ... to give the corresponding (user-defined) unit::String for 
                  the current computations.

  + `("standard grid")`  ... to give the (current standard) grid::Array{Float64,1} to which all radial orbital functions usually refer.

  + `("speed of light: c")`  ... to give the speed of light in atomic units.

  + `("summary flag/stream")`  ... to give the logical flag and stream for printing a summary file; a tupel (flag, iostream) is returned.
"""
function give(sa::String)
    global JAC_FRAMEWORK, JAC_CONTINUUM, JAC_ENERGY_UNIT, JAC_CROSS_SECTION_UNIT, JAC_RATE_UNIT
    global ELECTRON_MASS_SI, ELECTRON_MASS_U, FINE_STRUCTURE_CONSTANT, JAC_STANDARD_GRID, JAC_SUMMARY_IOSTREAM, JAC_PRINT_SUMMARY, 
           JAC_TEST_IOSTREAM, JAC_TEST_SUMMARY

    if        sa in ["alpha", "fine-structure constant alpha"]          return (JAC.FINE_STRUCTURE_CONSTANT)   
    elseif    sa == "electron mass: kg"                                 return (JAC.ELECTRON_MASS_SI)   
    elseif    sa == "electron mass: amu"                                return (JAC.ELECTRON_MASS_U)   
    elseif    sa == "electron g-factor"                                 return (2.00232)   
    elseif    sa == "framework"                                         return (JAC.JAC_FRAMEWORK)
    elseif    sa == "method: continuum"                                 return (JAC.JAC_CONTINUUM)
    elseif    sa in ["electron rest energy", "mc^2"]                    return (1/(JAC.FINE_STRUCTURE_CONSTANT^2))   
    elseif    sa == "unit: energy"                                      return (JAC.JAC_ENERGY_UNIT)
    elseif    sa == "unit: cross section"                               return (JAC.JAC_CROSS_SECTION_UNIT)
    elseif    sa == "unit: rate"                                        return (JAC.JAC_RATE_UNIT)
    elseif    sa == "unit: time"                                        return (JAC.JAC_TIME_UNIT)
    elseif    sa == "standard grid"                                     return (JAC.JAC_STANDARD_GRID)
    elseif    sa == "speed of light: c"                                 return (1.0/JAC.FINE_STRUCTURE_CONSTANT)
    elseif    sa == "summary flag/stream"                               return ( (JAC.JAC_PRINT_SUMMARY, JAC.JAC_SUMMARY_IOSTREAM) )
    elseif    sa == "test flag/stream"                                  return ( (JAC.JAC_PRINT_TEST,    JAC.JAC_TEST_IOSTREAM) )
    else      error("Unsupported keystring:: $sa")
    end
end


"""
  + `("ordered shell list: non-relativistic", n_max::Int64)`  ... to give an ordered list of non-relativistic shells::Array{Shell,1} up to 
                                                                  the (maximum) principal number n_max.
  + `("ordered subshell list: relativistic", n_max::Int64)`   ... to give an ordered list of relativistic subshells::Array{Subshell,1} up to the 
                                                                  (maximum) principal number n_max.
"""
function give(sa::String, n_max::Int64)
   
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

