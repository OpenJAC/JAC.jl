
export  integrate

"""
`JAC.integrate()`  ... integrates a one- or two-dimensional function by different numerical methods, either on a given grid or for a general
                       function.

  + `("function: on radial grid, Newton-Cotes", F::Array{Float64,1}, grid::Radial.Grid)`  ... to integrate the (radial) function F over the given
                                 (radial) grid, int_0^infinity dr F(r) by using a 5-point Newton-Cotes integration formula; a value::Float64 
                                 is returned.   

  + `("function: on radial grid, Simpson rule", F::Array{Float64,1}, grid::Radial.Grid)`  ... to integrate the (radial) function F over the given
                                 (radial) grid, int_0^infinity dr F(r) by using Simpson's rule; a value::Float64 is returned. 

  + `("function: on radial grid, trapez rule", F::Array{Float64,1}, grid::Radial.Grid)`  ... to integrate the (radial) function F over the given
                                 (radial) grid, int_0^infinity dr F(r), by using a simple trapez rule; a value::Float64 is returned. 
"""
function integrate(sa::String, F::Array{Float64,1}, grid::Radial.Grid)
    if       sa == "function: on radial grid, Newton-Cotes"     wa = integrateOnGridNewtonCotes(F, grid)   
    elseif   sa == "function: on radial grid, Simpson rule"     wa = integrateOnGridSimpsonRule(F, grid)
    elseif   sa == "function: on radial grid, trapez rule"      wa = integrateOnGridTrapezRule( F, grid)   
    else     error("Unsupported keystring = $sa.")
    end
    
end


"""
`JAC.integrateOnGridNewtonCotes()`  ... integrates by using a 5-point Newton-Cotes integration formula; see JAC.integrate().
"""
function integrateOnGridNewtonCotes(F::Array{Float64,1}, grid::Radial.Grid)
    coefficients = Array{Float64}([2 * 7, 32, 12, 32, 7]) * 2 / 45 * grid.h
    n = size(coefficients, 1)
  
    result = 0.;    i0     = 1
  
    while  abs(F[i0]) == 0
        i0 += 1
    end
    gamma  = log(F[i0]/F[i0+1] * grid.rp[i0+1]/grid.rp[i0]) / log(grid.r[i0]/grid.r[i0+1])
    result = 1/(gamma + 1) * grid.r[i0] * F[i0] / grid.rp[i0]
  
    result += F[i0] * coefficients[n]
    for i = i0+1:size(F, 1)
        result += coefficients[(i - i0) % (n - 1) + 1] * F[i]
    end
  
    if (size(F, 1) - i0) % (n - 1) + 1 == 1
        result -= coefficients[n] * F[size(F, 1)]
    end
  
    return( result )
end


"""
`JAC.integrateOnGridSimpsonRule()`  ... integrates by using Simpson's rule; see JAC.integrate().
"""
function integrateOnGridSimpsonRule(F::Array{Float64,1}, grid::Radial.Grid)
    coefficients = Array{Float64}([2 * 1, 4, 1]) / 3 * grid.h
    n = size(coefficients, 1)
  
    result = 0.
    result += F[1] * coefficients[n]
    for i = 2:size(F, 1)
        result += coefficients[(i - 1) % (n - 1) + 1] * F[i]
    end
  
    if (size(F, 1) - 1) % (n - 1) + 1 == 1
        result -= coefficients[n] * F[size(F, 1)]
    end
  
    return( result )
end


"""
`JAC.integrateOnGridTrapezRule()`  ... integrates by using a simple trapez rule; see JAC.integrate().
"""
function integrateOnGridTrapezRule(F::Array{Float64,1}, grid::Radial.Grid)
    coefficients = Array{Float64}([2 * 1, 1]) / 2 * grid.h
    n = size(coefficients, 1)
  
    result = 0.
  
    result += F[1] * coefficients[n]
    for i = 2:size(F, 1)
        result += coefficients[(i - 1) % (n - 1) + 1] * F[i]
    end
  
    if (size(F, 1) - 1) % (n - 1) + 1 == 1
        result -= coefficients[n] * F[size(F, 1)]
    end
  
    return( result )
end

