
"""
`Basics.interpolate()`  ... interpolates a one- or two-dimensional function by different numerical methods, either by value or on a whole grid.

  + `("function: to new radial grid, Grasp92", (F::Array{Float64,1}, oldGrid::Radial.Grid), (G::Void, newGrid::Radial.Grid) )`  
    ... to interpolate the (radial) function F from the oldgrid to newgrid, by a call to the Grasp92 Fortran procedures; 
        an function G::Array{Float64,1} is returned. **Not yet implemented !**    

  + `("function: to new radial grid, trapez rule", (F::Array{Float64,1}, oldGrid::Radial.Grid), (G::Void, newGrid::Radial.Grid)` 
    ... to integrate the (radial) function F from the oldgrid to newgrid, by using a trapez rule from the Numerical Recipies. 
        A function G::Array{Float64,1} is returned. **Not yet implemented !**
"""
function Basics.interpolate(sa::String, from::Tuple{Array{Float64,1},JAC.Radial.Grid}, to::Tuple{DataType,JAC.Radial.Grid})
    if       sa == "function: to new radial grid, Grasp92"          wa = interpolateOnGridGrasp92(from, to)   
    elseif   sa == "function: to new radial grid, trapez rule"      wa = interpolateOnGridTrapezRule(from, to)   
    else     error("Unsupported keystring = $sa.")
    end
    
end


"""
`Basics.interpolateOnGridGrasp92()`  ... interpolates ... by using Grasp92 Fortran procedures; see Basic.interpolate().
"""
function Basics.interpolateOnGridGrasp92(from::Tuple{Array{Float64,1},JAC.Radial.Grid}, to::Tuple{DataType,JAC.Radial.Grid})
    # First prepare all fields and settings to call quad_grasp92()
    error("Not yet implemented !")
end


"""
`Basics.interpolateOnGridTrapezRule()`  ... interpolates ... by using trapez rule from Numerical Recipies; see Basic.interpolate().
"""
function Basics.interpolateOnGridTrapezRule(from::Tuple{Array{Float64,1},JAC.Radial.Grid}, to::Tuple{DataType,JAC.Radial.Grid})
    error("Not yet implemented !")
end

