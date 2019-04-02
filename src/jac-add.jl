

"""
`JAC.add()`  ... 'adds' vaious quantities in an obvious and appropriate manner.

  + `(pota::Radial.Potential, potb::Radial.Potential)`  ... to add two radial potentials together that are defined on the same grid.
                                                             A potential::RadialPotential is returned that inherits its radial size from the
                                                             potential that is defined in a larger range of r-values.
"""
function  add(pota::Radial.Potential, potb::Radial.Potential)
    if  pota.grid.NoPoints != potb.grid.NoPoints  ||   pota.grid.rnt != potb.grid.rnt  ||  pota.grid.h != potb.grid.h  ||  
        pota.grid.hp != potb.grid.hp                   error("stop a")   end
    name = pota.name * "+" * potb.name;   nmax = max(length(pota.Zr), length(potb.Zr));   nmin = min(length(pota.Zr), length(potb.Zr))
    Zr   = zeros(nmax);   
    nx = length(pota.Zr);    Zr[1:nx] = Zr[1:nx] + pota.Zr[1:nx] 
    nx = length(potb.Zr);    Zr[1:nx] = Zr[1:nx] + potb.Zr[1:nx]
    Zr[nmin+1:nmax] .= Zr[nmin] 

    potential = Radial.Potential(name, Zr, pota.grid)
    return( potential )
end



"""
`JAC.addZerosToCsfR(nz::Int64, csf::CsfR)`  ... 'adds' a number of zeros to the occupation, seniority, etc. of a given CsfR.
"""
function  addZerosToCsfR(nz::Int64, csf::CsfR)
    !csf.useStandardSubshells   &&   error("Zeros can only be added to CSF with standard subshell order.")
    
    println("nz = $nz")
    occupation = csf.occupation;   seniority = csf.seniority;    subshellJ = csf.subshellJ;    subshellX = csf.subshellX 
    
    for  i = 1:nz
        push!(occupation, 0);    push!(seniority, 0);    push!(subshellJ, AngularJ64(0));    push!(subshellX, subshellJ[end])
    end
    
    newCsf = CsfR( csf.useStandardSubshells, csf.J, csf.parity, occupation, seniority, subshellJ, subshellX, Subshell[])
    return(newCsf)
end

