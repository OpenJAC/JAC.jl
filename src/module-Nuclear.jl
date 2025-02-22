
"""
`module JAC.Nuclear`  
... a submodel of JAC that contains procedures for defining a Nuclear.Model and for calculating various nuclear
    potentials.
"""
module Nuclear

using  QuadGK, ..Basics,  ..Defaults, ..Radial, ..Math
export Model


"""
`struct  Nuclear.Model`  ... defines a type for the nuclear model, i.e. for its form and parameters.

    + Z        ::Float64         ... nuclear charge
    + model    ::String          ... identifier of the nuclear model: {"Fermi", "Point", "Uniform"}
    + mass     ::Float64         ... atomic mass
    + radius   ::Float64         ... (root-mean square) radius of a uniform or Fermi-distributed nucleus
    + spinI    ::AngularJ64      ... nuclear spin I, must be >= 0
    + mu       ::Float64         ... magnetic dipole moment in Bohr magnetons
    + Q        ::Float64         ... electric quadrupole moment
"""
struct  Model
    Z          ::Float64
    model      ::String
    mass       ::Float64
    radius     ::Float64
    spinI      ::AngularJ64      
    mu         ::Float64          
    Q          ::Float64
end


"""
`Nuclear.Model(Z::Real)`  
    ... to specify a Fermi-type nucleus with charge Z, and where the nuclear spin and nuclear moments are all set to zero.
"""
function Model(Z::Real)
    Z < 0.1  &&  error("Z must be >= 0.1")
    model    = "Fermi"
    mass     = 2*Z + 0.005*Z^2
    radius   = Nuclear.Rrms(mass)
    spinI    = AngularJ64(0)
    mu       = 0.
    Q        = 0.

    Model(Z, model, mass, radius, spinI, mu, Q) 
end


"""
`Nuclear.Model(Z::Real, M::Float64)`  
    ... to specify a Fermi-type nucleus with charge Z and mass M, and where the nuclear spin and nuclear moments are all set to zero.
"""
function Model(Z::Real, M::Float64)
    Z < 0.1  &&  error("Z must be >= 0.1")
    model    = "Fermi"
    mass     = M
    radius   = Nuclear.Rrms(mass)
    spinI    = AngularJ64(0)
    mu       = 0.
    Q        = 0.

    Model(Z, model, mass, radius, spinI, mu, Q) 
end


"""
`Nuclear.Model(Z::Real, model::String)`  
    ... to specify a nucleus with charge Z, model = {"Fermi", "point", "uniform"}, and where the nuclear spin and 
        nuclear moments are all set to zero.
"""
function Model(Z::Real, model::String)
    Z < 0.1                                    &&  error("Z must be >= 0.1")
    !(model in ["Fermi", "point", "uniform"])  && error("Inappropriate model = $model")
    if  model == "point"
        mass     = 1000*Z
        radius   = 0.
    else
        mass     = 2*Z + 0.005*Z^2
        radius   = Nuclear.Rrms(mass)
    end
    spinI    = AngularJ64(0)
    mu       = 0.
    Q        = 0.

    Model(Z, model, mass, radius, spinI, mu, Q) 
end


"""
`Nuclear.Model(nm::Nuclear.Model;`
    
            Z=..,         model=..,         mass=..,        radius=..,     
            spinI=..,     mu=..,            Q=..)
    ... constructor for re-defining a nuclear model nm::Nuclear.Model.
"""
function Model(nm::Nuclear.Model;            Z::Union{Nothing,Float64}=nothing,          model::Union{Nothing,String}=nothing,         
    mass::Union{Nothing,Float64}=nothing,    radius::Union{Nothing,Float64}=nothing,     spinI::Union{Nothing,AngularJ64}=nothing,  
    mu::Union{Nothing,Float64}=nothing,      Q::Union{Nothing,Float64}=nothing)

    if  Z         == nothing   Zx          = nm.Z           else   Zx          = Z          end 
    if  model     == nothing   modelx      = nm.model       else   modelx      = model      end 
    if  mass      == nothing   massx       = nm.mass        else   massx       = mass       end 
    if  radius    == nothing   radiusx     = nm.radius      else   radiusx     = radius     end 
    if  spinI     == nothing   spinIx      = nm.spinI       else   spinIx      = spinI      end 
    if  mu        == nothing   mux         = nm.mu          else   mux         = mu         end 
    if  Q         == nothing   Qx          = nm.Q           else   Qx          = Q          end 
    
    Model(Zx, modelx, massx, radiusx, spinIx, mux, Qx)
end


# `Base.show(io::IO, m::Model)`  ... prepares a proper printout of the variable  m::Model.
function Base.show(io::IO, m::Model) 
    if      m.model == "Fermi"   
        print(io, "Fermi nuclear model for Z = $(m.Z) with mass = $(m.mass), radius R = $(m.radius) fm and ")
    elseif  m.model == "point"   
        print(io, "Point nuclear model for Z = $(m.Z) with mass = $(m.mass), radius R = $(m.radius) fm and ")
    elseif  m.model == "uniform"   
        print(io, "Uniform nuclear model for Z = $(m.Z) with mass = $(m.mass), radius R = $(m.radius) fm and ")
    else
        error("stop a")
    end

    print(io, "nuclear spin I = $(m.spinI), dipole moment mu = $(m.mu) and quadrupole moment Q = $(m.Q).")
end


"""
`struct  Nuclear.Compound`  
    ... defines a type for modeling isomeric (compound) nuclei with two nuclear symmetries, one common nuclear model
        and additional (nuclear) information. This type is of interest for studying hyperfine-induced transitions
        which involve different nuclear states.

    + Z             ::Float64      ... nuclear charge
    + model         ::String       ... identifier of the nuclear model: {"Fermi", "Point", "Uniform"}
    + mass          ::Float64      ... atomic mass
    + radius        ::Float64      ... (root-mean square) radius of a uniform or Fermi-distributed nucleus
    + lowerI        ::AngularJ64   ... nuclear spin I >= 0 of the lower (ground) nuclear level
    + lowerParity   ::Parity       ... parity of the lower (ground) nuclear level
    + upperI        ::AngularJ64   ... nuclear spin I >= 0 of the upper (isomeric) nuclear level
    + upperParity   ::Parity       ... parity of the upper (isomeric) nuclear level
    + upperEnergy   ::Float64      ... energy of the upper (isomeric) relative to the ground level [in user-specified units]
    + WE3element    ::Float64      ... <lowerI || W^(E3) || upperI> matrix element in [a.u.]
"""
struct  Compound
    Z               ::Float64
    model           ::String
    mass            ::Float64
    radius          ::Float64
    lowerI          ::AngularJ64 
    lowerParity     ::Parity
    upperI          ::AngularJ64
    upperParity     ::Parity
    upperEnergy     ::Float64 
    WE3element      ::Float64
end


"""
`Nuclear.Compound()`  ... constructor for an `empty` instance of Nuclear.Compound.
"""
function Compound()
    nm = Nuclear.Model(1.0)
    Compound(nm.Z, nm.model, nm.mass, nm.radius, nm.spinI, Basics.plus, nm.spinI, Basics.plus, 0., 0.)
end


# `Base.show(io::IO, compound::Nuclear.Compound)`  ... prepares a proper printout of the variable compound::Nuclear.Compound.
function Base.show(io::IO, compound::Nuclear.Compound) 
    println(io, "Z:              $(compound.Z)  ")
    println(io, "model:          $(compound.model)  ")
    println(io, "mass:           $(compound.mass)  ")
    println(io, "radius:         $(compound.radius)  ")
    println(io, "lowerI:         $(compound.lowerI)  ")
    println(io, "lowerParity:    $(compound.lowerParity)  ")
    println(io, "upperI:         $(compound.upperI)  ")
    println(io, "upperParity:    $(compound.upperParity)  ")
    println(io, "upperEnergy:    $(compound.upperEnergy)  ")
    println(io, "WE3element:     $(compound.WE3element)  ")
end

        
#################################################################################################################################
#################################################################################################################################



"""
`Nuclear.fermiA`  ... provides a value::Float64 for the fermi_a parameter.
"""
fermiA   = 2.3/(4 * log(3))


"""
`Nuclear.Rrms(A)`  
    ... provides a value::Float64 for the root-mean-squared radius R of a uniformly-distributed charge
        density of a nuclues with mass A.
"""
function Rrms(A)
    return( 0.836 * A^(1/3) + 0.57 )
end


"""
`Nuclear.fermiDensity(r, b)`  
    ... provides a value::Float64 for a fermi-distributed charge density at r and for given fermi_a
        and fermi_b parameters.
"""
function fermiDensity(r, b)
    return( 1/(1 + exp((r-b)/fermiA)) )
end


"""
`Nuclear.fermiRrms(b::Float64)`  
    ... provides a value::Float64 for the root-mean-squared radius of a fermi-distributed charge 
        density for given fermi_a and fermi_b parameters.
"""
function fermiRrms(b::Float64)

    return( sqrt(12 * fermiA^2 * Math.polylogExp(b/fermiA, 5) / Math.polylogExp(b/fermiA, 3)) )
end


"""
`Nuclear.computeFermiBParameter(R::Float64)`  
    ... computes a value::Float64 for the fermi_b parameter for a fermi-distributed nuclear with 
        root-mean square radius R.
"""
function computeFermiBParameter(R::Float64)

    function eq(b :: Float64)
        return fermiRrms(b) - R
    end

    b = R
    diff = eq(b)

    eps = 1E-14
    count = 0
    max = 100

    while  abs(diff) > eps   &&   count < max
        b = b -  diff
        
        diff = eq(b)
        count = count + 1
    end
    # println("Found b = $b, Δ = $(eq(b)) after $count iterations")
    
    return( b )
end


"""
`Nuclear.fermiDistributedNucleus(Rrms::Float64, Z::Float64, grid::Radial.Grid)`  
    ... computes the effective, radial-dependent charge Z(r) for a Fermi-distributed nucleus with rms 
        radius R and nuclear charge Z. The full nuclear potential is then given by V_nuc = - Z(r)/r; 
        a potential::Radial.Potential is returned.
"""
function fermiDistributedNucleus(Rrms::Float64, Z::Float64, grid::Radial.Grid)

    zz = zeros(Float64, grid.NoPoints);   zznew = zeros(Float64, grid.NoPoints);   dx = zeros(Float64, grid.NoPoints)
    ##x function  rho(r::Float64)  1.0 / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )   end
    function  r_rho(r::Float64)  r / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )   end
    function  rr_rho(r::Float64)  r^2 / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )  end

    if  Z < 1.2
        b = computeFermiBParameter(1.89);    b < 0  &&  error("Inappropriate R_rms radius.")
    else
        b = computeFermiBParameter(Rrms);    b < 0  &&  error("Inappropriate R_rms radius.")
    end
    
    fermiA_au = Defaults.convertUnits("length: from fm to atomic", fermiA)
    fermiC_au = Defaults.convertUnits("length: from fm to atomic", b)
    N = 1. / QuadGK.quadgk(rr_rho, 0., 1.3, rtol=1.0e-10)[1]
    #
    #
    for  i = 2:length(zznew)   
        zznew[i] = 1 / grid.r[i] * quadgk(rr_rho, 0.0, 1.0e-3, grid.r[i], rtol=1.0e-8)[1] 
        zznew[i] = zznew[i] + quadgk(r_rho, grid.r[i], 1.0, 1.0e+1, 1.0e+3, rtol=1.0e-8)[1]
        zznew[i] = Z * N * zznew[i] * grid.r[i]
    end 

    potential = Radial.Potential("nuclear-potential: Fermi-distributed", zznew, deepcopy(grid))
    return( potential )
end


"""
`Nuclear.pointNucleus(Z::Float64, grid::Radial.Grid)`  
    ... computes the effective, radial-dependent charge Z(r) for a point-like nucleus with nuclear charge Z. 
        The full nuclear potential is then given by V_nuc = - Z(r)/r; a potential::Radial.Potential is returned.
"""
function pointNucleus(Z::Float64, grid::Radial.Grid)

    zz = Z * ones(Float64, grid.NoPoints)

    potential = Radial.Potential("nuclear-potential: point-like", zz, deepcopy(grid))
    return( potential )
end


"""
`Nuclear.uniformNucleus(R::Float64, Z::Float64, grid::Radial.Grid)`  
    ... computes the effective, radial-dependent charge Z(r) for a uniformly-distributed nucleus with radius 
        R [in fm] and nuclear charge Z. The full nuclear potential is then given by V_nuc = - Z(r)/r outside 
        the nucleus; a potential::Radial.Potential is returned.
"""
function uniformNucleus(R::Float64, Z::Float64, grid::Radial.Grid)

    zz = zeros(Float64, grid.NoPoints);   R_au = Defaults.convertUnits("length: from fm to atomic", R)
    ##x println("uniformNucleus()::  R_au = $R_au")

    for i = 1:grid.NoPoints
        if     grid.r[i] <= R_au    zz[i] = Z / (2 * R_au) * (3. - grid.r[i]^2/R_au^2) * grid.r[i]
        ## if     grid.r[i] <= R_au    zz[i] = Z / (2 * R_au) * (grid.r[i]^2/R_au^2 - 3.) * grid.r[i] * grid.r[i]
        else   zz[i] = Z
        end
    end

    potential = Radial.Potential("nuclear-potential: uniformly-distributed", zz, deepcopy(grid))
    return( potential )
end


"""
`Nuclear.nuclearPotential(model::Nuclear.Model, grid::Radial.Grid)`  
    ... computes the effective, radial-dependent charge Z(r) for a nucleus with radius R and nuclear charge Z. 
        The full nuclear potential is then given by V_nuc = - Z(r)/r; a potential::Radial.Potential 
        is returned.
"""
function nuclearPotential(nm::Nuclear.Model, grid::Radial.Grid)
    zz = zeros(Float64, grid.NoPoints)

    if      nm.model == "point"      potential = Nuclear.pointNucleus(nm.Z, grid)
    elseif  nm.model == "uniform"    potential = Nuclear.uniformNucleus(nm.radius, nm.Z, grid)
    elseif  nm.model == "Fermi"      potential = Nuclear.fermiDistributedNucleus(nm.radius, nm.Z, grid)
    else    error("stop a")
    end

    return( potential )
    end


"""
`Nuclear.nuclearPotentialDH(model::Nuclear.Model, grid::Radial.Grid, lambda::Float64)`  
    ... computes the effective, radial-dependent charge Z(r) for a nucleus with radius R, nuclear charge Z 
        and for a Debye-Hueckel screening exp(-lambda r). The full Debye-Hueckel nuclear potential is then 
        given by V_nuc = - Z(r)/r; a potential::Radial.Potential is returned.
""" 
function nuclearPotentialDH(nm::Nuclear.Model, grid::Radial.Grid, lambda::Float64)
    pot = nuclearPotential(nm, grid)
    Zr  = pot.Zr
    for  i = 1:length(Zr)   Zr[i] = Zr[i] *exp(-lambda * pot.grid.r[i])    end
    
    potential = Radial.Potential(pot.name* "+ Debey-Hueckel screening", Zr, deepcopy(pot.grid))
    return( potential )
    end  

end # module

