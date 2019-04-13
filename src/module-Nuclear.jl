
"""
`module JAC.Nuclear`  ... a submodel of JAC that contains procedures for defining a Nuclear.Model and for calculating various nuclear
                          potentials; it is using JAC, JAC.Radial   
                          ##x  ./../deps/bin/libnuc.so.
"""
module Nuclear
  
    using Interact, QuadGK, JAC, JAC.Radial
    ##x const LIBNUC = joinpath(@__DIR__, "..", "deps", "bin", "libnuc.so")


    export Model


    """
    `struct  Nuclear.Model`  ... defines a type for the nuclear model, i.e. for its form and parameters.

        + Z        ::Float64         ... nuclear charge
        + model    ::String          ... identifier of the nuclear model: {"Fermi", "Point", "Uniform"}
        + mass     ::Float64         ... atomic mass
        + radius   ::Float64         ... (root-mean square) radius of a uniform or Fermi-distributed nucleus
        + spinI    ::AngularJ64      ... nuclear spin I, must be >= 0
        + mu       ::Float64         ... magnetic dipole moment in Bohr magnetons
        + Q        ::Float64         ... electric qadrupole moment
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
    `JAC.Nuclear.Model(Z::Real)`  ... constructor just for a given nuclear charge Z, and where a Fermi model is defined 
         with a = 1.0 and c = 1.0 for the moment. Both, the nuclear spin and moments are all set to zero in this case.
    """
    function Model(Z::Real)
        Z < 0.1  &&  error("Z must be >= 0.1")
        model    = "Fermi"
        mass     = 2*Z + 0.005*Z^2
        radius   = JAC.Nuclear.Rrms(mass)
        spinI    = AngularJ64(0)
        mu       = 0.
        Q        = 0.

        Model(Z, model, mass, radius, spinI, mu, Q) 
    end


    """
    `JAC.Nuclear.Model(Z::Real, model::String)`  ... constructor just for a given nuclear charge Z and
         model = {"Fermi", "point", "uniform"}, and where further parameters are defined approximately. 
         Both, the nuclear spin and moments are all set to zero in this case.
    """
    function Model(Z::Real, model::String)
        Z < 0.1                                    &&  error("Z must be >= 0.1")
        !(model in ["Fermi", "point", "uniform"])  && error("Inappropriate model = $model")
        if  model == "point"
            mass     = 1000*Z
            radius   = 0.
        else
            mass     = 2*Z + 0.005*Z^2
            radius   = JAC.Nuclear.Rrms(mass)
        end
        spinI    = AngularJ64(0)
        mu       = 0.
        Q        = 0.

        Model(Z, model, mass, radius, spinI, mu, Q) 
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
    `JAC.Nuclear.Model(gui::Guint; model::Nuclear.Model=Model(36.0))`  ... constructor that is defined by a graphical user interface.
    """
    function Model(gui::Guint; model::Nuclear.Model=Model(36.0))
        nmd = model
        
        if      gui == Gui
            t1 = "Nuclear model with: "
            b1 = slider(1:110, label = "charge", value = nmd.Z)
            b2 = dropdown([nmd.model, "Fermi", "point", "uniform"])
            b3 = spinbox(label="mass  ";   value=nmd.mass)
            b4 = spinbox(label="radius";   value=nmd.radius)
            b5 = spinbox(label="2*spin";   value=nmd.spinI.num)
            b6 = spinbox(label="mu    ";   value=nmd.mu)
            b7 = spinbox(label="Q     ";   value=nmd.Q)
            update = button("Update")
            ui = vbox( hbox( pad(0em, t1) ),
                       hbox( pad(1em, b1) ), 
                       hbox( pad(1em, b2), pad(1em, b3), pad(1em, b4) ), 
                       hbox( pad(1em, b5), pad(1em, b6), pad(1em, b7), pad(1em, update) )
                      )
            Interact.display(ui) 
            output = Interact.@map  (&update; Nuclear.Model( observe(b1)[], observe(b2)[], observe(b3)[], observe(b4)[], 
                                                             AngularJ64( Int64(observe(b5)[])//2 ), observe(b6)[], observe(b7)[] )  )
            return( output )

        else  error("Unsupported Guint = $gui.")
        end
    end
  

    """
   `JAC.Nuclear.fermiA`  ... provides a value::Float64 for the fermi_a parameter.
    """
    fermiA   = 2.3/(4 * log(3))
  
  
    """
    `JAC.Nuclear.Rrms(A)`  ... provides a value::Float64 for the root-mean-squared radius R of a uniformly-distributed charge density of 
                               a nuclues with mass A.
    """
    function Rrms(A)
        return( 0.836 * A^(1/3) + 0.57 )
    end
  
  
    """
    `JAC.Nuclear.fermiDensity(r, b)`  ... provides a value::Float64 for a fermi-distributed charge density at r and for given fermi_a
                                          and fermi_b parameters.
    """
    function fermiDensity(r, b)
        return( 1/(1 + exp((r-b)/fermiA)) )
    end
  

    """
    `JAC.Nuclear.fermiRrms(b::Float64)`  ... provides a value::Float64 for the root-mean-squared radius of a fermi-distributed charge 
                                             density for given fermi_a and fermi_b parameters.
    """
    function fermiRrms(b::Float64)
    
        return( sqrt(12 * fermiA^2 * JAC.Math.polylogExp(b/fermiA, 5) / JAC.Math.polylogExp(b/fermiA, 3)) )
    end
  

    """
    `JAC.Nuclear.computeFermiBParameter(R::Float64)`  ... computes a value::Float64 for the fermi_b parameter for a fermi-distributed 
                                                          nuclear with root-mean square radius R.
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
    `JAC.Nuclear.fermiDistributedNucleus(Rrms::Float64, Z::Float64, grid::Radial.Grid)`  ... computes the effective, radial-dependent charge Z(r) 
         for a Fermi-distributed nucleus with rms radius R and nuclear charge Z. The full nuclear potential is then given by V_nuc = - Z(r)/r; 
         a potential::Radial.Potential is returned.
    """
    function fermiDistributedNucleus(Rrms::Float64, Z::Float64, grid::Radial.Grid)
        # grid = JAC.JAC_STANDARD_GRID
    
        zz = zeros(Float64, grid.nr);   zznew = zeros(Float64, grid.nr);   dx = zeros(Float64, grid.nr)
        ##x function  rho(r::Float64)  1.0 / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )   end
        function  r_rho(r::Float64)  r / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )   end
        function  rr_rho(r::Float64)  r^2 / (1.0 + exp( (r-fermiC_au)/fermiA_au ) )  end
    
        b = computeFermiBParameter(Rrms);    b < 0  &&  error("Inappropriate R_rms radius.")
        
        fermiA_au = JAC.convert("length: from fm to atomic", fermiA)
        fermiC_au = JAC.convert("length: from fm to atomic", b)
        ##x fermi1_au = JAC.convert("length: from fm to atomic", 1.0)
        ##x N = 3/(4pi * fermiC_au^3) / (1 + pi^2 * fermiA_au^2 / fermiC_au^2)
        N = 1. / QuadGK.quadgk(rr_rho, 0., 1.3, rtol=1.0e-10)[1]
        ##x println("fermiA_au = $fermiA_au  fermiC_au = $fermiC_au  fermi1_au = $fermi1_au")
        ##x wa = N * QuadGK.quadgk(rr_rho, 0., 1.3, rtol=1.0e-10)[1]
        ##x println("N = $N  N-int = $wa " )
        #
        #
        for  i = 2:length(zznew)   
            zznew[i] = 1 / grid.r[i] * quadgk(rr_rho, 0.0, 1.0e-3, grid.r[i], rtol=1.0e-8)[1] 
            zznew[i] = zznew[i] + quadgk(r_rho, grid.r[i], 1.0, 1.0e+1, 1.0e+3, rtol=1.0e-8)[1]
            zznew[i] = Z * N * zznew[i] * grid.r[i]
        end 
    
        #== ccall((:__nuc_MOD_nucpotmain, LIBNUC), Nothing, 
              (Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
               Ref(Int32(grid.nr)), Ref(grid.rnt), Ref(grid.h), Ref(grid.hp), zz, Ref(b), Ref(fermiA), Ref(Z))
              
        for i = 1:3:length(zznew)
            dx[i] = zz[i] / zznew[i]
            println("++ $i )  old = $(zz[i])   new = $(zznew[i])   dx = $(dx[i])")
        end  ==#
        ##x error("stop -- nuclear potential") 
   
        potential = Radial.Potential("nuclear-potential: Fermi-distributed", zznew, deepcopy(grid))
        return( potential )
    end

  
    """
    `JAC.Nuclear.pointNucleus(Z::Float64, grid::Radial.Grid)`  ... computes the effective, radial-dependent charge Z(r) for a point-like nucleus 
         with nuclear charge Z. The full nuclear potential is then given by V_nuc = - Z(r)/r; a potential::Radial.Potential is returned.
    """
    function pointNucleus(Z::Float64, grid::Radial.Grid)
        # grid = JAC.JAC_STANDARD_GRID
    
        zz = Z * ones(Float64, grid.nr)
    
        potential = Radial.Potential("nuclear-potential: point-like", zz, deepcopy(grid))
        return( potential )
    end

  
    """
    `JAC.Nuclear.uniformNucleus(R::Float64, Z::Float64, grid::Radial.Grid)`  ... computes the effective, radial-dependent charge Z(r) for a 
         uniformly-distributed nucleus with radius R [in fm] and nuclear charge Z. The full nuclear potential is then given by V_nuc = - Z(r)/r outside 
         the nucleus; a potential::Radial.Potential is returned.
    """
    function uniformNucleus(R::Float64, Z::Float64, grid::Radial.Grid)
        # grid = JAC.JAC_STANDARD_GRID
    
        zz = zeros(Float64, grid.nr);   R_au = JAC.convert("length: from fm to atomic", R)
        ##x println("uniformNucleus()::  R_au = $R_au")
    
        for i = 1:grid.nr
            if     grid.r[i] <= R_au    zz[i] = Z / (2 * R_au) * (grid.r[i]^2/R_au^2 - 3.) * grid.r[i] * grid.r[i]
            else   zz[i] = Z
            end
        end
    
        potential = Radial.Potential("nuclear-potential: uniformly-distributed", zz, deepcopy(grid))
        return( potential )
    end

  
    """
    `JAC.Nuclear.nuclearPotential(model::Nuclear.Model, grid::Radial.Grid)`  ... computes the effective, radial-dependent charge Z(r) for 
         a nucleus with radius R and nuclear charge Z. The full nuclear potential is then given by V_nuc = - Z(r)/r; a potential::Radial.Potential 
         is returned.
    """
    function nuclearPotential(nm::Nuclear.Model, grid::Radial.Grid)
        zz = zeros(Float64, grid.nr)
    
        if      nm.model == "point"      potential = JAC.Nuclear.pointNucleus(nm.Z, grid)
        elseif  nm.model == "uniform"    potential = JAC.Nuclear.uniformNucleus(nm.radius, nm.Z, grid)
        elseif  nm.model == "Fermi"      potential = JAC.Nuclear.fermiDistributedNucleus(nm.radius, nm.Z, grid)
        else    error("stop a")
        end

        return( potential )
     end

    
    """
    `JAC.Nuclear.nuclearPotentialDH(model::Nuclear.Model, grid::Radial.Grid, lambda::Float64)`  ... computes the effective, 
         radial-dependent charge Z(r) for a nucleus with radius R, nuclear charge Z and for a Debye-Hueckel screening exp(-lambda r). 
         The full Debye-Hueckel nuclear potential is then given by V_nuc = - Z(r)/r; a potential::Radial.Potential is returned.
    """ 
    function nuclearPotentialDH(nm::Nuclear.Model, grid::Radial.Grid, lambda::Float64)
        pot = nuclearPotential(nm, grid)
        Zr  = pot.Zr
        for  i = 1:length(Zr)   Zr[i] = Zr[i] *exp(-lambda * pot.grid.r[i])    end
        
        potential = Radial.Potential(pot.name* "+ Debey-Hueckel screening", Zr, deepcopy(pot.grid))
        return( potential )
     end  

end # module

