
"""
`module  JAC.Radial`  ... a submodel of JAC that contains all structs and methods to deal with the radial grid, radial orbitals and potentials; 
                          it is using JAC.
"""
module Radial

    using  JAC, QuadGK, Printf
    export Grid, Potential, Orbital, SingleElecSpectrum


    """
    `struct  Radial.Grid`  ... defines a type for the radial grid which contains all information about the grid parameters, the genration 
                               of the B-spline basis as well as for performing radial integrations.

        ** Physical grid parameter **
        + rnt            ::Float64           ... smalles grid point > 0.
        + h              ::Float64           ... stepsize in the construction of the exponential grid.
        + hp             ::Float64           ... asymptotic stepsize of the log-lin grid.
        + NoPoints       ::Int64             ... No. of grid points.
        ** B-spline grid parameters and break points **
        + t              ::Array{Float64,1}  ... radial break points for the B-splines
        + nt             ::Int64             ... number of break points in the t-grid.
        + nth            ::Int64             ... take each nth-point from the 'physical grid' as break point.
        + orderL         ::Int64             ... B-spline order of large components
        + orderS         ::Int64             ... B-spline order of small components
        + nsL            ::Int64             ... number of B-splines for large components
        + nsS            ::Int64             ... number of B-splines for small components
        + orderGL        ::Int64             ... order of the Gauss-Lengendre integration if mesh == MeshGL
        ** Radial mesh points **
        + mesh           ::RadialMesh
        + nr             ::Int64             ... number of mesh points in the r-grid.
        + r              ::Array{Float64,1}  ... radial grid points
        + rp             ::Array{Float64,1}  ... derivative of the radial grid at the grid points
        + rpor           ::Array{Float64,1}  ... rp over r
        + wr             ::Array{Float64,1}  ... integration weights for all grid points, for instance, GL weights.
    """
    struct  Grid
        rnt              ::Float64
        h                ::Float64
        hp               ::Float64
        NoPoints         ::Int64
        #
        t                ::Array{Float64,1}
        nt               ::Int64 
        nth              ::Int64 
        orderL           ::Int64
        orderS           ::Int64
        nsL              ::Int64
        nsS              ::Int64
        orderGL          ::Int64
        #
        mesh             ::RadialMesh
        nr               ::Int64 
        r                ::Array{Float64,1}
        rp               ::Array{Float64,1}
        rpor             ::Array{Float64,1}
        wr               ::Array{Float64,1}
    end


    """
    `JAC.Radial.Grid()`  ... constructor to define an 'empty' instance of the radial grid.

    `JAC.Radial.Grid("grid: exponential")`  ... constructor to define an instance of the (Grasp) standard exponential grid.

    `JAC.Radial.Grid("grid: by given parameters"; rnt = rntx, h = hx, hp = hpx, NoPoints = NoPointx)`  ... constructor to define an instance 
                                                 of the radial grid with the given parameters.
    """
    function Grid()
        Grid(0., 0., 0., 0,   Float64[], 0, 0, 0, 0, 0, 0, 0,   MeshGrasp, 0, Float64[], Float64[], Float64[], Float64[] )
    end


    function Grid(sa::String; rnt::Float64 = 2.0e-6, h::Float64 = 5.0e-2, hp::Float64 = 0., NoPoints::Int64 = 390)
        
        # Define the general parameters for the B-spline grid and radial mesh
        mesh = MeshGL  ## MeshGL   MeshGrasp
        nth = 7;    orderL = 7;    orderS = 8;   orderGL = 7;    nsL = 0;   nsS = 0
        nr = NoPoints + 10;    r = zeros(nr);    rp = zeros(nr);    rpor = zeros(nr)

        # Now define the physical grid due to the given parameters: either an exponential or exponential-linear grid
        if  hp == 0.    &&    sa in [ "grid: exponential", "grid: by given parameters"]
            rp[1] = rnt;   eph = exp(h);    ett = 1.0
           
            for  i  in 2:nr
                ett     = eph * ett
                ettm1   = ett - 1.0
                r[i]    = rnt * ettm1
                rp[i]   = rnt * ett
                rpor[i] = ett/ettm1
            end 

        elseif  hp != 0.    &&    sa in [ "grid: by given parameters"]
            function rprime(r :: Float64)
                return( 1/ (1/(r + rnt) + h/hp) )
            end
        
            function f(r :: Float64, i :: Int)
                return( log( r/rnt + 1) + h/hp * r - (i - 1) * h )
            end
        
            function fprime(r :: Float64)
                return( 1. / (r + rnt) + h/hp )
            end
        
            r[1]    = 0.
            rp[1]   = rprime(0.) 
            rpor[1] = 0.
        
            rc = 0.
            rn = 0.
        
            for i = 2:nr
                rn += rnt
                while ((abs((rc - rn)/rn) > 100 * eps(Float64))) 
                    rc = rn
                    rn = rc - f(rc, i) / fprime(rc)
                end
                r[i]    = rn
                rp[i]   = rprime(rn)
                rpor[i] = rp[i] / r[i]
            end
        else
            error("Unsupported keystring = $sa for given parameters.")
        end
        
        # Define the radial break points and the number of B-splines
        nt = 0;    nrmax = 0;    t = zeros(orderL)
        for  i = 1:nth:nr
            nrmax = i;    nt = nt + 1;    push!( t, r[i])
        end
        nsL = nt + orderL - 1 - 8;    nsS = nt + orderS - 1 - 8;    nr = nrmax - 10
        for   ip = nt+1:nt+max(orderL, orderS)-1    push!( t, r[nt])   end
        ##x ns = nv + k - 1       # No. of splines
        ##x for  ip = ii+1:ii+k-1    push!( tlist, grid.r[i])   end
        
        # Define the radial grid due to the definition
        wr = zeros(nr)
        if  mesh == MeshGL
            nr = 1;   r = Float64[ 0.];   rp = Float64[];   rpor = Float64[];   wr = Float64[0.];    rlow = 0.
            wax = gauss(orderGL);   ra = wax[1];   wa = wax[2] 
            for  i = 1:length(t)
                rstep = t[i] - rlow
                if  rstep > 0
                    bma = t[i] - rlow;   bpa = t[i] + rlow 
                    for  j = 1:orderGL
                       nr = nr + 1;    push!(r, ra[j] * bma/2. + bpa/2. );    push!(wr, wa[j] * bma/2. )
                    end
                end
                rlow = t[i]
            end
            #
            println("Define a radial grid of type $mesh with $nr grid points")
            println(" [rnt=" * @sprintf("%.3e",rnt) * ", h=" * @sprintf("%.3e",h) * 
                    ", hp=" * @sprintf("%.3e",hp) * ", NoPoints=$NoPoints, r_max=" * @sprintf("%.3e",r[nr]) * ";")
            println("  B-splines wit break points at every $(nth)th point, nsL=$nsL, nsS=$nsS, orderL=$orderL, orderS=$orderS, orderGL=$orderGL] ")
            ## for j = 1:nr   println("j = $j,  r = $(r[j]),  wr = $(wr[j])")   end
        end

        Grid(rnt, h, hp, NoPoints,  t, nt, nth, orderL, orderS, nsL, nsS, orderGL,   mesh, nr, r, rp, rpor, wr)
    end


    """
    `Base.show(io::IO, grid::Radial.Grid)`  ... prepares a proper printout of the variable grid::Radial.Grid.
    """
    function Base.show(io::IO, grid::Radial.Grid) 
    
        sa = Base.string(grid::Radial.Grid);    print(io, sa * "\n")
    
        nr = grid.nr;    nt = grid.nt;    if  nr < 6   return( nothing )    end
        print(io, "r:    ", grid.r[1:3],    "  ...  ", grid.r[nr-2:nr],     "\n") 
        if  grid.mesh == MeshGrasp
            print(io, "rp:   ", grid.rp[1:3],   "  ...  ", grid.rp[nr-2:nr],    "\n") 
            print(io, "rpor: ", grid.rpor[1:3], "  ...  ", grid.rpor[nr-2:nr],  "\n") 
        end
        print(io, "wr:   ", grid.wr[1:3],   "  ...  ", grid.wr[nr-2:nr],    "\n") 
        print(io, "t:    ", grid.t[1:3],    "  ...  ", grid.t[nt-2:nt]          ) 
    end


    """
    `Base.string(grid::Radial.Grid)`  ... provides a String notation for the variable grid::Radial.Grid.
    """
    function Base.string(grid::Radial.Grid) 
        if   grid.NoPoints == 0
            sa = "Radial grid not defined;  NoPoints = $(grid.NoPoints) ..."
        else
            sa = "Radial grid:  rnt = $(grid.rnt),  h = $(grid.h),  hp = $(grid.hp),  NoPoints = $(grid.NoPoints), "
            sa = sa * "orderL = $(grid.orderL),  orderS = $(grid.orderS),  nsL = $(grid.nsL),  nsS = $(grid.nsS),  mesh = $(grid.mesh), ...  "
        end 
    end


    """
    `struct  Radial.GridGL`  ... defines a type for Gauss-Legendre grid of given order but where the interval is hard-coded due to
                                 the given keystring

        + nt           :Int64                 ... number of mesh points in the grid.
        + t            ::Array{Float64,1}     ... mesh points in the variable t.
        + wt           ::Array{Float64,1}     ... weights to the mesh points in t.
    """
    struct  GridGL
        nt             ::Int64  
        t              ::Array{Float64,1}
        wt             ::Array{Float64,1}
    end


    """
    `JAC.Radial.GridGL("QED", orderGL::Int64; printout::Bool=false)`  ... constructor to define Gauss-Legendre grid 
         for the typical QED computation in the interval [1.0, infinity].
    """
    function GridGL(sa::String, orderGL::Int64; printout::Bool=false)
        !(sa == "QED")  && error("Unrecognized keystring; sa = $sa")
        txlow = 1.;    t = Float64[];    wt = Float64[];    nt = 0
        for i = 1:100000
            # Define the exponential increase (1.5) and the maximum size (infinity=150.)
            txup = txlow * 1.5;    if  txup > 150.   break   end
            #
            wax = QuadGK.gauss(orderGL);   tx = wax[1];   wtx = wax[2] 
            bma  = txup - txlow;   bpa = txup + txlow
            for  j = 1:orderGL
                nt = nt + 1;    push!(t, tx[j] * bma/2. + bpa/2. );    push!(wt, wtx[j] * bma/2. )
            end
            txlow = txup
        end
        if  printout   println("Gauss-Legendre grid with $nt mesh points from t = 1 ... $txlow.")   end
           
        GridGL(nt, t, wt)
    end


    """
    `Base.show(io::IO, grid::GridGL)`  ... prepares a proper printout of the variable grid::GridGL.
    """
    function Base.show(io::IO, grid::GridGL) 
        print(io, "Gauss-Legendre grid with $(grid.nt) mesh points: \n")
        print(io, "t:   ", grid.t[1:5],    "  ...  ", grid.t[grid.nt-4:grid.nt],    "\n") 
        print(io, "w:   ", grid.wt[1:5],   "  ...  ", grid.wt[grid.nt-4:grid.nt],    "\n") 
    end
    


    """
    `struct  Radial.Potential`  ... defines a struct for the radial potential which contains all information about its physical content.

        + name           ::String            ... A name for the potential.
        + V              ::Array{Float64,1}  ... radial potential function.
        + grid           ::RadialGrid        ... radial grid on which the potential is defined.
    """
    struct  Potential
        name             ::String
        V                ::Array{Float64,1}
        grid             ::Radial.Grid
    end


    """
    `JAC.Radial.Potential()`  ... constructor to define an 'empty' instance of the radial potential.
    """
    function Potential()
        Potential("", Float64[], JAC.Radial.Grid())
    end


    """
    `Base.show(io::IO, potential::Radial.Potential)`  ... prepares a proper printout of the variable potential::Radial.Potential.
    """
    function Base.show(io::IO, potential::Radial.Potential) 
    
        sa = Base.string(potential);    print(io, sa * "\n")
    
        n = length(potential.V);                    if  n < 6   return( nothing )    end
        print(io, "V:    ", potential.V[1:23],    "  ...  ", potential.V[n-22:n],    "\n") 
        print(io, potential.grid)
    end


    """
    `Base.string(potential::Radial.Potential)`  ... provides a String notation for the variable potential::Radial.Potential.
    """
    function Base.string(potential::Radial.Potential) 
        if  length(potential.V)  == 0
            sa = "Radial potential not yet defined; kind = $(potential.name) ..."
        else
            sa = "$(potential.name) (radial) potential ... defined on $(length(potential.V)) grid points ..."
        end 
    end


    """
    `struct  Orbital`  ... defines a type for a single-electron radial orbital function with a large and small component, and which can refer to
             either the standard or an explicitly given grid due to the logical flag isStandardGrid. Bound-state orbitals with energy < 0 are 
             distinguished from free-electron orbitals by the flag isBound. -- Note that the arrays P, Q and grid cannot be defined only by the 
             standard constructor but are typically set explicitly after an instance of this type has been created.

        + subshell        ::Subshell          ... Relativistic subshell.
        + isBound         ::Bool              ... Logical flag to distinguish between bound (true) and free-electron orbitals (false).
        + useStandardGrid ::Bool              ... Logical flag for using the standard grid (true) or an explicitly given grid (false).
        + energy          ::Float64           ... Single-electron energies of bound orbitals are always negative.
        + P               ::Array{Float64,1}  ... Large and ..
        + Q               ::Array{Float64,1}  ... small component of the radial orbital.
        + Pprime          ::Array{Float64,1}  ... dP/dr.
        + Qprime          ::Array{Float64,1}  ... dQ/dr.
        + grid            ::Array{Float64,1}  ... explic. defined radial grid array for P, Q, if StandardGrid = false.
    """
    mutable struct Orbital
        subshell          ::Subshell
        isBound           ::Bool             
        useStandardGrid   ::Bool
        energy            ::Float64 
        P                 ::Array{Float64,1} 
        Q                 ::Array{Float64,1}
        Pprime            ::Array{Float64,1}
        Qprime            ::Array{Float64,1}
        grid              ::Radial.Grid
    end


    """
    `JAC.Orbital(subshell::Subshell, energy::Float64)`  ... constructor for given subshell and energy, and where isStandardGrid is set to true;
         the grid must be defined explicitly and the large and small components are not yet defined in this case.
    """
    function Orbital(subshell::Subshell, energy::Float64)
        if energy < 0    isBound = true    else    isBound = false    end
        isStandardGrid = true
        P = Array{Float64,1}[];    Q = Array{Float64,1}[];    grid = RadialGrid()

        Orbital(subshell, isBound, isStandardGrid, energy, P, Q, grid)
    end


    """
    `JAC.Orbital(label::String, energy::Float64)`  ... constructor for given string identifier and energy, and where isStandardGrid is set 
         to true; the grid must be defined explicitly and the large and small components are not yet defined in this case.
    """
    function Orbital(label::String, energy::Float64)
        if energy < 0    isBound = true    else    isBound = false    end
    
        subshell = Subshell(label);    isStandardGrid = true
        P = Array{Float64,1}[];        Q = Array{Float64,1}[];    grid = RadialGrid()

        Orbital(subshell, isBound, isStandardGrid, energy, P, Q, grid)
    end


    """
    `Base.show(io::IO, orbital::Orbital)`  ... prepares a proper printout of the variable orbital::Orbital.
    """
    function Base.show(io::IO, orbital::Orbital) 
        n = length(orbital.P)

        if   orbital.useStandardGrid
            stdgrid = JAC.give("standard grid")
            stdgrid.NoPoints == 0     &&   return( print("Standard grid has not yet been defined.") )

            n > stdgrid.NoPoints      &&   error("length of P does not match to standard grid")    
            n != length(orbital.Q)    &&   error("P and Q have different length")  
        else  
            !(n == length(orbital.Q) == length(orbital.Q))   &&    error("P, Q, grid have different length")
        end
    
        sa = Base.string(orbital::Orbital);    print(io, sa * "\n")

        if n <= 6
            print(io, "Large component P: ", orbital.P[1:end], "\n") 
            print(io, "Small component Q: ", orbital.Q[1:end], "\n") 
        else 
            n = length(orbital.P)
            print(io, "Large component P: ", orbital.P[1:25], "  ...  ", orbital.P[n-24:n], "\n") 
            print(io, "Small component Q: ", orbital.Q[1:25], "  ...  ", orbital.P[n-24:n], "\n")
        end

        if orbital.useStandardGrid 
           n <= 6   &&   print(io, "Defined on Grid:   ", stdgrid.r[1:n], "\n") 
           n >  6   &&   print(io, "Defined on Grid:   ", stdgrid.r[1:3], "  ...  ", stdgrid.r[n-3:n], "\n") 
        else
           n <= 6   &&   print(io, "Defined on Grid:   ", orbital.grid.r[1:n], "\n") 
           n > 6    &&   print(io, "Defined on Grid:   ", orbital.grid.r[1:5], "  ...  ", orbital.grid.r[n-4:n], "\n") 
        end
    end


    """
    `Base.string(orbital::Orbital)`  ... provides a String notation for the variable orbital::Orbital with just basic information.
    """
    function Base.string(orbital::Orbital) 
        if  orbital.isBound           sa = "Bound-state orbital "     else    sa = "Free-electron orbital "         end
        if  orbital.useStandardGrid   sb = "the standard grid: "      else    sb = "an explicitly-defined grid: "   end
        energy = orbital.energy
        sc = string(orbital.subshell) * " with energy $energy a.u. ";   n = length(orbital.P)

        return( sa * sc * "is defined with $n (grid) points on " * sb )
    end


    """
    `struct  SingleSymOrbitals`  ... defines a type for a (relativistic and quasi-complete) single-electron spectrum for symmetry kappa
             with N positive and/or N negative states. All these states are defined with regard to the same grid.

        + name         ::String                ... A name for this spectrum, may contain information about the original type of basis functions.
        + kappa        ::Int64                 ... symmetry of the one-electron spectrum
        + NoStates     ::Int64                 ... Number of positive and negative states (if onlyPositive = true)
        + onlyPositive ::Bool                  ... True if only the positive part is kept.
        + pOrbitals    ::Array{Orbital,1}      ... Positive-energy orbitals states, in increasing order.
        + nOrbitals    ::Array{Orbital,1}      ... Negative-energy orbitals states, in increasing order.
        + grid         ::RadialGrid            ... radial grid on which the states are represented.
    """
    struct SingleSymOrbitals
        name           ::String
        kappa          ::Int64
        NoStates       ::Int64
        onlyPositive   ::Bool
        pOrbitals      ::Array{Orbital,1} 
        nOrbitals      ::Array{Orbital,1} 
        grid           ::Radial.Grid
    end


    """
     `JAC.Radial.SingleSymOrbitals()`  ... constructor for providing an 'empty' instance of this struct.
    """
    function SingleSymOrbitals()
        SingleSymOrbitals("", 0, 0, true, Orbital[], Orbital[], Radial.Grid() )
    end


    """
    `Base.show(io::IO, spectrum::SingleElecSpectrum)`  ... prepares a proper printout of the variable spectrum::SingleElecSpectrum.
    """
    function Base.show(io::IO, spectrum::SingleSymOrbitals) 
        sa = "Single-electron spectrum $(spectrum.name) for symmetry kappa=$(spectrum.kappa) with $(spectrum.NoStates) states "
        if  spectrum.onlyPositive   sa = sa * "from just the positive spectrum."
        else                        sa = sa * "from the positive and negative spectrum."
        end
        println(io, sa) 
   end


    """
    `JAC.Radial.OrbitalBunge1993(subshell::Subshell, Z::Int64)`  ... to calculate the radial orbital on the standard grid.
    """
    function OrbitalBunge1993(subshell::Subshell, Z::Int64, grid::Radial.Grid)
        #x grid = JAC.give("standard grid");    
        #x qn = JAC.SubshellQuantumNumbers( string(subshell) );    n = qn[1];    kappa = qn[2];   l = qn[3]
        n = subshell.n;    kappa = subshell.kappa;   l = JAC.subshell_l(subshell)
        ##x println("n, l = $n, $l")
    
        wa = JAC.store("orbital functions: NR, Bunge (1993)", Z)
        if      l == 0  &&  n > wa[5][1]+0     error("orbitals with s symmetry for Z = $Z only available up to $(wa[5][1]+0) only.")
        elseif  l == 1  &&  n > wa[6][1]+1     error("orbitals with p symmetry for Z = $Z only available up to $(wa[6][1]+1) only.")
        elseif  l == 2  &&  n > wa[7][1]+2     error("orbitals with d symmetry for Z = $Z only available up to $(wa[7][1]+2) only.")
        elseif  l == 3  &&  n > wa[8][1]+3     error("orbitals with f symmetry for Z = $Z only available up to $(wa[8][1]+3) only.")
        elseif  l == 4  &&  n > wa[9][1]+4     error("orbitals with g symmetry for Z = $Z only available up to $(wa[9][1]+4) only.")
        end
    
        P = Float64[]
        for  k in 1:grid.NoPoints-5
            Rnl = 0.; r = grid.r[k]
            for i in 4:length(wa[5+l])
                njl = wa[5+l][i][1] 
                Zjl = wa[5+l][i][2]
                Cjl = wa[5+l][i][2+n-l]
                Njl = (2*Zjl)^(njl+0.5) / sqrt( factorial(2*njl) ) 
                Rnl = Rnl + Cjl*Njl*r^(njl-1)*exp(-Zjl*r);
            end
            ##x println("OrbitalBunge1993-aa: k, r, Rnl = $k  $r  $Rnl ")
            push!(P, Rnl*r )
        end
 
        # Now define small component by apllying kinetic-balance condition
        Q = zeros(size(P, 1))
        #x println("small component not yet defined properly; P = $P ")
    
        dP(i) = JAC.Math.derivative(P, i)
        for i = 2:size(Q, 1)
            Q[i] = -1/(2 * JAC.INVERSE_FINE_STRUCTURE_CONSTANT) * (dP(i) / grid.h / grid.rp[i] + kappa/grid.r[i]) * P[i]
        end
        ##x println("OrbitalBunge1993-ab: ",Q)
    
        o     = Orbital( subshell, true, true, -0.5, P, Q, grid)
        norma = JAC.RadialIntegrals.overlap(o, o, grid)
        o.P = o.P/sqrt(norma)
        o.Q = o.Q/sqrt(norma)

        normb = JAC.RadialIntegrals.overlap(o, o, grid)
        ##x println("OrbitalBunge1993-ac:  for subshell $subshell : norm-before = $norma, norm-after = $normb")
    
        return( o )
    end


    """
    `JAC.Radial.OrbitalMcLean1981(subshell::Subshell,Z::Int64)`  ... to calculate the radial orbital on the standard grid.  
    """
    function compute_McLean1981(subshell::Subshell,Z::Int64)
        grid = JAC.give("standard grid");    
        #x qn = JAC.SubshellQuantumNumbers( string(subshell) );    n = qn[1];   l = qn[3]
        n = subshell.n;   l = subshell_l( subshell )
    
        wa = JAC.store("orbital functions: NR, McLean (1981)", Z)

        !(55 <= Z <= 92)   &&                     error("Nuclear charge must be 55 <= (Z=$Z) <= 92.")            
        wa == Any[]        &&                     error("No data available for neutral Z = $Z.")

        if      l == 0     &&  n > wa[5][1]+0     error("Orbitals with s symmetry for Z = $Z only available up to $(wa[5][1]+0) only.")
        elseif  l == 1     &&  n > wa[6][1]+1     error("orbitals with p symmetry for Z = $Z only available up to $(wa[6][1]+1) only.")
        elseif  l == 2     &&  n > wa[7][1]+2     error("orbitals with d symmetry for Z = $Z only available up to $(wa[7][1]+2) only.")
        elseif  l == 3     &&  n > wa[8][1]+3     error("orbitals with f symmetry for Z = $Z only available up to $(wa[8][1]+3) only.")
        elseif  l == 4     &&  n > wa[9][1]+4     error("orbitals with g symmetry for Z = $Z only available up to $(wa[9][1]+4) only.")
        end
    
        P = Float64[]
        for  k in 1:grid.NoPoints-5
            Rnl = 0.; r = grid.r[k]
            for i in 4:length(wa[5+l])
                njl = wa[5+l][i][1] 
                Zjl = wa[5+l][i][2]
                Cjl = wa[5+l][i][2+n-l]
                Njl = (2*Zjl)^(njl+0.5) / sqrt( factorial(2*njl) ) 
                Rnl = Rnl + Cjl*Njl*r^(njl-1)*exp(-Zjl*r);
            end
            println("k, r, Rnl = $k  $r  $Rnl ")
            push!(P, Rnl*r )
        end
 
        println("small component not yet defined properly; P = $P ")

        return( Orbital( subshell, -0.5) )
    end


    """
    `JAC.Radial.OrbitalPrimitiveSlater(subshell::Subshell, N::Int64, alpha0::Float64, beta0::Float64, grid::Radial.Grid)`  ... to calculate a 
         list of (non-relativistic) radial Slater primitives::Array{Orbital,1} for the subshell-symmetry on the given grid. All orbitals have 
         the same subshell, isBound = true, useStandardGrid = true and energy = 0. The small components are constructed by applying kinetic 
         balance to the large component, and all the primitives are properly normalized.
    """
    function OrbitalPrimitiveSlater(subshell::Subshell, N::Int64, alpha0::Float64, beta0::Float64, grid::Radial.Grid) 
        function P(i::Int64, l::Int64, eta::Float64)   return( grid.r[i]^(l+1) * exp(-eta*grid.r[i]) )   end
        dP(i) = JAC.Math.derivative(Px, i)

        orbitals = Orbital[];    ll = subshell_l(subshell)
        for  i = 1:N
            eta = alpha0 * beta0^(i-1)
            # Px  = [for i = 1:grid.NoPoints P(i, ll, eta)]
            Qx  = zeros(size(Px, 1))
            for j = 2:size(Q, 1)
                Qx[j] = -1/(2 * JAC.INVERSE_FINE_STRUCTURE_CONSTANT) * (dP(j) / grid.h / grid.rp[j] + kappa/grid.r[j]) * Px[j]
            end
    
            o     = Orbital( subshell, true, true, 0., Px, Qx, grid() )
            norma = JAC.RadialIntegrals.overlap(o, o, grid)
            o.P = o.P/sqrt(norma)
            o.Q = o.Q/sqrt(norma)

            normb = JAC.RadialIntegrals.overlap(o, o, grid)
            println("OrbitalPrimitiveSlater-aa:  for i = $i : norm-before = $norma, norm-after = $normb")
            
            push!( orbitals, o)
        end

        return( SingleElecSpectrum() )
    end



    """
    `JAC.Radial.determineZbar(pot::Radial.Potential)`  ... determines the effective charge that is asymptotically seen by the electron
                                                           in the potential pot. A Zbar::Float64 is returned.
    """
    function determineZbar(pot::Radial.Potential) 
        nx = 5   # Number of grid points for determining the effective charge Zbar
        mtp = size( pot.V, 1);    grid = pot.grid;    meanZbar = zeros(nx);    devsZbar = zeros(nx);   ny = 0
        ##x println("mtp = $mtp,  pot.V = $(pot.V[mtp-nx+1:mtp])")
        for  i = mtp-nx+1:mtp
            ny  = ny + 1;    meanZbar[ny] = pot.V[mtp]
        end
        mZbar   = sum(meanZbar) / nx;   for  i = 1:nx    devsZbar[i] = (meanZbar[i] - mZbar)^2    end
        stdZbar = sqrt( sum(devsZbar) / nx )
        println("Radial potential with effective charge Zbar=" * @sprintf("%.4e",mZbar) *
                " (Delta-Zbar=" * @sprintf("%.4e",stdZbar) * ")." )
        
        return( mZbar )
    end

end # module
