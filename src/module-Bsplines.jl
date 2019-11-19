
"""
`module  JAC.Bsplines`  
    ... a submodel of JAC that contains all structs and methods to deal with the B-spline generation and the self-consistent field.
"""
module Bsplines

    using  Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..Nuclear, JAC
    ##x global JAC_counter = 0


    """
    `struct  Bsplines.Bspline`  ... defines a type for a single B-spline that is defined with regard to a well-defined grid.

        + i            ::Int64                     ... index of the B-spline
        + k            ::Int64                     ... order k of the B-spline
        + lower        ::Int64                     ... lower radial index (on the grid) from where the functions is generally nonzero.
        + upper        ::Int64                     ... upper radial index up to which the functions is generally nonzero.
        + bs           ::Array{Float64,1}          ... radial B-spline functions as defined on the predefined grid.
        + bp           ::Array{Float64,1}          ... derivative of bs on the predefined grid.
    """
    struct Bspline
        i              ::Int64    
        k              ::Int64     
        lower          ::Int64 
        upper          ::Int64 
        bs             ::Array{Float64,1}   
        bp             ::Array{Float64,1}   
    end


    """
    `struct  Bspline.Primitives`  ... defines a type for a set of primitive functions which typically belongs to a well-defined grid.

        + grid         ::Radial.Grid               ... radial grid on which the states are represented.
        + bsplinesL    ::Array{Bspline,1}          ... set of B-splines for the large components on the given grid.
        + bsplinesS    ::Array{Bspline,1}          ... set of B-splines for the small components on the given grid.
    """
    struct Primitives
        grid           ::Radial.Grid
        bsplinesL      ::Array{Bspline,1}
        bsplinesS      ::Array{Bspline,1}
    end


    # `Base.show(io::IO, primitives::Bsplines.Primitives)`  ... prepares a proper printout of the variable Bsplines.Primitives.
    function Base.show(io::IO, primitives::Bsplines.Primitives) 
        println(io, "grid:               $(primitives.grid)  ")
        println(io, "bsplinesL:          $(primitives.bsplinesL)  ")
        println(io, "bsplinesS:          $(primitives.bsplinesS)  ")
    end


    """
    `Bsplines.generateMatrix!(kappa::Int64, TTp::String, primitives::Bsplines.Primitives, storage::Dict{Array{Any,1},Array{Float64,2}})`  
         ... generates some (specified part of the Hamiltonian) matrix for the generalized eigenvalue problem and stores them into 
            the dictionary storage. The methods first tests of whether this matrix has been computed before alreay. Apart from the 
            overlap matrices (kappa = 0), all matrices belong to a particular single-electron angular momentum and to some matrix-block 
            of the overall eigenvalue problem. TTp = { 'LL', 'LS', 'SL', 'SS'} specifies the particular block, for which the B-splines 
            are given by primitivesT and primitivesTp. All B-splines are supposed to be defined with regard to the same (radial) grid; 
            an error message is issued if this is not the case. -- A matrix::Array{Float64,2} is returned which is quadratic for 
            'LL' and 'SS' and whose dimension depends on the number of B-splines for the large and small component, otherwise.
    """
    function generateMatrix!(kappa::Int64, TTp::String, primitives::Bsplines.Primitives, storage::Dict{Array{Any,1},Array{Float64,2}})
        # Look up the dictioniary of whether the requested matrix was calculated before
        wa = get( storage, [kappa, TTp], zeros(2,2) )
        if  wa != zeros(2,2)        
            return( wa )    
        end
        
        # Now calculate and store the requested matrix
        if      TTp == "LL-overlap"
            nsL = primitives.grid.nsL;    wa = zeros( nsL, nsL ) 
            for  i = 1:nsL
                for  j = 1:nsL
                    wa[i,j] = JAC.RadialIntegrals.overlap(primitives.bsplinesL[i], primitives.bsplinesL[j], primitives.grid)
                end
            end
        elseif  TTp == "SS-overlap"
            nsS = primitives.grid.nsS;    wa = zeros( nsS, nsS ) 
            for  i = 1:nsS
                for  j = 1:nsS
                    wa[i,j] = JAC.RadialIntegrals.overlap(primitives.bsplinesS[i], primitives.bsplinesS[j], primitives.grid)
                end
            end
        elseif  TTp == "LS-D_kappa^-"
            nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    wa = zeros( nsL, nsS ) 
            for  i = 1:nsL
                for  j = 1:nsS
                    wa[i,j] = Defaults.getDefaults("speed of light: c") *
                              JAC.RadialIntegrals.nondiagonalD(-1, kappa, primitives.bsplinesL[i], primitives.bsplinesS[j], primitives.grid)
                end
            end
        elseif  TTp == "SL-D_kappa^+"
            nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    wa = zeros( nsS, nsL ) 
            for  i = 1:nsS
                for  j = 1:nsL
                    wa[i,j] = Defaults.getDefaults("speed of light: c") *
                              JAC.RadialIntegrals.nondiagonalD( 1, kappa, primitives.bsplinesS[i], primitives.bsplinesL[j], primitives.grid)
                end
            end
        else   println("TTp = $TTp ");    error("stop a")
        end
        
        storage[ [kappa, TTp] ] = copy(wa)
        return( wa )
    end


    """
    `Bsplines.generateOrbitalFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::Bsplines.Primitives)`  
        ... generates the large and small components for the subshell sh from the primitives and their eigenvalues & eigenvectors. 
            A (normalized) orbital::Orbital is returned.
    """
    function generateOrbitalFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::Bsplines.Primitives)
        nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
        l  = Basics.subshell_l(sh);   ni = nsS + sh.n - l;          if   sh.kappa > 0   ni = ni + 1  end
        en = wc.values[ni];        if  en < 0.   isBound = true  else   isBound = false           end
        ev = wc.vectors[ni];       if  length(ev) != nsL + nsS    error("stop a")                 end
        
        P = zeros(primitives.grid.nr);    Pprime = zeros(primitives.grid.nr)
        Q = zeros(primitives.grid.nr);    Qprime = zeros(primitives.grid.nr)
        for  i = 1:nsL   
            for  j = primitives.bsplinesL[i].lower:primitives.bsplinesL[i].upper  P[j]      = P[j] + ev[i] * primitives.bsplinesL[i].bs[j]      end
            for  j = primitives.bsplinesL[i].lower:primitives.bsplinesL[i].upper  Pprime[j] = Pprime[j] + ev[i] * primitives.bsplinesL[i].bp[j] end
            ## for  j = 1:primitives.grid.nr  P[j]      = P[j] + ev[i] * primitives.bsplinesL[i].bs[j]      end
            ## for  j = 1:primitives.grid.nr  Pprime[j] = Pprime[j] + ev[i] * primitives.bsplinesL[i].bp[j] end
        end 
        for  i = 1:nsS   
            for  j = primitives.bsplinesS[i].lower:primitives.bsplinesS[i].upper  Q[j]      = Q[j] + ev[nsL+i] * primitives.bsplinesS[i].bs[j]      end
            for  j = primitives.bsplinesS[i].lower:primitives.bsplinesS[i].upper  Qprime[j] = Qprime[j] + ev[nsL+i] * primitives.bsplinesS[i].bp[j] end
            ## for  j = 1:primitives.grid.nr     Q[j]      = Q[j] + ev[nsL+i] * primitives.bsplinesS[i].bs[j]      end
            ## for  j = 1:primitives.grid.nr     Qprime[j] = Qprime[j] + ev[nsL+i] * primitives.bsplinesS[i].bp[j] end
        end 
        
        # Determine the maximum number of grid points for this orbital and normalized it propery
        mtp = 0;   for j = primitives.grid.nr:-1:1    if  abs(P[j]) + abs(Q[j]) > 2.0e-10   mtp = j;   break   end     end
        Px = zeros(mtp);    Qx = zeros(mtp);    Px[1:mtp] = P[1:mtp];    Qx[1:mtp] = Q[1:mtp]    
        Pprimex = zeros(mtp);    Qprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp];    Qprimex[1:mtp] = Qprime[1:mtp]    
        for  j = 1:mtp   if  abs(Px[j]) < 1.0e-10    Px[j] = 0.   end
                         if  abs(Qx[j]) < 1.0e-10    Qx[j] = 0.   end 
                         if  abs(Pprimex[j]) < 1.0e-10    Pprimex[j] = 0.   end
                         if  abs(Qprimex[j]) < 1.0e-10    Qprimex[j] = 0.   end      end
                         
        # Ensure that the large component of all orbitals start 'positive'
        wSign     = sum( Px[1:20] )
        if  wSign < 0.   Px[1:mtp] = -Px[1:mtp];   Qx[1:mtp] = -Qx[1:mtp]   
                         Pprimex[1:mtp] = -Pprimex[1:mtp];   Qprimex[1:mtp] = -Qprimex[1:mtp]   end
        orbital   = Orbital(sh, isBound, true, en, Px, Qx, Pprimex, Qprimex, Radial.Grid())
        
        # Renormalize the radial orbital   
        wN        = sqrt( JAC.RadialIntegrals.overlap(orbital, orbital, primitives.grid) )
        Px[1:mtp] = Px[1:mtp] / wN;    Qx[1:mtp] = Qx[1:mtp] / wN
        Pprimex[1:mtp] = Pprimex[1:mtp] / wN;    Qprimex[1:mtp] = Qprimex[1:mtp] / wN    
        
        return( Orbital(sh, isBound, true, en, Px, Qx, Pprimex, Qprimex, Radial.Grid()) )   
    end


    """
    `Bsplines.generateOrbitalFromPrimitives(energy::Float64, sh::Subshell, mtp::Int64, ev::Array{Float64,1}, primitives::Bsplines.Primitives)`  
         ... generates the large and small components for the subshell sh from the primitives and the eigenvector ev. 
            A (non-normalized) orbital::Orbital is returned.
    """
    function generateOrbitalFromPrimitives(energy::Float64, sh::Subshell, mtp::Int64, ev::Array{Float64,1}, primitives::Bsplines.Primitives)
        nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
        P = zeros(primitives.grid.nr);    Pprime = zeros(primitives.grid.nr)
        Q = zeros(primitives.grid.nr);    Qprime = zeros(primitives.grid.nr)
        for  i = 1:nsL   
            for  j = primitives.bsplinesL[i].lower:primitives.bsplinesL[i].upper  P[j]      = P[j] + ev[i] * primitives.bsplinesL[i].bs[j]      end
            for  j = primitives.bsplinesL[i].lower:primitives.bsplinesL[i].upper  Pprime[j] = Pprime[j] + ev[i] * primitives.bsplinesL[i].bp[j] end
            ## for  j = 1:primitives.grid.nr  P[j]      = P[j] + ev[i] * primitives.bsplinesL[i].bs[j]      end
            ## for  j = 1:primitives.grid.nr  Pprime[j] = Pprime[j] + ev[i] * primitives.bsplinesL[i].bp[j] end
        end 
        for  i = 1:nsS   
            for  j = primitives.bsplinesS[i].lower:primitives.bsplinesS[i].upper  Q[j]      = Q[j] + ev[nsL+i] * primitives.bsplinesS[i].bs[j]      end
            for  j = primitives.bsplinesS[i].lower:primitives.bsplinesS[i].upper  Qprime[j] = Qprime[j] + ev[nsL+i] * primitives.bsplinesS[i].bp[j] end
            ## for  j = 1:primitives.grid.nr     Q[j]      = Q[j] + ev[nsL+i] * primitives.bsplinesS[i].bs[j]      end
            ## for  j = 1:primitives.grid.nr     Qprime[j] = Qprime[j] + ev[nsL+i] * primitives.bsplinesS[i].bp[j] end
        end 
        
        ##x mtp = primitives.grid.nr - 40
        Px      = zeros(mtp);    Qx      = zeros(mtp);    Px[1:mtp]      = P[1:mtp];         Qx[1:mtp]      = Q[1:mtp]    
        Pprimex = zeros(mtp);    Qprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp];    Qprimex[1:mtp] = Qprime[1:mtp]    
        
        #== Determine the maximum number of grid points for this orbital and normalized it propery
        mtp = 0;   for j = primitives.grid.nr:-1:1    if  abs(P[j]) + abs(Q[j]) > 2.0e-10   mtp = j;   break   end     end
        Px = zeros(mtp);    Qx = zeros(mtp);    Px[1:mtp] = P[1:mtp];    Qx[1:mtp] = Q[1:mtp]    
        Pprimex = zeros(mtp);    Qprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp];    Qprimex[1:mtp] = Qprime[1:mtp]    
        for  j = 1:mtp   if  abs(Px[j]) < 1.0e-10    Px[j] = 0.   end
                         if  abs(Qx[j]) < 1.0e-10    Qx[j] = 0.   end 
                         if  abs(Pprimex[j]) < 1.0e-10    Pprimex[j] = 0.   end
                         if  abs(Qprimex[j]) < 1.0e-10    Qprimex[j] = 0.   end      
            println("j = $j     $(primitives.grid.r[j])   $(Px[j])   $(Qx[j]) ")                 
        end
                         
        # Ensure that the large component of all orbitals start 'positive'
        wSign     = sum( Px[1:20] )
        if  wSign < 0.   Px[1:mtp] = -Px[1:mtp];   Qx[1:mtp] = -Qx[1:mtp]   
                         Pprimex[1:mtp] = -Pprimex[1:mtp];   Qprimex[1:mtp] = -Qprimex[1:mtp]   end  ==#
        
        ##x println("nr = $(primitives.grid.nr),  mtp = $mtp")                 
      
        return( Orbital(sh, false, true, energy, Px, Qx, Pprimex, Qprimex, Radial.Grid()) )   
    end


    """
    `Bsplines.generatePrimitives(grid::Radial.Grid)`  
        ... generates the knots and the B-spline primitives of order k, both for the large and small components. 
            The function applies the given grid parameters; no primitive is defined beyond grid[n_max]. The definition of 
            the primitives follow the work of Zatsarinny and Froese Fischer, CPC 202 (2016) 287. --- A (set of) 
            primitives::Bsplines.Primitives is returned.
    """
    function generatePrimitives(grid::Radial.Grid)
        !(1 <= grid.orderL <= 11)   &&   error("Order should be 2 <= grid.orderL <= 11; obtained order = $(grid.orderL)")
        !(1 <= grid.orderS <= 11)   &&   error("Order should be 2 <= grid.orderS <= 11; obtained order = $(grid.orderS)")
        #
        # Now determined the B-splines on the grid for the large and small components
        primitivesL = Bsplines.Bspline[];   primitivesS = Bsplines.Bspline[];   lower = 0;   upper = 0
        bspline = zeros( grid.nr );    bprime = zeros( grid.nr )
        #
        for  i = 2:grid.nsL+1
            # Suppress a finite value at r = 0  by  j = 2:...
            for  j = 1:grid.nr     bspline[j] = Bsplines.generateSpline(i, grid.orderL, grid.t, grid.r[j])                  end
            if     bspline[1] > 0.   lower = 1
            else   for  j = 2:grid.nr       if  bspline[j] > 0.   lower = j-1;   break   end  end    
            end
            for  j = grid.nr-1:-1:1    
               if     bspline[j] <= 0.     continue
               else                        upper = j+1;   break
               end
            end
            # 
            for  j = lower:upper   bprime[j] = Bsplines.generateSplinePrime(i, grid.orderL, grid.t, grid.r[j])    end
            #
            push!( primitivesL, Bsplines.Bspline( i, grid.orderL, lower, upper, copy(bspline), copy(bprime) ) )
        end
        #
        for  i = 2:grid.nsS+1
            # Suppress a finite value at r = 0  by  j = 2:...
            for  j = 1:grid.nr     bspline[j] = Bsplines.generateSpline(i, grid.orderS, grid.t, grid.r[j])                  end
            if     bspline[1] > 0.   lower = 1
            else   for  j = 2:grid.nr       if  bspline[j] > 0.   lower = j-1;   break   end  end    
            end
            for  j = grid.nr-1:-1:1    
               if     bspline[j] <= 0.     continue
               else                        upper = j+1;   break
               end
            end
            # 
            for  j = lower:upper   bprime[j] = Bsplines.generateSplinePrime(i, grid.orderS, grid.t, grid.r[j])    end
            #
            push!( primitivesS, Bsplines.Bspline( i, grid.orderS, lower, upper, copy(bspline), copy(bprime) ) )
        end
        
        return( Bsplines.Primitives(grid, primitivesL, primitivesS) )
    end



    """
    `Bsplines.generateSpline(i::Int64, k::Int64, tlist::Array{Float64,1}, r::Float64)`  
        ... generates the i-th B-spline of order k at the radial point r; this generation is carried out recursively. 
            A value::Float64 is returned.
    """
    function generateSpline(i::Int64, k::Int64, tlist::Array{Float64,1}, r::Float64)
        if      i == 1   &&   0  ==  r                      return( 1.0 )    
        elseif  k == 1
                if       tlist[i] <  r <= tlist[i+1]        return( 1.0 )   else    return( 0. )    end
        elseif  k > 1
                wa = 0.;    wci = (tlist[i+k-1] - tlist[i]);    wcip = (tlist[i+k] - tlist[i+1])
                if wci  != 0    wa = wa + (r - tlist[i])/(tlist[i+k-1] - tlist[i])   * generateSpline(i, k-1, tlist, r)     end
                if wcip != 0    wa = wa + (tlist[i+k] - r)/(tlist[i+k] - tlist[i+1]) * generateSpline(i+1, k-1, tlist, r)   end
                return( wa )
        else    error("stop a")
        end
    end



    """
    `Bsplines.generateSplinePrime(i::Int64, k::Int64, tlist::Array{Float64,1}, r::Float64)` 
        ... generates the first derivative of the i-th B-spline of order k at the radial point r; it uses the values of related
            B-splines at the same point. A value::Float64 is returned.
    """
    function generateSplinePrime(i::Int64, k::Int64, tlist::Array{Float64,1}, r::Float64)
        wa = 0.
        if  (tlist[i+k-1] - tlist[i]) == 0. 
            Defaults.warn(AddWarning, "DB for i=$i, k=$k, tlist[i+k-1]=$(tlist[i+k-1]), tlist[i]=$(tlist[i]), bspline=$(generateSpline(i, k-1, tlist, r)), wa=$wa")  
        else
            wa = wa + (k - 1)/(tlist[i+k-1] - tlist[i]) * generateSpline(i,   k-1, tlist, r)
        end

        
        if  (tlist[i+k] - tlist[i+1]) == 0.  
            Defaults.warn(AddWarning, "DB for i=$i, k=$k, tlist[i+k]=$(tlist[i+k]), tlist[i+1]=$(tlist[i+1]), bspline=$(generateSpline(i+1, k-1, tlist, r)), wa=$wa")   
        else
            wa = wa - (k - 1)/(tlist[i+k] - tlist[i+1]) * generateSpline(i+1, k-1, tlist, r)
        end

        return( wa )
    end



    """
    `Bsplines.setupLocalMatrix(kappa::Int64, primitives::Bsplines.Primitives, pot::Radial.Potential, storage::Dict{Array{Any,1},Array{Float64,2}}) 
         ...setp-up the local parts of the generalized eigenvalue problem for the symmetry block kappa and the given (local) potential. 
            The B-spline (basis) functions are defined by primitivesL for the large component and primitivesS for the small one, respectively.
    """
    function setupLocalMatrix(kappa::Int64, primitives::Bsplines.Primitives, pot::Radial.Potential, storage::Dict{Array{Any,1},Array{Float64,2}})
        nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
        wa  = zeros( nsL+nsS, nsL+nsS );   wb = zeros( nsL+nsS, nsL+nsS )
        # (1) Compute or fetch the diagonal 'overlap' blocks
        wb[1:nsL,1:nsL]                 = generateMatrix!(0, "LL-overlap", primitives, storage)
        wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = generateMatrix!(0, "SS-overlap", primitives, storage)
        # (2) Re-compute the diagonal blocks for the local potential
        for  i = 1:nsL
            for  j = 1:nsL   
                wa[i,j] = JAC.RadialIntegrals.Vlocal(primitives.bsplinesL[i], primitives.bsplinesL[j], pot, primitives.grid)
            end
        end
        for  i = 1:nsS    
            for  j = 1:nsS    
                wa[nsL+i,nsL+j] = JAC.RadialIntegrals.Vlocal(primitives.bsplinesS[i], primitives.bsplinesS[j], pot, primitives.grid)
                ##x if  i >= j  && abs(wa[nsL+i,nsL+j]) > 1.0e-5   println("i = $i, j = $j, wa = $(wa[nsL+i,nsL+j]) ")   end
            end
        end
        ##x error("stop for test reasons")
        # (3) Substract the rest mass from the 'SS' block
        wa[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = wa[nsL+1:nsL+nsS,nsL+1:nsL+nsS] - 2 * Defaults.getDefaults("speed of light: c")^2 * wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS]
        # (4) Compute or fetch the diagonal 'D_kappa' blocks
        wa[1:nsL,nsL+1:nsL+nsS] = wa[1:nsL,nsL+1:nsL+nsS] + generateMatrix!(kappa, "LS-D_kappa^-", primitives, storage)
        wa[nsL+1:nsL+nsS,1:nsL] = wa[nsL+1:nsL+nsS,1:nsL] + generateMatrix!(kappa, "SL-D_kappa^+", primitives, storage)
        
        # Test for 'real-symmetric matrix' ... and symmetrize otherwise
        nx = 0
        for  i = 1:nsL+nsS    
            for  j = i+1:nsL+nsS    
                if  abs(  (wa[i,j] - wa[j,i])/(wa[i,j] + wa[j,i]) ) > 1.0e-7   nx = nx + 1    end
                wij = (wa[i,j] + wa[j,i]) / 2
                wa[i,j] = wa[j,i] = wij
            end
        end
        ny = (nsL+nsS)^2/2 - (nsL+nsS)
        if  nx > 0    Defaults.warn(AddWarning, "setupLocalMatrix:: $nx (from $(ny)) non-symmetric H-matrix integrals " *
                                           "for kappa = $kappa with relative deviation > 1.0e-7.")                 end
       
        return( wa )
    end



    """
    `Bsplines.generateGalerkinMatrix(primitives::Bsplines.Primitives, energy::Float64, sh::Subshell, pot::Radial.Potential)`  
        ... generates the Galerkin-A matrix for the given potential and B-spline primitives; a matrix::Array{Float64,2} is returned.
    """
    function generateGalerkinMatrix(primitives::Bsplines.Primitives, energy::Float64, sh::Subshell, pot::Radial.Potential)

        # Define the storage for the calculations of matrices; this is necessary to use the Bsplines.generateMatrix!() function
        println("(Re-) Define a storage array for various B-spline matrices:")
        storage  = Dict{Array{Any,1},Array{Float64,2}}()
        # Set-up the overlap matrix
        nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    wb  = zeros( nsL+nsS, nsL+nsS )
        wb[1:nsL,1:nsL]                 = generateMatrix!(0, "LL-overlap", primitives, storage)
        wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = generateMatrix!(0, "SS-overlap", primitives, storage)
        # Set-up the local Hamiltonian matrix
        wa = setupLocalMatrix(sh.kappa, primitives, pot::Radial.Potential, storage::Dict{Array{Any,1},Array{Float64,2}})
        wa[1:end,1:end] = wa[1:end,1:end] - energy * wb[1:end,1:end]
        ##x wa[1:nsL+nsS,1:nsL+nsS] = wa[1:nsL+nsS,1:nsL+nsS] - energy * wb[1:end,1:end]

        return( wa )
    end


    #==  """
    `Bsplines.generateOrbitalContinuum(subshell::Subshell, energy::Float64, potential::Radial.Potential, grid::Radial.Grid, mtp::Int64, 
                                           basis::Basis)`  
        ... generates a continuum orbital for given subshell, energy and potential (effective charge Z(r)) within a B-spline basis. 
            This orbital is constructed by setting up and solving a linear system of equations and orthogonal to all (bound) 
            orbitals in the given basis. For this, the one-particle Dirac equation is re-written as matrix equation in 
            order to determine coefficients in a B-spline basis. A continuum orbital::Radial.Orbital is returned.
    """
    function generateOrbitalContinuum(subshell::Subshell, energy::Float64, potential::Radial.Potential, grid::Radial.Grid, mtp::Int64, basis::Basis)
        println("Now prepare a continuum orbital for $sh ")
        N = 19;   K = 0
        waL = Bsplines.generatePrimitives(4, N, grid )
        waS = Bsplines.generatePrimitives(5, N, grid )
        whe = zeros(N-1, N-1);   wd = zeros(N-1)
        # First set up the matrix and overlap for the bound orbitals of the same symmetry
        for   sh in basis.subshells
            if  sh.kappa == subshell.kappa 
               K = K + 1
               ## for  n = 1:N-1   whe[K, n] =  
            end
        end 

        ## solution = LAPACK.gesv!(whe, wed)
        return( nothing )
    end  ==#



    """
    `Bsplines.generateOrbitalsForPotential(primitives::Bsplines.Primitives, kappa::Int64, subshells::Array{Subshell,1}, 
                                           pot::Radial.Potential; printout::Bool=true)`  
        ... generates all single-electron orbitals from subhsells, which are of symmetry kappa. 
            A set of orbitals::Dict{Subshell, Orbital} is returned.
    """
    function generateOrbitalsForPotential(primitives::Bsplines.Primitives, kappa::Int64, subshells::Array{Subshell,1}, 
                                          pot::Radial.Potential; printout::Bool=true)
        orbitals = Dict{Subshell, Orbital}()
        if  kappa == 0     return( orbitals )     end       # This should not cause any problem.
        
        # Define the storage for the calculations of matrices; this is necessary to use the Bsplines.generateMatrix!() function.
        if  printout    println("(Re-) Define a storage array for various B-spline matrices:")    end
        storage  = Dict{Array{Any,1},Array{Float64,2}}()
        # Set-up the overlap matrix
        nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
        wb = zeros( nsL+nsS, nsL+nsS )
        
        # Compute or fetch the diagonal 'overlap' blocks
        wb[1:nsL,1:nsL]                 = generateMatrix!(0, "LL-overlap", primitives, storage)
        wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = generateMatrix!(0, "SS-overlap", primitives, storage)
        
        # Compute the local Hamiltonian matrix and diagonalize it
        wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
        w2 = Basics.diagonalize("generalized eigenvalues: Julia, eigfact", wa, wb)
        nsi = nsS;    if kappa > 0   nsi = nsi + 1   end
        if  printout   Basics.tabulateKappaSymmetryEnergiesDirac(kappa, w2.values, nsi, Nuclear.Model(4.))    end
        
        # Collect all the computed single-electron orbitals
        for  sh in subshells
            orbitals = Base.merge( orbitals, Dict( sh => generateOrbitalFromPrimitives(sh, w2, primitives) ))
        end
        
        return( orbitals )
    end



    """
    `Bsplines.generateOrbitalsHydrogenic(primitives::Bsplines.Primitives, nuclearModel::Nuclear.Model, 
                                             subshells::Array{Subshell,1}; printout::Bool=true)`  
        ... generates the hydrogenic orbitals for the nuclear model and for the given subshells. A set of 
            orbitals::Dict{Subshell, Orbital}() is returned.
    """
    function generateOrbitalsHydrogenic(primitives::Bsplines.Primitives, nuclearModel::Nuclear.Model, subshells::Array{Subshell,1}; 
                                        printout::Bool=true) 

        newOrbitals = Dict{Subshell, Orbital}()
        # Define the storage for the calculations of matrices; this is necessary to use the Bsplines.generateMatrix!() function
        if  printout    println("(Re-) Define a storage array for various B-spline matrices:")    end
        storage  = Dict{Array{Any,1},Array{Float64,2}}()
        # Set-up the overlap matrix
        nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
        wb = zeros( nsL+nsS, nsL+nsS )
        
        # Compute or fetch the diagonal 'overlap' blocks
        wb[1:nsL,1:nsL]                 = generateMatrix!(0, "LL-overlap", primitives, storage)
        wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = generateMatrix!(0, "SS-overlap", primitives, storage)
        
        # Test for 'real-symmetric matrix' ... and symmetrize otherwise
        for  i = 1:nsL+nsS    
            for  j = i+1:nsL+nsS    
                if  abs(  (wb[i,j] - wb[j,i])/wb[i,j] ) > 1.0e-8  
                   Defaults.warn(AddWarning, "Nonsymmetric overlap at $i, $j with rel. D = $( abs((wb[i,j] - wb[j,i])/wb[i,j]) ) ")  end
                wij = (wb[i,j] + wb[j,i]) / 2
                wb[i,j] = wb[j,i] = wij
            end
        end
        
        # Define the nuclear potential for the give nuclear model
        if  printout    println("Nuclear model = $(nuclearModel) ")    end
        if       nuclearModel.model == "point"
            pot = Nuclear.pointNucleus(nuclearModel.Z, primitives.grid)
        elseif   nuclearModel.model == "Fermi"
            pot = Nuclear.fermiDistributedNucleus(nuclearModel.radius, nuclearModel.Z, primitives.grid) 
        elseif   nuclearModel.model == "uniform"
            pot = Nuclear.uniformNucleus(nuclearModel.radius, nuclearModel.Z, primitives.grid)
        else
            error("stop a")
        end

        alreadyDone = falses( length(subshells) )
        for  i = 1:length(subshells)
            if  alreadyDone[i]   continue   end
            # 
            sh = subshells[i]
            if  printout    println("Generate hydrogenic orbital for subshell $sh ")    end
            wa = Bsplines.setupLocalMatrix(sh.kappa, primitives, pot, storage)
            w2 = Basics.diagonalize("generalized eigenvalues: Julia, eigfact", wa, wb)
            nsi = nsS;    if sh.kappa > 0   nsi = nsi + 1   end
            if  printout   Basics.tabulateKappaSymmetryEnergiesDirac(sh.kappa, w2.values, nsi, nuclearModel)    end
            newOrbitals[sh] = generateOrbitalFromPrimitives(sh, w2, primitives)
            # Take over orbitals of the same symmetry
            for  j = 1:length(subshells)
                if  alreadyDone[j]   continue
                elseif  sh.kappa == subshells[j].kappa   
                    newOrbitals[subshells[j]] = generateOrbitalFromPrimitives(subshells[j], w2, primitives)
                    if  printout   println(">> Use hydrogenic orbital from this symmetry block also for $(subshells[j]).")    end
                    alreadyDone[j] = true
                end
            end
        end

        return( newOrbitals )
    end


    """
    `Bsplines.solveSelfConsistentFieldOptimized(primitives::Bsplines.Primitives, nuclearModel::Nuclear.Model, basis::Basis,
                                                    settings::AsfSettings; printout::Bool=true)` 
        ... solves the self-consistent field for ...     A basis::Basis is returned.
    """
    function solveSelfConsistentFieldOptimized(primitives::Bsplines.Primitives, nuclearModel::Nuclear.Model, basis::Basis, 
                                               settings::AsfSettings; printout::Bool=true) 
        Defaults.setDefaults("standard grid", primitives.grid)
        
        #
        #
        #

        newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, previousOrbitals)
        return( newBasis )
    end
    


    """
    `Bsplines.solveSelfConsistentFieldMean(primitives::Bsplines.Primitives, nuclearModel::Nuclear.Model, basis::Basis, 
                                               settings::AsfSettings; printout::Bool=true)` 
        ... solves the self-consistent field for a given local mean-field potential as specified by the settings::AsfSettings. 
            A basis::Basis is returned.
    """
    function solveSelfConsistentFieldMean(primitives::Bsplines.Primitives, nuclearModel::Nuclear.Model, basis::Basis, 
                                          settings::AsfSettings; printout::Bool=true) 
        Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
        # Define the storage for the calculations of matrices
        if  printout    println("(Re-) Define a storage array for various B-spline matrices:")    end
        storage  = Dict{Array{Any,1},Array{Float64,2}}()
        
        # Set-up the overlap matrix; compute or fetch the diagonal 'overlap' blocks
        grid = primitives.grid;   nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
        wb = zeros( nsL+nsS, nsL+nsS )
        wb[1:nsL,1:nsL]                 = generateMatrix!(0, "LL-overlap", primitives, storage)
        wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = generateMatrix!(0, "SS-overlap", primitives, storage)
        # Determine the symmetry block of this basis and define storage for the kappa blocks and orbitals from the last iteration
        wNotYet = trues( length(basis.subshells) );   kappas = Int64[]
        for  k = 1:length(basis.subshells)
            sh = basis.subshells[k];   
            if  wNotYet[k]   
                push!(kappas, sh.kappa);    wNotYet[k] = false
                for kx = 1:length(basis.subshells)   if wNotYet[kx]   &&   sh.kappa == basis.subshells[kx].kappa   wNotYet[kx] = false  end   end
            end
        end 
        bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(basis.orbitals)
        for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end
        # Determine te nuclear potential once at the beginning
        nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, grid)
        
        # Start the SCF procedure for all symmetries
        isNotSCF = true;   NoIteration = 0
        while  isNotSCF
            NoIteration = NoIteration + 1;   go_on = false 
            if  NoIteration >  settings.maxIterationsScf
                println("Maximum number of SCF iterations = $(settings.maxIterationsScf) is reached ... computations proceed.")
                break
            end
            if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
            #
            for kappa in kappas
                # (1) First re-define an (arbitrary) 'level' that represents the mean occupation for the local potential
                wBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, previousOrbitals)
                NoCsf  = length(wBasis.csfs)
                wmc    = zeros( NoCsf );   wN = 0.
                for i = 1:NoCsf   wmc[i] = AngularMomentum.twoJ(wBasis.csfs[i].J) + 1.0;   wN = wN + abs(wmc[i])^2    end
                for i = 1:NoCsf   wmc[i] = wmc[i] / sqrt(wN)   end
                wLevel = Level( AngularJ64(0), AngularM64(0), Basics.plus, 0, -1., 0., true, wBasis, wmc)
                # (2) Re-compute the local potential
                ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
                ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
                ## wp3 = compute("radial potential: Kohn-Sham", grid, wLevel)
                ## wp4 = compute("radial potential: Dirac-Fock-Slater", grid, wLevel)
                ##x println("\n Compare potentials \n   core-Hartree    Hartree-Slater    Kohn-Sham    Dirac-Fock-Slater")
                ##x for  mx = 1:10:length(wp1.V)
                ##x     @printf("%.8e  %.8e  %.8e  %.8e  \n", wp1.V[mx], wp2.V[mx], wp3.V[mx],wp4.V[mx])
                ##x end
                if       settings.methodScf == "meanHS"     wp = compute("radial potential: Hartree-Slater",    grid, wLevel)
                elseif   settings.methodScf == "meanDFS"    wp = compute("radial potential: Dirac-Fock-Slater", grid, wLevel)
                else     error("stop potential")
                end
                pot = Basics.add(nuclearPotential, wp)
                # (3) Set-up the diagonal part of the Hamiltonian matrix
                wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
                # (4) Solve the generalized eigenvalue problem
                wc = Basics.diagonalize("generalized eigenvalues: Julia, eigfact", wa, wb)
                # (5) Analyse and print information about the convergence of the symmetry blocks and the occupied orbitals
                wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
                if  wcBlock > 1.000 * settings.accuracyScf   go_on = true   end
                for  sh in basis.subshells
                    if      sh in settings.frozenSubshells   ## do nothing
                    elseif  sh.kappa == kappa
                        newOrbital = generateOrbitalFromPrimitives(sh, wc, primitives)
                        wcOrbital  = Basics.analyzeConvergence(previousOrbitals[sh], newOrbital)
                        if  wcOrbital > settings.accuracyScf   go_on = true   end
                           sa = "  $sh::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "  
                           sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                           sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                           if  printout    println(sa)    end
                        ## println("  $sh  en [a.u.] = $(newOrbital.energy)   self-consistency = $(wcOrbital), $(wcBlock) [kappa=$kappa] ") 
                        previousOrbitals[sh] = newOrbital
                    end
                end
                # (6) Re-define the bsplineBlock
                bsplineBlock[kappa] = wc
            end
            if  go_on   nothing   else   break   end
        end

        newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, previousOrbitals)
        return( newBasis )
    end


end # module
