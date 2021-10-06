
using    GSL, HypergeometricFunctions, Plots, Printf, SpecialFunctions, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength, ..Radial, ..ManyElectron, 
            ..Nuclear, ..Pulse, ..TableStrings

"""
`StrongField.HydrogenPnl(epsiloni::Float64, n::Int, l::Int, rGrid::Array{Float64,1})`  
    ... returns the non-relativistic hydrogen-like radial wave function with modified binding energy epsiloni
            at the radial grid points rGrid
"""
function  HydrogenPnl(epsiloni::Float64, n::Int, l::Int, rGrid::Array{Float64,1})
    Pnl = Float64[]
    Z = n*sqrt(-2*epsiloni)
    
    for j = 1:length(rGrid)
        r = rGrid[j]
        p = r * sqrt( (2*Z/n)^3 * factorial(n-l-1)/(2*n*factorial(n+l)) ) * exp(-Z*r/n) * (2*Z*r/n)^l * GSL.sf_laguerre_n(n-l-1,2*l+1,2*Z*r/n)
        push!( Pnl, p )
    end
    
    return( Pnl )
end

"""
`StrongField.HydrogenDPnlDr(epsiloni::Float64, n::Int, l::Int, rGrid::Array{Float64,1})`  
    ... returns the r-derivative of the non-relativistic hydrogen-like radial wave function with modified binding energy epsiloni
            at the radial grid points rGrid
"""
function  HydrogenDPnlDr(epsiloni::Float64, n::Int, l::Int, rGrid::Array{Float64,1})
    DPnl = Float64[]
    Z = n*sqrt(-2*epsiloni)
    
    for j = 1:length(rGrid)
        r = rGrid[j]
        if  n-l-2 >= 0
            p = 1/n * 2^(l+1) * exp(-r*Z/n) * (r*Z/n)^l * sqrt( Z^3 * factorial(n-l-1) / (n^4*factorial(l+n)) ) * ( (n+l*n-r*Z) * GSL.sf_laguerre_n(n-l-1,2*l+1,2*r*Z/n) - 2*r*Z * GSL.sf_laguerre_n(n-l-2,2+2*l,2*r*Z/n) )
        elseif  n == 1 && l == 0
            p = -2 * exp(-r*Z) * sqrt(Z^3) * (r*Z - 1)
        else #not implemented
            p = 0
        end
        push!( DPnl, p )
    end
    
    return( DPnl )
end


#---------------------------------------------------------------------------------------------------
#--------------------------------------FORMULATION in l-ml-BASIS:-----------------------------------
#---------------------------------------------------------------------------------------------------

"""
`StrongField.pReducedMEHydrogenicUncoupled(epsilonp::Float64, lp::Int, n::Int, l::Int, epsiloni::Float64, volkov::AbstractVolkovState)`  
    ... computes the reduced matrix elements of the momentum operator <epsilonp lp ||p||n l> in the one-particle picture
"""
function  pReducedMEHydrogenicUncoupled(epsilonp::Float64, lp::Int, n::Int, l::Int, epsiloni::Float64, volkov::AbstractVolkovState)
    rmax = 100.
    orderGL = 1000
    gaussLegendre = Radial.GridGL("Finite",0.0,rmax,orderGL)
    rgrid = gaussLegendre.t
    weights = gaussLegendre.wt
    
    Pnl = HydrogenPnl( epsiloni, n, l, rgrid )
    if  typeof(volkov) == FreeVolkov
        Pepsplp = VolkovP( epsilonp, lp, rgrid )
    elseif  typeof(volkov) == CoulombVolkov
        Pepsplp = CoulombVolkovP( epsilonp, lp, volkov.Z, rgrid )
    end
    
    DPnl = HydrogenDPnlDr( epsiloni, n, l, rgrid )
    
    integral = 0. * im
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGL
        r = rgrid[j]
        integrand = conj( Pepsplp[j] )/r * ( r*DPnl[j] - ((lp-l)*(lp+l+1))/2 * Pnl[j] )

        #Gauss-Legendre sum
        integral = integral + weights[j] * integrand
    end

    #Note that GSL.sf_coupling_3j takes takes the input (2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)
    integral = integral * (-im)^(lp+1) * (-1)^lp * GSL.sf_coupling_3j( 2*lp, 2*1, 2*l, 0, 0, 0 )
    
    return( integral )
end


"""
`StrongField.scalarProdBoundContHydrogenicUncoupled(epsilonp::Float64, lp::Int, n::Int, l::Int, epsiloni::Float64, volkov::AbstractVolkovState)`  
    ... computes the scalar product of the bound (hydrogenic) and continuum (Volkov) states in the one-particle picture
"""
function  scalarProdBoundContHydrogenicUncoupled(epsilonp::Float64, lp::Int, n::Int, l::Int, m::Int, epsiloni::Float64, volkov::AbstractVolkovState)
    rmax = 1000.
    orderGL = 1000
    gaussLegendre = Radial.GridGL("Finite",0.0,rmax,orderGL)
    rgrid = gaussLegendre.t
    weights = gaussLegendre.wt
    
    Pnl = HydrogenPnl( epsiloni, n, l, rgrid )
    
    if  typeof(volkov) == FreeVolkov
        Pepsplp = VolkovP( epsilonp, lp, rgrid )
    elseif  typeof(volkov) == CoulombVolkov
        Pepsplp = CoulombVolkovP( epsilonp, lp, volkov.Z, rgrid )
    end
    
    DPnl = HydrogenDPnlDr( epsiloni, n, l, rgrid )
    
    integral = 0. * im
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGL
        r = rgrid[j]
        integrand = conj( Pepsplp[j] ) * Pnl[j]

        #Gauss-Legendre sum
        integral = integral + weights[j] * integrand
    end
    
    return((-im)^l * integral)
end


#---------------------------------------------------------------------------------------------------
#--------------------------------------FORMULATION in j-l-mj-BASIS:---------------------------------
#---------------------------------------------------------------------------------------------------

"""
`StrongField.pReducedMEHydrogenic(Pepsplp::Array{ComplexF64,1} , lp::Int, jp::Float64, epsiloni::Float64, n::Int, l::Int, j::Int, rGrid::Radial.Grid )`  
    ... computes the reduced matrix elements of the momentum operator <epsilonp lp jp||p||n l j> in the one-particle picture for hydrogenic initial states, 
        where the radial wave function Pepsplp of |epsilonp lp jp> needs to be provided as an argument on the grid rGrid
"""
function  pReducedMEHydrogenic(Pepsplp::Array{ComplexF64,1}, lp::Int, jp::Float64, epsiloni::Float64, n::Int, l::Int, j::Float64, rGrid::Radial.Grid )

    if rGrid.meshType != Radial.MeshGL()
        println("StrongField module needs a radial grid of meshType Radial.MeshGL() to perform the radial integrals.")
        return(0.)
    end

    #---- Compute the factor in <epsilonp lp jp||p||n l j> = fac * <epsilonp lp||p||n l>-----
    fac = 0.
    
    for    kms = -1:0 #ms = 1/2 + kms
        ms = 1/2 + kms
        for    ml = -l:l
            mj = ml + ms 
            for    mlp = -lp:lp
                mjp = mlp + ms 
                fac = fac + AngularMomentum.ClebschGordan_old(  AngularJ64(Rational(lp)),   AngularM64(Rational(mlp)) , 
                                                                AngularJ64(1//2),           AngularM64(Rational(ms)), 
                                                                AngularJ64(Rational(jp)),   AngularM64(Rational(mjp))       )       * 
                            AngularMomentum.ClebschGordan_old(  AngularJ64(Rational(l)),    AngularM64(Rational(ml)),  
                                                                AngularJ64(1//2),           AngularM64(Rational(ms)), 
                                                                AngularJ64(Rational(j)),    AngularM64(Rational(mj))       )    
            end
        end
    end
    #----------------------------------------------------------------------------------------
    

    #--------------------------------Compute <epsilonp lp||p||n l>---------------------------
    rgrid = rGrid.r
    orderGL = size(rgrid)[1]
    weights = rGrid.wr
    
    Pnl = HydrogenPnl( epsiloni, n, l, rgrid )
    DPnl = HydrogenDPnlDr( epsiloni, n, l, rgrid )
    
    integral = 0. * im
    
    #Sum over grid and compute Gauss-Legendre sum
    for    k = 1:orderGL
        r = rgrid[k]
        integrand = conj( Pepsplp[k] )/r * ( r*DPnl[k] - ((lp-l)*(lp+l+1))/2 * Pnl[k] )

        #Gauss-Legendre sum
        integral = integral + weights[k] * integrand
    end

    #Note that GSL.sf_coupling_3j takes takes the input (2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)
    integral = integral * (-im)^(lp+1) * (-1)^lp * GSL.sf_coupling_3j( 2*lp, 2*1, 2*l, 0, 0, 0 )
    
    #----------------------------------------------------------------------------------------
    
    return( fac * integral )
end
    

"""
`StrongField.scalarProdBoundContHydrogenic(Pepsplp::Array{ComplexF64,1}, epsiloni::Float64, n::Int, l::Int, rGrid::Radial.Grid)`  
    ... computes the scalar product of the bound and continuum states in the one-particle picture for hydrogenic initial states
        where the radial wave function Pepsplp of |epsilonp lp jp> needs to be provided as an argument on the grid rGrid
"""
function  scalarProdBoundContHydrogenic(Pepsplp::Array{ComplexF64,1}, epsiloni::Float64, n::Int, l::Int, rGrid::Radial.Grid)
    
    if rGrid.meshType != Radial.MeshGL()
        println("StrongField module needs a radial grid of meshType Radial.MeshGL() to perform the radial integrals.")
        return(0.)
    end
    
    rgrid = rGrid.r
    orderGL = size(rgrid)[1]
    weights = rGrid.wr
    
    #Prepare radial wave functions
    Pnl = HydrogenPnl( epsiloni, n, l, rgrid )
    DPnl = HydrogenDPnlDr( epsiloni, n, l, rgrid )       
    
    integral = 0. * im
    
    #Sum over grid and compute Gauss-Legendre sum
    for    k = 1:orderGL
        r = rgrid[k]
        integrand = conj( Pepsplp[k] ) * Pnl[k]

        #Gauss-Legendre sum
        integral = integral + weights[k] * integrand
    end
    
    return((-im)^l * integral)
end

