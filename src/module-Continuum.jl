
"""
`module  JAC.Continuum`  
    ... a submodel of JAC that contains all functions and methods for computing, either in a given potential or with
        regard to a well-defined bound state.
"""
module Continuum

    using  GSL, Printf, SpecialFunctions
    using  ..Basics, ..Bsplines, ..Defaults, ..ManyElectron, ..Radial, ..Nuclear
    ##x global JAC_counter = 0


    """
    `struct  Continuum.Settings`  ... defines a type for the parameters for computing continuum orbitals.

        + includeExchange         ::Bool    ... True, if the exchange is to be included and false otherwise.
        + mtp                     ::Int64   ... No of grid points for which the continuum orbital(s) are to be computed.
    """
    struct Settings 
        includeExchange           ::Bool  
        mtp                       ::Int64  
    end 


    """
    `Continuum.generateOrbitalAsymptoticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`   
        ... to generate a simple relativistic and non-normalized continuum orbital but by assuming an asymptotic sin/cos behaviour 
            for both components at all mesh points. For this behaviour, both components and their derivatives can be calculated 
            analytically. The effective charge Zbar is determined for the given potential and the non-Coulombic phase shifts is set 
            to zero. An non-normalized cOrbital::Orbital is returned.
         
            Warning: At present, only the real part eta.re is taken into account; check for consistency !!
    """
    function generateOrbitalAsymptoticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)  
        P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
        
        phi0 = 0.3
        wc   = Defaults.getDefaults("speed of light: c");    q = sqrt( 2*energy + energy/(wc*wc));    kappa = sh.kappa
        Zbar = Radial.determineZbar(pot)
        y    = Zbar * (energy + wc^2) / (wc^2 * q);   gammaBar = sqrt(kappa^2 - Zbar^2 / wc^2)
        eta  = - (kappa - im*y) / (energy + wc^2) / (gammaBar + im*y) / (2im)
        NP   = sqrt( (energy + 2*wc^2)/(pi*wc^2*q) );   NQ = -sqrt( energy / (pi*wc^2*q) )
        
        ## println("Continuum-asymptotic Coulomb orbital for q = $q,  kappa = $kappa,  Zbar=$Zbar  y = $y,   gammaBar = $gammaBar,  eta = $eta")
        if  Defaults.GBL_PRINT_DEBUG  println("Continuum-asymptotic Coulomb orbital for q = $q,  kappa = $kappa,  Zbar=$Zbar,  eta = $eta")    end

        
        for  i = 2:settings.mtp  
            theta = q * pot.grid.r[i]  +  y * log(2q * pot.grid.r[i])  - angle( SpecialFunctions.gamma(gammaBar + im*y) )  -  pi*gammaBar/2  +  eta.re
            P[i] = NP * cos(theta + phi0);      Pprime[i] = - NP * (q + y/pot.grid.r[i]) * sin(theta + phi0)       
            Q[i] = - NQ * sin(theta + phi0);    Qprime[i] = - NP * (q + y/pot.grid.r[i]) * cos(theta + phi0)       
        end

        cOrbital = Orbital( sh, false, true, energy, P, Q, Pprime, Qprime, Radial.Grid())
        
        return( cOrbital )
    end


    """
    `Continuum.generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid, basis::Basis,
                                       settings::Continuum.Settings)`   
        ... to generate a continuum orbital for the (continuum) subshell sh, the energy and the effective charge within the given 
            potential. The continuum orbital is generated orthogonal with regard to all subshells of the same symmetry in the basis. 
            All further specifications about this generations are made by proper settings. A tupel of a (continuum) 
            (orbital::Orbital, phase::Float64) is returned.
    """
    function generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid, settings::Continuum.Settings)  
        #
        # Generate a (local) potential for the given level
        nuclearPotential  = Nuclear.nuclearPotential(nm, grid)
        ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
        ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
        ## wp3 = compute("radial potential: Kohn-Sham", grid, wLevel)
        ## wp4 = compute("radial potential: Dirac-Fock-Slater", grid, wLevel)
        wp = compute("radial potential: Kohn-Sham", grid, level)   
        pot = Basics.add(nuclearPotential, wp)
        Defaults.warn(AddWarning, "All continuum orbitals are generated in a local (DFS) potential.")  
        
        # Generate a continuum orbital due to the given solution method
        ##x println("Defaults.GBL_CONT_SOLUTION = $(Defaults.GBL_CONT_SOLUTION)")
        if      Defaults.GBL_CONT_SOLUTION  ==  ContBessel 
            cOrbital = Continuum.generateOrbitalBessel(energy, sh, grid::Radial.Grid, settings)
        elseif  Defaults.GBL_CONT_SOLUTION  ==  AsymptoticCoulomb 
            cOrbital = Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)
        elseif  Defaults.GBL_CONT_SOLUTION  ==  NonrelativisticCoulomb 
            cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, pot.grid, settings)
        elseif  Defaults.GBL_CONT_SOLUTION  ==  BsplineGalerkin
            cOrbital = Continuum.generateOrbitalGalerkin(energy, sh, pot, settings)
        else    error("stop a")
        end
        #
        #
        # Normalize the continuum orbital and determine its phase
        if      Defaults.GBL_CONT_NORMALIZATION  ==  PureSine
            cOrbital, phase = Continuum.normalizeOrbitalPureSine(cOrbital, pot.grid, settings)
        elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombSine
            cOrbital, phase = Continuum.normalizeOrbitalCoulombSine(cOrbital, pot, settings)
        elseif  Defaults.GBL_CONT_NORMALIZATION  ==  OngRussek
            cOrbital, phase = Continuum.normalizeOrbitalOngRussek(cOrbital, pot, settings)
        else    error("stop b")
        end
        
        return( cOrbital, phase )
    end


    """
    `Continuum.generateOrbitalLocalPotential(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`   
        ... to generate a continuum orbital for the (continuum) subshell sh, the energy and the given potential. The effective 
            charge for the normalization of the continuum orbital is derived from the potential. All further specifications about 
            this generations are made by proper settings; however, the function termintates if the settings.includeExchange = true. 
            A tupel of a (continuum) (orbital::Orbital, phase::Float64, normFactor::Float64) is returned.
    """
    function generateOrbitalLocalPotential(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)
        settings.includeExchange   &&   error("Continuum orbital for local potential does not allow 'exchange'.")
        
        # Generate a continuum orbital due to the given solution method
        if      Defaults.GBL_CONT_SOLUTION  ==  ContBessel 
            cOrbital = Continuum.generateOrbitalBessel(energy, sh, pot.grid, settings)
        elseif  Defaults.GBL_CONT_SOLUTION  ==  ContSine 
            cOrbital = Continuum.generateOrbitalPureSine(energy, sh, pot.grid, settings)
        elseif  Defaults.GBL_CONT_SOLUTION  ==  AsymptoticCoulomb 
            cOrbital = Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)
        elseif  Defaults.GBL_CONT_SOLUTION  ==  NonrelativisticCoulomb 
            cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, pot.grid, settings)
        elseif  Defaults.GBL_CONT_SOLUTION  ==  BsplineGalerkin
            cOrbital = Continuum.generateOrbitalGalerkin(energy, sh, pot, settings)
        else    error("stop a")
        end
        #
        ## return( cOrbital, 0., 0. )
        #
        # Normalize the continuum orbital and determine its phase
        if      Defaults.GBL_CONT_NORMALIZATION  ==  PureSine
            cOrbital, phase, normFactor = Continuum.normalizeOrbitalPureSine(cOrbital, pot.grid, settings)
        elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombSine
            cOrbital, phase, normFactor = Continuum.normalizeOrbitalCoulombSine(cOrbital, pot, settings)
        elseif  Defaults.GBL_CONT_NORMALIZATION  ==  OngRussek
            cOrbital, phase, normFactor = Continuum.normalizeOrbitalOngRussek(cOrbital, pot, settings)
        else    error("stop b")
        end
        #
        ##x println("*** phase = $phase   normFactor = $normFactor")
        return( cOrbital, phase, normFactor )
    end


    """
    `Continuum.generateOrbitalBessel(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)`   
        ... to generate a simple, non-relativistic and non-normalized Bessel wave for the large component of the continuum
            orbital and a small component due to the kinetic-balance condition. A (non-normalized) orbital::Orbital is returned.
            While Pprime is obtained from the analytical expression of the large component, Qprime is set to zero here.
    """
    function generateOrbitalBessel(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)  
        P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
        
        q = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c")
        for  i = 2:settings.mtp   P[i] = GSL.sf_bessel_jl(l, q * grid.r[i]) * grid.r[i]    end
        for  i = 2:settings.mtp   qz = q * grid.r[i];   Pprime[i] = (l/qz * GSL.sf_bessel_jl(l, qz) - GSL.sf_bessel_jl(l+1, qz)) + P[i]/grid.r[i]     end
        ## for  i = 1:settings.mtp   P[i]      = sin(q * grid.r[i] + 5.1)    end
        ## for  i = 1:settings.mtp   Pprime[i] = q * cos(q * grid.r[i] + 5.1)    end
        for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]         end
        
        println("Continuum-Bessel orbital for q = $q,  l = $l")

        cOrbital = Orbital( sh, false, false, energy, P, Q, Pprime, Qprime, grid)
        
        return( cOrbital )
    end


    """
    `Continuum.generateOrbitalGalerkin(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`   
        ... to generate a non-normalized continuum orbital within the given local potential by using the Galerkin method and a 
            given B-spline basis. A (non-normalized) orbital::Orbital is returned.
    """
    function generateOrbitalGalerkin(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)  
        P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
        
        wa = Bsplines.generatePrimitives(pot.grid)
        wb = Bsplines.generateGalerkinMatrix(wa, energy, sh, pot)
        wc = wb * adjoint(wb)
        wd  = Basics.diagonalize("matrix: Julia, eigfact", wc)
        wx  =  adjoint(wd.vectors[1]) * wd.vectors[1]
        ##x println("*** wx-norm = $wx")
        ##x println("Galerkin: A-eigenvalues = $(wd.values)")
        
        cOrbital = Bsplines.generateOrbitalFromPrimitives(energy, sh, settings.mtp, wd.vectors[1], wa)  
        mtp = size(cOrbital.P,1)
        println(">> Continuum B-spline-Galerkin orbital for energy=" * @sprintf("%.4e",energy) * ",  kappa=$(sh.kappa) " *
                "[mpt=$mtp, r[mtp]=" * @sprintf("%.4e",pot.grid.r[mtp]) * ", smallest eigenvalue=" * @sprintf("%.4e",wd.values[1]) * "].")
        
        return( cOrbital )
    end

    """
    `Continuum.generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, Zeff::Float64, grid::Radial.Grid, 
                                                     settings::Continuum.Settings)`   
        ... to generate a simple, non-relativistic and non-normalized Coulomb wave for the (continuum) subshell sh, 
            the small component is simply obtained by the kinetic-balance condition. A orbital::Orbital is returned.
            *****
            This function does not yet work because there is no Julia implementation of the hypergeometric function with 
            complex arguments available.
    """
    function generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, Zeff::Float64, grid::Radial.Grid, settings::Continuum.Settings) 
        println("\n error: presently, NO Julia implementation of the hypergeometric function with complex arguments available.")
        
        # Define the derivative of the confluent hypergeometric function
        function hgPrime(a,b,z)    wa = a/z * (GSL.hypergeom(a+1,b,z) - GSL.hypergeom(a,b,z));   return(wa)   end
        
        P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
        
        q  = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c");   np = Zeff/q;   ws = 1
        wa = 2 * sqrt(Zeff) / sqrt(1 - Base.MathConstants.e^(-2pi*np)) / factorial(2l + 1)
        for  s = 1:l    ws = ws * sqrt(s^2 + np^2)    end
        
        for  i = 2:settings.mtp   P[i] = wa * ws * (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * 
                                         GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i])                     end
        for  i = 2:settings.mtp   
            Pprime[i] = 2q*l * (2q*grid.r[i])^(l-1) * Base.MathConstants.e^(-im*q*grid.r[i]) * GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i]) -
                        im*q * (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i])     +
                        (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * hgPrime(l+1+im*np, 2l+2., 2im*q*grid.r[i])
            Pprime[i] = wa * ws * Pprime[i]          
        end
        for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]                end

        cOrbital = Orbital( sh, false, true, energy, P, Q, Pprime, Qprime, Radial.Grid())
        
        return( cOrbital )
    end


    """
    + `(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`   
        ... to generate a simple, non-relativistic and non-normalized Coulomb wave for the (continuum) subshell sh, 
            the small component is simply obtained by the kinetic-balance condition. A orbital::Orbital is returned.
    """
    function generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings) 
        Zbar     = Radial.determineZbar(pot)
        cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, Zbar, pot.grid, settings)
        
        return( cOrbital )
    end


    """
    `Continuum.gridConsistency(maxEnergy::Float64, grid::Radial.Grid)`   
        ... to check the consistency of the given grid with the maximum energy of the required continuum electrons; 
            an error message is issued if the grid.hp = 0.   or  15 * grid.hp < wavelength(maxEnergy)   or  if the grid has 
            less than 600 grid points. The function also returns the recommended grid point where the normalization and phase 
            is to be determined. This number is currently set to nrContinuum = grid.nr - 200  ... to correct for the wrong 
            'phase behaviour' at large r-values.
    """
    function gridConsistency(maxEnergy::Float64, grid::Radial.Grid)
        Defaults.warn(AddWarning, "Continuum.gridConsistency(): Improve the grid point for normalization; currently 600.")
        
        wavenb      = sqrt( 2maxEnergy + maxEnergy * Defaults.getDefaults("alpha")^2 )
        wavelgth    = 2pi / wavenb
        nrContinuum = grid.nr - 200
        
        if      grid.hp == 0.               error("Improper grid for continuum processes with grid.hp = 0.")
        elseif  15 * grid.hp > wavelgth     error("Improper grid for continuum processes with 15*grid.hp = $(15*grid.hp) > " *
                                                  "wavelength = $wavelgth .")
        elseif  grid.nr < 600               error("Improper No of grid points for continuum processes; grid.nr = $(grid.nr).")
        elseif  grid.r[nrContinuum] < 2.0   error("Improper grid extent for continuum processes; grid.r[nrContinuum] = " *
                                                  "$(grid.r[nrContinuum]) < 2.")
        end
        
        return( nrContinuum )
    end


    """
    `Continuum.generateOrbitalPureSine(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)`   
        ... to generate a simple, non-relativistic and non-normalized Bessel wave for the large component of the continuum
            orbital and a small component due to the kinetic-balance condition. A (non-normalized) orbital::Orbital is returned.
            While Pprime is obtained from the analytical expression of the large component, Qprime is set to zero here.
    """
    function generateOrbitalPureSine(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)  
        P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
        
        phi0 = 0.3
        q = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c");    N = sqrt(2/(pi*q))
        for  i = 1:settings.mtp   P[i]      = N * sin(q * grid.r[i] - l*pi/2 + phi0)                  end
        for  i = 1:settings.mtp   Pprime[i] = N * q * cos(q * grid.r[i] - l*pi/2 + phi0)              end
        for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]        end
        
        println(">> Continuum-pure sine orbital for q = $q,  l = $l")

        cOrbital = Orbital( sh, false, false, energy, P, Q, Pprime, Qprime, grid)
        
        return( cOrbital )
    end


    """
    `Continuum.normalizeOrbitalPureSine(cOrbital::Orbital, grid::Radial.Grid, settings::Continuum.Settings)`   
        ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function. 
            An (on-energy-scale-normalized) orbital::Orbital and its relative phase phi w.r.t.  sin(kr + phi) is returned.
    """
    function normalizeOrbitalPureSine(cOrbital::Orbital, grid::Radial.Grid, settings::Continuum.Settings) 
        mtp = size( cOrbital.P, 1) - 20;    q = sqrt( 2*cOrbital.energy );        l = Basics.subshell_l(cOrbital.subshell)

        if  true
            i = mtp;    PPprime = cOrbital.P[i] / cOrbital.Pprime[i];    kr  = q * grid.r[i];   at = atan( PPprime * q )
            phi = atan( PPprime * q ) - kr + l*pi/2;         phi = rem(phi, 2*pi) + pi   ##x to bring phi in the intervaö 0 <= phi < pi
            A   = cOrbital.P[i] / sin(kr - l*pi/2 + phi);    N   = sqrt(2/(pi*q)) / A
            mPhi = phi;                                      mN  = N
            println(">> Pure-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",mN) *
                    " and phase phi=" * @sprintf("%.4e",mPhi) *
                    " at r=" * @sprintf("%.4e",grid.r[mtp])   * " a.u." )
        else
            # Compute a mean phase, normalization factor and their mean deviations
            nx = 1   # Number of grid points for determining the phase and normalization constants
            meanPhi = zeros(nx);    devsPhi = zeros(nx);   meanN = zeros(nx);    devsN = zeros(nx);   ny = 0
            for i = mtp-nx+1:mtp
                PPprime = cOrbital.P[i] / cOrbital.Pprime[i];    kr  = q * grid.r[i];   at = atan( PPprime * q )
                phi = atan( PPprime * q ) - kr + l*pi/2;         phi = rem(phi, 2*pi) + pi   ##x to bring phi in the intervaö 0 <= phi < pi
                A   = cOrbital.P[i] / sin(kr - l*pi/2 + phi);    N   = sqrt(2/(pi*q)) / A
                if  abs(N) > 1.0e1
                    println("** Skip normalization at i = $i, r[i] = $(grid.r[i])  A = $A   phi = $phi  corbital = $(cOrbital.P[i])");   continue   end
                # Collect data to form the mean and standard deviations
                ny  = ny + 1;    meanPhi[ny] = phi;     meanN[ny] = N
            end
            #
            # Determine the mean and standard deviations
            mPhi = sum(meanPhi) / nx;    mN = sum(meanN) / nx
            for  i = 1:nx    devsPhi[i] = (meanPhi[i] - mPhi)^2;    devsN[i] = (meanN[i] - mN)^2    end
            stdPhi = sqrt( sum(devsPhi) / nx );   stdN = sqrt( sum(devsN) / nx );
            println(">> Pure-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",mN) *
                    " (Delta-N=" * @sprintf("%.4e",stdN) *
                    ") and phase phi=" * @sprintf("%.4e",mPhi) *
                    " (Delta-phi=" * @sprintf("%.4e",stdPhi) * 
                    ") at r=" * @sprintf("%.4e",grid.r[mtp]) * " a.u." )
        end
        
        P = mN .* cOrbital.P;   Q = mN .* cOrbital.Q;   Pprime = mN .* cOrbital.Pprime;   Qprime = mN .* cOrbital.Qprime;   
        newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                              P, Q, Pprime, Qprime, cOrbital.grid)
        mN   = maximum( P[end-300:end] )   ## determine the maximum value of the P component as asymptotic measure for the amplitude
                              
        return( newOrbital, mPhi, mN )
    end


    """
    `Continuum.normalizeOrbitalCoulombSine(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)`   
        ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function. 
            An (on-energy-scale-normalized) orbital::Orbital is returned.
    """
    function normalizeOrbitalCoulombSine(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings) 
        mtp = size( cOrbital.P, 1) - 2;   energy = cOrbital.energy    
        wc   = Defaults.getDefaults("speed of light: c");    q = sqrt( 2*energy + energy/(wc*wc));      kappa = cOrbital.subshell.kappa
        Zbar = Radial.determineZbar(pot);    y = Zbar * (energy + wc^2) / (wc^2 * q);   gammaBar = sqrt(kappa^2 - Zbar^2 / wc^2)
        eta  = - (kappa - im*y) / (energy + wc^2) / (gammaBar + im*y) / (2im)
        NP   = sqrt( (energy + 2*wc^2)/(pi*wc^2*q) );   NQ = -sqrt( energy / (pi*wc^2*q) )
        
        if  true
            i = mtp;     thetaprime = q + y/pot.grid.r[i]        
            # Select the 'ratio' on which the thetaprimephase is determined
            if  false   PprimeP = cOrbital.Pprime[i] / cOrbital.P[i];   at = atan( - PprimeP / thetaprime );   son = "P'/P"    
            else        PprimeP = cOrbital.Q[i] / cOrbital.P[i];        at = atan( - PprimeP / NQ * NP );      son = "Q/P"     end
            #
            theta = q * pot.grid.r[i]  +  y * log(2q * pot.grid.r[i])  - angle( SpecialFunctions.gamma(gammaBar + im*y) )  -  pi*gammaBar/2  +  eta.re
            phi = at - theta;       phi = rem(phi, 2*pi) + pi
            A   = cOrbital.P[i] / cos(theta + phi);             N   = NP / A
            mPhi = phi;                                         mN  = N
            println(">> Coulomb-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",mN) *
                    " and phase phi=" * @sprintf("%.4e",mPhi) *
                    " at r=" * @sprintf("%.4e",pot.grid.r[mtp])   * " a.u." )
        else
            # Compute a mean phase, normalization factor and their mean deviations
            nx = 1   # Number of grid points for determining the phase and normalization constants
            meanPhi = zeros(nx);    devsPhi = zeros(nx);   meanN = zeros(nx);    devsN = zeros(nx);   ny = 0;   son = "xxx";   at = 0.
            #
            for i = mtp-nx+1:mtp
                thetaprime = q + y/pot.grid.r[i] 
                # Select the 'ratio' on which the thetaprimephase is determined
                if  false   PprimeP = cOrbital.Pprime[i] / cOrbital.P[i];   at = atan( - PprimeP / thetaprime );   son = "P'/P"    
                else        PprimeP = cOrbital.Q[i] / cOrbital.P[i];        at = atan( - PprimeP / NQ * NP );      son = "Q/P"     end
                #
                theta = q * pot.grid.r[i]  +  y * log(2q * pot.grid.r[i])  - angle( SpecialFunctions.gamma(gammaBar + im*y) )  -  pi*gammaBar/2  +  eta.re
                phi = at - theta;       phi = rem(phi, 2*pi) + pi
                A   = cOrbital.P[i] / cos(theta + phi);             N   = NP / A
                ##x println("***** thetaprime = $thetaprime    theta = $theta    phi = $phi    A = $A")
                if  abs(N) > 1.0e1   
                    println("**Skip normalization at i = $i, r[i] = $(pot.grid.r[i])  A = $A   phi = $phi  corbital $(cOrbital.P[i])");   continue   end
                # Collect data to form the mean and standard deviations
                ny  = ny + 1;    meanPhi[ny] = phi;     meanN[ny] = N
            end
            #
            # Determine the mean and standard deviations
            mPhi = sum(meanPhi) / nx;    mN = sum(meanN) / nx
            for  i = 1:nx    devsPhi[i] = (meanPhi[i] - mPhi)^2;    devsN[i] = (meanN[i] - mN)^2    end
            stdPhi = sqrt( sum(devsPhi) / nx );   stdN = sqrt( sum(devsN) / nx );
            println(">> Coulomb-sine normalized continuum orbital (on " *son* ") with normalization N=" * @sprintf("%.4e",mN) *
                    " (Delta-N=" * @sprintf("%.4e",stdN) *
                    ") and phase phi=" * @sprintf("%.4e",mPhi) *
                    " (Delta-phi=" * @sprintf("%.4e",stdPhi) * 
                    ") at r=" * @sprintf("%.4e",pot.grid.r[mtp]) * " a.u." )
        end
        
        P = mN .* cOrbital.P;   Q = mN .* cOrbital.Q;   Pprime = mN .* cOrbital.Pprime;   Qprime = mN .* cOrbital.Qprime;   
        newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                              P, Q, Pprime, Qprime, cOrbital.grid)
        mN = maximum( P[end-300:end] )   ## determine the maximum value of the P component as asymptotic measure for the amplitude
        
        return( newOrbital, mPhi, mN )
    end


    """
    `Continuum.normalizeOrbitalOngRussek(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)`   
        ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function. 
            An (on-energy-scale-normalized) orbital::Orbital is returned.
    """
    function normalizeOrbitalOngRussek(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings) 
        mtp = size( cOrbital.P, 1) - 30;                    energy = cOrbital.energy;      kappa = cOrbital.subshell.kappa
        wc   = Defaults.getDefaults("speed of light: c");   E = wc^2 + cOrbital.energy;    q = sqrt( 2*energy + energy/(wc*wc))     
        qq = sqrt(E^2/wc^2 - wc^2)
        @show q, qq
        
        if  true
            i = mtp;     dVdr = pot.Zr[i] / (pot.grid.r[i]^2);    V  = - pot.Zr[i] / pot.grid.r[i]
            Pe = cOrbital.P[i];   Pep = cOrbital.Pprime[i];       r0 = pot.grid.r[i]
            # Compute U-function
            ##x @show wc /2.
            ##x @show ((E-V)*dVdr - wc^2 *kappa*(kappa+1) /r0^3)
            ##x @show (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1) / r0^2))
            @show Pe, Pep
            U = wc /2. * ((E-V)*dVdr - wc^2 *kappa*(kappa+1) /r0^3) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1) / r0^2)^(3/2)) * Pe
            ## U = U - wc * (Pep + dVdr * Pe / (2*(E-V+ 2* wc^2)) ) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2))
            U = U - wc * (Pep + dVdr * Pe / (2*(E-V+ 2* wc^2)) ) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2))
            # Compute A-function
            A = (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2)) / (E-V+wc^2) * (Pe^2 + U^2)
            A = sqrt(A);              phir0 = atan(U, Pe)
            phi = phir0 - q*r0;       phi = rem(phi, 2*pi) + pi
            mPhi = phi;               mN  = 1 / (A * sqrt(pi*wc))
            println(">> WKB normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",mN) *
                    " and phase phi=" * @sprintf("%.4e",mPhi) *
                    " at r=" * @sprintf("%.4e",pot.grid.r[mtp])   * " a.u.  U/Pe = $(U/Pe)" )

        else  error("stop a")
        end
        
        P = mN .* cOrbital.P;   Q = mN .* cOrbital.Q;   Pprime = mN .* cOrbital.Pprime;   Qprime = mN .* cOrbital.Qprime;   
        newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                              P, Q, Pprime, Qprime, cOrbital.grid)
        mN = maximum( P[end-300:end] )   ## determine the maximum value of the P component as asymptotic measure for the amplitude

        return( newOrbital, mPhi, mN )
    end

end # module
