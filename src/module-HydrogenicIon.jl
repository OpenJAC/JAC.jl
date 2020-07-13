
"""
`module  JAC.HydrogenicIon`  
    ... a submodel of JAC that contains methods for computing various one-electron energies, matrix elements, etc..
"""
module HydrogenicIon 

    using Printf, ..Basics, ..Defaults, ..Math, ..Nuclear, ..Radial, ..RadialIntegrals,  GSL


    """
    `HydrogenicIon.energy(sh::Shell, Z::Float64)` 
        ... to compute the (non-relativistic) energy for the given shell and for a point-like nuclear charge Z; the energy is 
            printed in the current energy units to screen but is returned in Hartree. A energy::Float64 is returned.
    """
    function energy(sh::Shell, Z::Float64)
        energy = - Z^2 / sh.n^2 / 2.
        energx = Defaults.convertUnits("energy: from atomic to eV", energy)
        sa     = "  Energy for shell $sh is [in $(Defaults.getDefaults("unit: energy"))]: " * @sprintf("%.8e", energx)
        println(sa)
        return( energy )
    end


    """
    `HydrogenicIon.energy(sh::Subshell, Z::Float64)` 
        ... to computes the Dirac energy for the hydrogenic subshell sh and for point-like nucleus with nuclear charge Z; 
            a energy::Float64 in atomic units and without the rest energy of the electron is returned. That is the binding 
            energy of a 1s_1/2 electron for Z=1 is -0.50000665.
    """
    function energy(sh::Subshell, Z::Float64)
        if  Z <= 0.1    error("Requires nuclear charge Z >= 0.1")    end
        # Compute the energy from the Dirac formula
        jPlusHalf = (Basics.subshell_2j(sh) + 1) / 2;   nr = sh.n - jPlusHalf;    alpha = Defaults.getDefaults("alpha") 
        wa = sqrt(jPlusHalf^2 - Z^2 * alpha^2 )  
        wa = sqrt(1.0 + Z^2 * alpha^2 / (nr + wa)^2)
        wa = Defaults.getDefaults("speed of light: c")^2 * (1/wa - 1.0)
        wb = Defaults.convertUnits("energy: from atomic to eV", wa)
        sa = "  Energy for subshell $sh is [in $(Defaults.getDefaults("unit: energy"))]: " * @sprintf("%.8e", wb)
        println(sa)
        return( wa )
    end


    """
    `HydrogenicIon.radialOrbital(sh::Shell, Z::Float64, r::Float64)`
        ... to compute the (non-relativistic) orbital function P(r) for the given shell and nuclear charge Z; 
            a value::Float64 is returned.
    """
    function radialOrbital(sh::Shell, Z::Float64, r::Float64)
        n = sh.n;    l = sh.l
        value = r^(l+1) / factorial(2l+1) * sqrt( factorial(n+l) / (factorial(n-l-1) *2n)) *
                (2Z / n)^(l+1.5) * exp(-Z*r/n) 
        
        if  abs(value) < 1.e-15  &&  2Z*r/n > 99.    return(0.)   end
        
        value = value * GSL.hypergeom(-(n-l-1.), 2l+2., 2Z*r/n)
        return( value )
    end


    """
    `HydrogenicIon.radialOrbital(sh::Shell, Z::Float64, grid::Radial.Grid)`
        ... to compute the same but for all r-values as specified by the given grid; a PList::Array{Float64,1} is returned.
    """
    function radialOrbital(sh::Shell, Z::Float64, grid::Radial.Grid)
        plist = zeros( grid.NoPoints )
        for i = 1:length(plist)    plist[i] = HydrogenicIon.radialOrbital(sh, Z, grid.r[i])    end  
        return( plist )
    end


    """
    `HydrogenicIon.radialOrbital(sh::Shell, Z::Float64, rlist::Array{Float64,1})`
        ... to compute the same but for an array of r-values [r_1, r_2, ...]; a PList::Array{Float64,1} is returned.
    """
    function radialOrbital(sh::Shell, Z::Float64, rlist::Array{Float64,1})
        plist = similar(rlist)
        for i = 1:length(plist)    plist[i] = HydrogenicIon.radialOrbital(sh, Z, rlist[i])    end
        return( plist )
    end


    """
    `HydrogenicIon.radialOrbital(sh::Subshell, Z::Float64, grid::Radial.Grid)`
        ... to compute a relativstic hydrogenic Dirac orbital on the given grid by applying the kinetic-balance to a 
            corresponding non-relavistic orbital; an orbital::Radial.Orbital is returned.
    """
    function radialOrbital(sh::Subshell, Z::Float64, grid::Radial.Grid)
        en  = HydrogenicIon.energy(sh, Z)
        P   = HydrogenicIon.radialOrbital( Shell(sh.n, Basics.subshell_l(sh)), Z, grid)
        
        Q   = zeros(size(P, 1));   Pprime   = zeros(size(P, 1));    Qprime   = zeros(size(P, 1))
        if  grid.meshType == Radial.MeshGrasp()
            dP(i) = Math.derivative(P, i)
            for i = 2:size(Q, 1)
                Q[i] = -1/(2 * Defaults.INVERSE_FINE_STRUCTURE_CONSTANT) * (dP(i) / grid.h / grid.rp[i] + sh.kappa/grid.r[i]) * P[i]
                @warn("radialOrbital():: P' and Q' not yet defined.")
            end
        elseif  grid.meshType == Radial.MeshGL()
            @warn("radialOrbital():: Q[:] = zero everywhere; kinetic-balance not yet defined for Gauss-Legendre grids.")
        else
            error("stop a")
        end
    
        orb   = Radial.Orbital( sh, true, false, en, P, Q, Pprime, Qprime, grid)
        norma = RadialIntegrals.overlap(orb, orb, grid)
        orb.P = orb.P/sqrt(norma)
        orb.Q = orb.Q/sqrt(norma)
        normb = RadialIntegrals.overlap(orb, orb, grid)
        println("HydrogenicIon.radialOrbital():  for subshell $sh : norm-before = $norma, norm-after = $normb")
        
        return( orb )
    end


    """
    `HydrogenicIon.radialOrbital(sh::Subshell, nm::Nuclear.Model, grid::Radial.Grid)`
        ... to compute a relativstic hydrogenic Dirac orbital for the given nuclear model by using an explicit diagonalization
            of the Dirac Hamiltonian in a B-spline basis; an orbital::Radial.Orbital is returned.
    """
    function radialOrbital(sh::Subshell, nm::Nuclear.Model, grid::Radial.Grid)
        error("HydrogenicIon.radialOrbital() Not yet implemented.")
        return( orb )
    end
    
end # module
