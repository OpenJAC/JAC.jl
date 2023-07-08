  
    using  QuadGK, ..Basics,  ..Defaults, ..Nuclear, ..Radial, ..Math
    export Model
  

    """
    `Plasma.computeElectronNumberDensity(grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64)`  
        ... computes the electron number density ne from the one-particle energies of the electrons, the
            chemical potential chemMu [a.u.] and the temperature [a.u.]. 
    """
    function computeElectronNumberDensity(grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64)
        ne = zeros(length(grid.r))
        for (k,v) in orbitals
            wa = Basics.FermiDirac(v.energy, chemMu, temp)
            for (ir, r)  in  enumerate(grid.r)     
                ne[ir] = ne[ir] + wa * (v.P[ir]^2 + v.Q[ir]^2)
            end
        end
        
        return ( ne )
    end
 

    """
    `Plasma.determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)`  
        ... determines the chemical potential so that Sum_i f(epsilon_i, mu, temp) = Z. 
            The Newton-Raphson methods is used to iterate to the chemical potential; a chemMu::Float64 is returned.
    """
    function determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
        function g(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
            wa = - nm.Z
            for  (k,v)  in orbitals
                occ = Basics.twice( Basics.subshell_j(k)) + 1
                wb  = (v.energy - mu) /temp
                if  wb > 300.   wb = 300.   end
                ## @show occ, wb, wa
                wa  = wa + occ / (exp(wb) + 1)    
            end
            return( wa )
        end
        ## for  mu = -18.1:0.005:-18.0  @show  mu, g(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)  end
        ## error("stop x")
        #
        function gprime(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
            wa = 0.
            for  (k,v)  in orbitals
                occ = Basics.twice( Basics.subshell_j(k)) + 1
                wb  = (v.energy - mu) /temp
                if  wb > 300.   wb = 300.   end
                wc  = exp( wb )
                wa  = wa + occ * wc^2 / temp / (wc+1)^2    
                ## @show occ, wb, wc, wa
            end
            return( wa )
        end
        # Iterate for the chemical potential
        chemMu = -0.1;     newMu = 0.;     nx = 0
        while true
            nx = nx + 1
            newMu = chemMu - g(chemMu, orbitals, temp, nm) / gprime(chemMu, orbitals, temp, nm)
            if  abs(newMu - chemMu) < 1.0e-6  break
            else  ## println(">>> Newton-Raphson: $nx)  chemMu = $chemMu  g = $(g(chemMu, orbitals, temp, nm)) ")
                  chemMu = newMu
            end
        end
        
        println(">>> Newton-Raphson: $nx)  chemMu = $chemMu  g = $(g(chemMu, orbitals, temp, nm)) ")
        return ( chemMu )
    end
  

    """
    `Plasma.determineWignerSeitzRadius(rho::Float64, nm::Nuclear.Model)`  
        ... determines the Wigner-Seitz radius R^(WS) from the plasma density rho [g/cm^3] ne and the nuclear charge Z.
    """
    function determineWignerSeitzRadius(rho::Float64, nm::Nuclear.Model)
        
        # Convert the density into rho [A_Z * u/a_o^3]
        wa = Defaults.convertUnits("density: from [g/cm^3] to atomic", rho) / nm.mass
        wb = 4pi * wa
        wr = (3.0 / wb)^(1/3)
        
        println(">>> Wigner-Seitz radius R^(WS) = $wr  [a_o] for density $rho [g/cm^3] and atomic mass $(nm.mass).")
        
        return ( wr )
    end

        
    """
    `Plasma.perform(scheme::Plasma.AverageAtomScheme, computation::Plasma.Computation; output::Bool=tru)`  
        ... to perform an average-atom plasma computation for an atoms with nuclear charge Z which generates a self-consistent set 
            of orbitals. For output=true, dictionary is returned from which the relevant results can be can easily accessed by 
            proper keys.
    """
    function  perform(scheme::Plasma.AverageAtomScheme, computation::Plasma.Computation; output::Bool=true)
        if  output    results = Dict{String, Any}()    else    results = nothing    end
            
        nm   = computation.nuclearModel
        wa   = Bsplines.generatePrimitives(computation.grid)
        temp = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", computation.settings.temperature)

        # Generate a subshell list that is taken into account into the average atom computations
        shells    = Basics.generateShellList(1, scheme.nMax, scheme.lMax)
        subshells = Basics.generateSubshellList(shells)
        println(">> Subshells included into the average-atom scheme: \n   $subshells")
        # Generate hydrogenic orbitals for all subshells and display comparison
        orbitals  = Bsplines.generateOrbitalsHydrogenic(wa, nm, subshells; printout=false)
        Basics.displayOrbitalProperties(stdout, orbitals, -nm.Z, temp, nm, computation.grid)
        chemMu    = Plasma.determineChemicalPotential(orbitals, temp, nm)
        # Solve the orbitals and chemical potential self-consistently in the average-atom model
        orbitals  = Bsplines.solveSelfConsistentAverageAtom(wa, nm, orbitals, temp, scheme.scField, printout=true)
        chemMu    = Plasma.determineChemicalPotential(orbitals, temp, nm)
    
        
        Defaults.warn(PrintWarnings())
        Defaults.warn(ResetWarnings())
        return( results )
    end
    
    

