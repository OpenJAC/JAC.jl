export  compute

"""
`JAC.compute()`  ... computes various (radial and angular) functions.

  + `("angular coefficients: e-e, Ratip2013", csfa::CsfR, csfb::CsfR)`  ... to compute the angular coefficients in the decomposition of a 
                             (reduced) many-electron matrix element with a general rank-0 electron-electron interaction operator 
                             ``⟨csf_a ||V(e-e)|| csfb⟩ = \\sum_t  T(a_t, b_t, c_t, d_t) * R^k (a_t, b_t, c_t, d_t)``
                             by a call to the Fortran procedure `anco_calculate_csf_pair` of the RATIP program; a 
                             Tuple{Array{AngularTcoeff,1},Array{AngularVcoeff,1}}` is returned.
"""
function compute(sa::String, csfa::CsfR, csfb::CsfR)
    if sa == "angular coefficients: e-e, Ratip2013"
        if csfa == csfb
            JAC.AngularCoefficientsRatip2013.load_csl(csfa)
            t_coeffs, v_coeffs = JAC.AngularCoefficientsRatip2013.angular_coefficients_pair(1,1)
        else
            JAC.AngularCoefficientsRatip2013.load_csl(csfa, csfb)
            t_coeffs, v_coeffs = JAC.AngularCoefficientsRatip2013.angular_coefficients_pair(1,2)
        end
    else 
        error("Unsupported keystring = ", sa)
    end
    return t_coeffs, v_coeffs
end


"""
  + `("angular coefficients: 1-p, Ratip2013", rank::Int64, csfa::CsfR, csfb::CsfR)`  ... to compute the the angular coefficients in the 
                             decomposition of a (reduced) many-electron matrix element with a general single-particle operator of the given rank 
                             ``⟨csf_a ||O^rank|| csfb⟩ = \\sum_t  T(a_t, b_t) * R (a_t, b_t)``  by a call to the Fortran procedure 
                             `anco_calculate_csf_pair_1p` of the RATIP program; an `Array{AngularTcoeff,1}` is returned.
"""
function compute(sa::String, rank, csfa::CsfR, csfb::CsfR)
    if      sa == "angular coefficients: 1-p, Ratip2013"
        if csfa == csfb
            JAC.AngularCoefficientsRatip2013.load_csl(csfa)
            t_coeffs = JAC.AngularCoefficientsRatip2013.angular_coefficients_pair_1p(rank,1,1)
        else
            JAC.AngularCoefficientsRatip2013.load_csl(csfa, csfb)
            t_coeffs = JAC.AngularCoefficientsRatip2013.angular_coefficients_pair_1p(rank,1,2)
        end
    else    error("Unsupported keystring = ", sa)
    end
    return t_coeffs
end


"""
  + `("angular coefficients: 1-p, Grasp92", parity, rank::Integer, csfa::CsfR, csfb::CsfR)`  ... to compute the the angular coefficients
               in the decomposition of a (reduced) many-electron matrix element with a general single-particle operator of the given parity 
               and rank by a call to the Fortran procedure `mct_generate_coefficients` of the RATIP program; an `Array{AngularTcoeff,1}` 
               is returned.
"""
function compute(sa::String, parity, rank::Integer, csfa::CsfR, csfb::CsfR)
    if      sa == "angular coefficients: 1-p, Grasp92"
        if csfa == csfb
            subshells  = JAC.AngularCoefficientsRatip2013.load_csl(csfa)
            mct_coeffs = JAC.AngularCoefficientsRatip2013.mct_generate_coefficients(1, 1, 1, 1, Int32(parity), rank)
        else
            # Add 'zeros' to the fields if the length of occupation does not agree
            if       ( nz = length(csfa.occupation) - length(csfb.occupation) ) <  0   csfa = Basics.addZerosToCsfR( -nz, csfa)
            elseif   ( nz = length(csfa.occupation) - length(csfb.occupation) ) >  0   csfb = Basics.addZerosToCsfR(  nz, csfb)    end
            subshells  = JAC.AngularCoefficientsRatip2013.load_csl(csfa, csfb)
            mct_coeffs = JAC.AngularCoefficientsRatip2013.mct_generate_coefficients(1, 1, 2, 2, Int32(parity), rank)
            ##x println("csfa.subshells = $(csfa.subshells)")
        end
    else    error("Unsupported keystring = ", sa)
    end
    return map(t -> JAC.AngularCoefficientsRatip2013.AngularTcoeff(t, subshells), mct_coeffs) # Convert Fmctcoefficient to AngularTcoeff
end



"""
 + `("matrix: CI, J^P symmetry", JP::LevelSymmetry, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid,
                  settings::AsfSettings; printout::Bool=true)`  ... to compute the CI matrix for a given J^P symmetry block of basis and 
                  by making use of the nuclear model and the grid; a matrix::Array{Float64,2} is returned.
"""
function compute(sa::String, JP::LevelSymmetry, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid,
                             settings::AsfSettings; printout::Bool=true)    
    !(sa == "matrix: CI, J^P symmetry")   &&   error("Not supported keystring")

    # Determine the dimension of the CI matrix and the indices of the CSF with J^P symmetry in the basis
    idx_csf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J == JP.J   &&   basis.csfs[idx].parity == JP.parity    push!(idx_csf, idx)    end
    end
    n = length(idx_csf)
    if printout    print("Compute CI matrix of dimension $n x $n for the symmetry $(string(JP.J))^$(string(JP.parity)) ...")    end

    # Generate an effective nuclear charge Z(r) on the given grid
    potential = JAC.Nuclear.nuclearPotential(nuclearModel, grid)
    if  settings.qedCI    meanPot = compute("radial potential: Dirac-Fock-Slater", grid, basis)   end   

    
    matrix = zeros(Float64, n, n)
    for  r = 1:n
        for  s = 1:n
            wa = compute("angular coefficients: e-e, Ratip2013", basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]])
            me = 0.
            for  coeff in wa[1]
                jj = subshell_2j(basis.orbitals[coeff.a].subshell)
                me = me + coeff.T * sqrt( jj + 1) * JAC.RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
                if  settings.qedCI    
                    me = me + JAC.InteractionStrengthQED.qedLocal(basis.orbitals[coeff.a], basis.orbitals[coeff.b], nuclearModel, meanPot, grid)  
                end
            end

            for  coeff in wa[2]
                if  settings.coulombCI    
                    me = me + coeff.V * JAC.InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                     basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)   end
                if  settings.breitCI
                    me = me + coeff.V * JAC.InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                         basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)               end
            end
            matrix[r,s] = me
        end
    end 
    if printout    println("   ... done.")    end

    return( matrix )
end



"""
 + `("matrix: CI for plasma, J^P symmetry", JP::LevelSymmetry, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid,
                  settings::AsfSettings, plasmaSettings::PlasmaShift.Settings; printout::Bool=true)`  
    ... to compute the CI matrix for a given J^P symmetry block of basis and by making use of the nuclear model and the grid; 
        a matrix::Array{Float64,2} is returned.
"""
function compute(sa::String, JP::LevelSymmetry, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid,
                             settings::AsfSettings, plasmaSettings::PlasmaShift.Settings; printout::Bool=true)    
    !(sa == "matrix: CI for plasma, J^P symmetry")   &&   error("Not supported keystring")

    # Determine the dimension of the CI matrix and the indices of the CSF with J^P symmetry in the basis
    idx_csf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J == JP.J   &&   basis.csfs[idx].parity == JP.parity    push!(idx_csf, idx)    end
    end
    n = length(idx_csf)
    
    if  settings.breitCI  error("No Breit interaction supported for plasma computations; use breitCI=false  in the asfSettings.")   end   
    if  settings.qedCI    error("No QED estimates supported for plasma computations; use qedCI=false  in the asfSettings.")   end   
    
    # Now distinguis the CI matrix for different plasma models
    if  plasmaSettings.plasmaModel == JAC.PlasmaShift.DebyeHueckel
        if printout    print("Compute DebyeHueckel-CI matrix of dimension $n x $n for the symmetry $(string(JP.J))^$(string(JP.parity)) ...")    end

        # Generate an effective nuclear charge Z(r) for a screened Debye-Hueckel potential on the given grid
        potential = JAC.Nuclear.nuclearPotentialDH(nuclearModel, grid, plasmaSettings.lambdaDebye)

        matrix = zeros(Float64, n, n)
        for  r = 1:n
            for  s = 1:n
                wa = compute("angular coefficients: e-e, Ratip2013", basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]])
                me = 0.
                for  coeff in wa[1]
                    jj = subshell_2j(basis.orbitals[coeff.a].subshell)
                    me = me + coeff.T * sqrt( jj + 1) * JAC.RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
                end

                for  coeff in wa[2]
                    if  settings.coulombCI    
                        me = me + coeff.V * JAC.InteractionStrength.XL_Coulomb_DH(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                            basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, 
                                                                                            plasmaSettings.lambdaDebye)   end
                end
                matrix[r,s] = me
            end
        end 
    else
        error("Unsupported plasma model = $(plasmaSettings.plasmaModel)")
    end
    
    if printout    println("   ... done.")    end

    return( matrix )
end



"""
  + `("radial orbital: NR, Bunge (1993)", subshell::Subshell, Z::Int64)`  ... to compute a radial orbital::Orbital for the given subshell
                           and nuclear charge by using the Roothan-Hartree-Fock data by Bunge et al., Atomic Data and Nuclear Data Tables 
                           53 (1993) 113, as obtained for a non-relativistic RHF computation of the neutral atom. These functions are used 
                           a large component, and the small component is obtained from the kinetic balance condition. Radial orbitals
                           can be obtained for the ground-state configuration for all elements with 2 <= Z <= 54. 
  + `("radial orbital: NR, McLean (1981)", subshell::Subshell, Z::Int64)`  ... to compute a radial orbital::Orbital for the given subshell
                           and nuclear charge by using the Roothan-Hartree-Fock data by McLean and McLean, Atomic Data and Nuclear Data Tables 
                           26 (1981) 197., as obtained for a non-relativistic RHF computation of the neutral atom. These functions are used 
                           a large component, and the small component is obtained from the kinetic balance condition. Radial orbitals
                           can be obtained for the ground-state configuration for all elements with 55 <= Z <= 92. 
"""
function compute(sa::String, subshell::Subshell, Z::Int64)
    if      sa == "radial orbital: NR, Bunge (1993)"          wa = JAC.Radial.OrbitalBunge1993(subshell,Z)
    elseif  sa == "radial orbital: NR, McLean (1981)"         wa = JAC.Radial.OrbitalMcLean1981(subshell,Z)
    elseif  sa == "radial orbital: hydrogenic"                wa = JAC.Radial.OrbitalHydrogenic(subshell,Z)
    elseif  sa == "radial orbital: Thomas-Fermi"              wa = JAC.Radial.OrbitalThomasFermi(subshell,Z)
    else    error("Unsupported keystring = $sa ")
    end

    return( wa )
end



"""
  + `("radial potential: core-Hartree", grid::Radial.Grid, level::Level)`  ... to compute a (radial) core-Hartree potential from the given list 
                         of orbitals; cf. JAC.computePotentialCoreHartree. A potential::RadialPotential is returned. 

  + `("radial potential: Kohn-Sham", grid::Radial.Grid, level::Level)`  ... to compute a (radial) Kohn-Sham potential from the given list of 
                         orbitals; cf. JAC.computePotentialKohnSham. A potential::RadialPotential is returned.

  + `("radial potential: Dirac-Fock-Slater", grid::Radial.Grid, level::Level)`  ... to compute a (radial) Dirac-Fock-Slater potential from the
                         given list of orbitals; this potential is rather simple but includes some undesired self-interaction and exhibits an 
                         asymptotic behaviour. Cf. JAC.computePotentialDFS. A potential::RadialPotential is returned.

  + `("radial potential: Dirac-Fock-Slater", grid::Radial.Grid, basis::Basis)`  ... to compute a (radial) Dirac-Fock-Slater potential but for 
                         mean occupation of the given basis. Cf. JAC.computePotentialDFS. A potential::RadialPotential is returned.

  + `("radial potential: extended Hartree", grid::Radial.Grid, level::Level)`  ... to compute a (radial) extended Hartree potential from the 
                         given list of orbitals; it is a local potential that is based on direct and exchange coefficients as obtained from 
                         the configuration-averaged energy. Cf. JAC.computePotentialExtendedHartree. A potential::RadialPotential is returned.
"""
function compute(sa::String, grid::Radial.Grid, level::Level)
    if      sa == "radial potential: core-Hartree"               wa = JAC.computePotentialCoreHartree(grid, level)
    elseif  sa == "radial potential: Hartree"                    wa = JAC.computePotentialHartree(grid, level)
    elseif  sa == "radial potential: Hartree-Slater"             wa = JAC.computePotentialHartreeSlater(grid, level)
    elseif  sa == "radial potential: Kohn-Sham"                  wa = JAC.computePotentialKohnSham(grid, level)
    elseif  sa == "radial potential: Dirac-Fock-Slater"          wa = JAC.computePotentialDFS(grid, level) 
    elseif  sa == "radial potential: extended-Hartree"           wa = JAC.computePotentialExtendedHartree(grid, level)
    else    error("Unsupported keystring = $sa ")
    end

    return( wa )
end

function compute(sa::String, grid::Radial.Grid, basis::Basis)
    if      sa == "radial potential: Dirac-Fock-Slater"          wa = JAC.computePotentialDFS(grid, basis) 
    else    error("Unsupported keystring = $sa ")
    end

    return( wa )
end



"""
`JAC.computeDiracEnergy(sh::Subshell, Z::Float64)`  ... computes the Dirac energy for the hydrogenic subshell sh and for a point-like nucleus
                        with nuclear charge Z; a energy::Float64 in atomic units and without the rest energy of the electron is returned. 
                        That is the binding energy of a 1s_1/2 electron for Z=1 is -0.50000665.
"""
function computeDiracEnergy(sh::Subshell, Z::Float64)
    if  Z <= 0.1    error("Requires nuclear charge Z >= 0.1")    end
    # Compute the energy from the Dirac formula
    jPlusHalf = (JAC.subshell_2j(sh) + 1) / 2;   nr = sh.n - jPlusHalf;    alpha = Constants.give("alpha") 
    wa = sqrt(jPlusHalf^2 - Z^2 * alpha^2 )  
    wa = sqrt(1.0 + Z^2 * alpha^2 / (nr + wa)^2)
    wa = Constants.give("speed of light: c")^2 * (1/wa - 1.0)
    return( wa )
end



"""
`JAC.computeMeanSubshellOccupation()`

    + (sh::Subshell, levels::Array{Level,1})`  ... computes the mean subshell occupation for the subshell sh and for the given levels; 
                                                   a q::Float64 is returned.
"""
function computeMeanSubshellOccupation(sh::Subshell, levels::Array{Level,1})
    q = 0.
    if  length(levels) < 1   error("stop a")    end
    for  level in levels
        subshells = level.basis.subshells;    nsh = 0
        for  i = 1:length(level.basis.subshells)    if  sh == subshells[i]   nsh = i;   break   end    end
        if   nsh == 0   error("Subshell not found in CSF basis.")   end
        #
        for  i = 1:length(level.basis.csfs)   q = q + abs(level.mc[i])^2 * level.basis.csfs[i].occupation[nsh]    end
    end
    return( q/length(levels) )
end



"""
`JAC.computeMeanSubshellOccupation()`

    + (sh::Subshell, basis::Basis)`  ... computes the mean subshell occupation for the subshell sh and for the given CSF in the basis; 
                                         a q::Float64 is returned.
"""
function computeMeanSubshellOccupation(sh::Subshell, basis::Basis)
    q = 0.
    for  csf in basis.csfs
        subshells = basis.subshells;    nsh = 0
        for  i = 1:length(basis.subshells)    if  sh == subshells[i]   nsh = i;   break   end    end
        if   nsh == 0   error("Subshell not found in CSF basis.")   end
        #
        for  i = 1:length(basis.csfs)   q = q + basis.csfs[i].occupation[nsh]    end
    end
    return( q/length(basis.csfs) )
end



"""
`JAC.computePotentialCoreHartree(grid::Radial.Grid, level::Level)`  ... to compute a (radial) core-Hartree potential for the given level; 
                                 a potential::RadialPotential is returned. The core-Hartree potential is defined by

                                 V_CH(r) = int_0^infty dr'  rho_c(r') / r_>                 with r_> = max(r,r')  
                             
                                 and  rho_c(r) =  sum_a (Pa^2(r) + Qa^2(r))                 the charge density of the core-electrons.  
                                 An Radial.Potential with -r * V_CH(r)  is returned to be consistent with an effective charge Z(r).
"""
function computePotentialCoreHartree(grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.nr
    println("coreSubshells = $(basis.coreSubshells);    subshells = $(basis.subshells)")
    rhoc = zeros( npoints );    wx = zeros( npoints );    wb = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.coreSubshells
    ## for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = computeMeanSubshellOccupation(sh, [level])    ## (JAC.subshell_2j(sh) + 1.0)  for closed-shells
        nrho = length(orb.P)
        for    i = 1:nrho   rhoc[i] = rhoc[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhoc[j]/rg   end
        wb[i] = JAC.RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
    end
    # Define the potential with regard to Z(r)
    for i = 1:npoints    wb[i] = - wb[i] * grid.r[i]   end
    #
    wc = JAC.Radial.Potential("CoreHartree", wb, grid)
    return( wc )
end



"""
`JAC.computePotentialHartree(grid::Radial.Grid, level::Level)`  ... to compute a (radial) Hartree potential for the given level; 
                                 a potential::RadialPotential is returned. Here, the Hartree potential is simply defined by

                                 V_H(r) = - Sum_a  q_a^bar  Y0_aa(r)

                                 An Radial.Potential with V_H(r)  is returned that is consistent with an effective charge Z(r).
"""
function computePotentialHartree(grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.nr;    wx = zeros( npoints )
    # Sum_a ...
    for  a in basis.subshells
        occa = computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a];   nrho = length(orba.P);      rhoaa  = zeros(nrho)
        for  i = 1:nrho    rhoaa[i] = orba.P[i]^2 + orba.Q[i]^2    end
        for  i = 1:nrho    wx[i]    = wx[i] - occa * JAC.RadialIntegrals.Yk_ab(0, grid.r[i], rhoaa, nrho, grid)   end
        Yk = JAC.RadialIntegrals.Yk_ab(0, grid.r[nrho], rhoaa, nrho, grid)
        for  i = nrho+1:npoints       wx[i] = wx[i] - occa * Yk    end
    end
    
    wc = JAC.Radial.Potential("Hartree", wx, grid)
    return( wc )
end


"""
`JAC.computePotentialHartreeSlater(grid::Radial.Grid, level::Level)`  ... to compute a (radial) Hartree-Slater potential for the given level; 
                                   a potential::RadialPotential is returned. Here, the Hartree-Slater potential is defined by
                                                     
                                                                            3  [     3          ]^1/3      r
                                 V_HS(r) = = - Sum_a  q_a^bar  Y0_aa(r)  +  -  [ ---------  rho ]       *  -
                                                                            2  [ 4pi^2 r^2      ]          2

                                 with rho(r) = Sum_a  q_a^bar * (Pa^2(r) + Qa^2(r))   ... charge density of all electrons.  
                                 An Radial.Potential with V_HS(r)  is returned that is consistent with an effective charge Z(r).
"""
function computePotentialHartreeSlater(grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.nr;    wx = zeros( npoints );    rho = zeros( npoints )
    # Sum_a ...
    for  a in basis.subshells
        occa = computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a];   nrho = length(orba.P);      rhoaa  = zeros(nrho)
        for  i = 1:nrho    rhoaa[i] = orba.P[i]^2 + orba.Q[i]^2    end
        for  i = 1:nrho    wx[i]    = wx[i]  - occa * JAC.RadialIntegrals.Yk_ab(0, grid.r[i], rhoaa, nrho, grid)   end
        for  i = 1:nrho    rho[i]   = rho[i] + occa * rhoaa[i]     end
        Yk = JAC.RadialIntegrals.Yk_ab(0, grid.r[nrho], rhoaa, nrho, grid)
        for  i = nrho+1:npoints       wx[i] = wx[i] - occa * Yk    end
    end
    ## for  i = 2:npoints     rho[i]   = rho[i] / (4pi * grid.r[i])   end;       rho[1] = 0.
    for  i = 2:npoints     wx[i]    = wx[i] + (3/2) * (3/(4pi^2 * grid.r[i]^2) * rho[i])^(1/3) * grid.r[i] / 2   end
    #
    wc = JAC.Radial.Potential("Hartree-Slater", wx, grid)
    return( wc )
end



"""
`JAC.computePotentialKohnSham(grid::Radial.Grid, level::Level)`  ... to compute a (radial) Kohn-Sham potential for the given level;
                              a potential::RadialPotential is returned. The Kohn-Sham potential is defined by

                                                                               2    [   81                   ]^(1/3)
                              V_KS(r) =  int_0^infty dr'  rho_t(r') / r_>  --  ---  [--------   r * rho_t(r) ]
                                                                               3 r  [ 32 pi^2                ]

                              with r_> = max(r,r')  and  rho_t(r) = sum_a (Pa^2(r) + Qa^2(r))   ... charge density of all electrons.  
                              An Radial.Potential with -r * V_KS(r)  is returned to be consistent with an effective charge Z(r).
"""
function computePotentialKohnSham(grid::Radial.Grid, level::Level)
    basis = level.basis;       npoints = grid.nr
    rhot = zeros( npoints );   wb = zeros( npoints );   wx = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = computeMeanSubshellOccupation(sh, [level])
        nrho = length(orb.P)
        for    i = 1:nrho   rhot[i] = rhot[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals (without alpha)
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
        wb[i] = JAC.RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        wb[i] = wb[i] - 2 / (3grid.r[i]) * (81 /(32 * pi^2) * grid.r[i] * rhot[i])^(1/3)
    end
    # Define the potential with regard to Z(r)
    for i = 2:npoints    wb[i] = - wb[i] * grid.r[i]   end;    wb[1] = 0.
    #
    wc = JAC.Radial.Potential("Kohn-Sham", wb, grid)
    return( wc )
end



"""
`JAC.computePotentialDFS()`

    + `(grid::Radial.Grid, level::Level)`  ... to compute a (radial) Dirac-Fock-Slater potential for the given level;
                         a potential::RadialPotential is returned. The Dirac-Fock-Slater potential is defined by

                                                                          [   3                   ]^(1/3)
                         V_DFS(r) = int_0^infty dr'  rho_t(r') / r_>  --  [-----------   rho_t(r) ]
                                                                          [ 4 pi^2 r^2            ]

                         with r_> = max(r,r')  and  rho_t(r) = sum_a (Pa^2(r) + Qa^2(r))   ... charge density of all electrons.  
                         An Radial.Potential with -r * V_DFS(r)  is returned to be consistent with an effective charge Z(r).
"""
function computePotentialDFS(grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.nr
    ##x println("coreSubshells = $(basis.coreSubshells);    subshells = $(basis.subshells)")
    rhot = zeros( npoints );    wb = zeros( npoints );    wx = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = computeMeanSubshellOccupation(sh, [level])
        nrho = length(orb.P)
        for    i = 1:nrho   rhot[i] = rhot[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
        wb[i] = JAC.RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        wb[i] = wb[i] - (3 /(4*pi^2 * grid.r[i]^2) * rhot[i])^(1/3)
    end
    # Define the potential with regard to Z(r)
    for i = 2:npoints    wb[i] = - wb[i] * grid.r[i]   end;    wb[1] = 0.
    #
    wc = JAC.Radial.Potential("DFS", wb, grid)
    return( wc )
end



"""
`JAC.computePotentialDFS()`

    + `(grid::Radial.Grid, basis::Basis)`  ... to compute the same but for the mean occupation of the orbitals in the given basis.
"""
function computePotentialDFS(grid::Radial.Grid, basis::Basis)
    npoints = grid.nr
    rhot = zeros( npoints );    wb = zeros( npoints );    wx = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = computeMeanSubshellOccupation(sh, basis)
        nrho = length(orb.P)
        for    i = 1:nrho   rhot[i] = rhot[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
        wb[i] = JAC.RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        wb[i] = wb[i] - (3 /(4*pi^2 * grid.r[i]^2) * rhot[i])^(1/3)
    end
    # Define the potential with regard to Z(r)
    for i = 2:npoints    wb[i] = - wb[i] * grid.r[i]   end;    wb[1] = 0.
    #
    wc = JAC.Radial.Potential("DFS for CSF basis", wb, grid)
    return( wc )
end



"""
`JAC.computePotentialExtendedHartree(grid::Radial.Grid, level::Level)`  ... to compute a (radial) extended-Hartree potential for the given level;
                                     a potential::RadialPotential is returned. The extended-Hartreepotential is defined by Gu (2008, Eq.(8))
                                     and is based on the Y^k_ab (r) function and various direct and exchange coefficients as the arise from
                                     average energy of the configuration.
"""
function computePotentialExtendedHartree(grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.nr
    rhoab = zeros( npoints );    wx = zeros( npoints );    wg = zeros( npoints )
    # First term Sum_ab ...
    for  a in basis.subshells
        occa = computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a]
        nrho = length(orba.P);      rhoa  = zeros(nrho)
        for  i = 1:nrho    rhoa[i] = orba.P[i]^2 + orba.Q[i]^2    end
        for  b in basis.subshells
            occb = computeMeanSubshellOccupation(b, [level])
            orbb = basis.orbitals[b]
            nrho = length(orbb.P);   rhobb = zeros(nrho)
            if  a != b  wc = occa*occb   else   wc = occa * (occb -1.0)   end
            if  wc < 0.  error("Improper potential for subshells with mean occupation < 1; stop a")    end
            for  i = 1:nrho    rhobb[i] = orbb.P[i]^2 + orbb.Q[i]^2       end
            for  i = 1:length(orba.P)    wx[i] = wx[i] + wc * JAC.RadialIntegrals.Yk_ab(0, grid.r[i], rhobb, nrho, grid) * rhoa[i]   end
        end
    end
    # Second term Sum_a ...
    for  a in basis.subshells
        occa = computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a]
        nrho = length(orba.P);      rhoa  = zeros(nrho)
        for  i = 1:nrho    rhoa[i] = orba.P[i]^2 + orba.Q[i]^2   end
        wc = occa * (occa - 1.0)
        if  wc < 0.  error("Improper potential for subshells with mean occupation < 1; stop b")    end
        for  k = 1:10
            fk = JAC.AngularMomentum.fk_ab(a,a);    if  fk == 0   cycle    end
            for  i = 1:length(orba.P)    wx[i] = wx[i] + wc * fk * JAC.RadialIntegrals.Yk_ab(k, grid.r[i], rhoa, nrho, grid) * rhoa[i]   end
        end
    end
    # Third term Sum_a != b ...
    for  a in basis.subshells,   b in basis.subshells
        if  a == b   continue   end 
        occa = computeMeanSubshellOccupation(a, [level])
        occb = computeMeanSubshellOccupation(b, [level])
        orba = basis.orbitals[a]
        orbb = basis.orbitals[b]
        nrho = min(length(orba.P), length(orbb.P));      rhoab  = zeros(nrho)
        for  i = 1:nrho    rhoab[i] = orba.P[i]*orbb.P[i] + orba.Q[i]*orbb.Q[i]   end
        wc = occa * occb
        for  k = 0:10
            gk = JAC.AngularMomentum.gk_ab(a,b);    if  gk == 0   continue    end
            for  i = 1:nrho    wx[i] = wx[i] + wc * gk * JAC.RadialIntegrals.Yk_ab(k, grid.r[i], rhoab, nrho, grid) * rhoab[i]   end
        end
    end
    # Weight factor
    for  a in basis.subshells
        occa = computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a]
        nrho = length(orba.P);      rhoa  = zeros(nrho)
        for  i = 1:nrho    rhoa[i] = orba.P[i]^2 + orba.Q[i]^2   end
        wc = occa
        for  i = 1:length(orba.P)    wg[i] = wg[i] + wc * rhoa[i] * grid.r[i]   end
    end
    #
    for  i = 1:npoints   if  wg[i] > 0.    wx[i] = wx[i]/wg[i]   else   wx[i] = 0.   end    end
    return( wx )
end


