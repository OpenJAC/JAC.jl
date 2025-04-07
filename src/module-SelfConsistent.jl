
"""
`module  JAC.SelfConsistent`  
	... a submodel of JAC that contains all structs and methods to generate self-consistent fields of different 
	    kind and complexity.
"""
module SelfConsistent

using  Printf, ..Basics, ..BsplinesN, ..Defaults, ..Hamiltonian, ..ManyElectron, ..Nuclear, ..Radial, JAC


"""
`SelfConsistent.initializeBasis(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, primitives::BsplinesN.Primitives, 
                                settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Initialized a many-electron basis from the given list of configurations, the nuclear model as well as ASF settings.
        It assumes that a proper set of primitives::Primitives has been initialized before. The initial set of orbitals in this
        basis is determined by the settings::AsfSettings.  A basis::Basis is returned.
"""
function initializeBasis(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, primitives::BsplinesN.Primitives, 
                         settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    NoElectrons = configs[1].NoElectrons;   subshells = Subshell[];   coreSubshells = Subshell[];     csfs = CsfR[] 
    orbitals    = Dict{Subshell, Orbital}()
    
    # Perform some simple tests: Number of electrons must be equal in all configurations
    for  conf in configs   if  conf.NoElectrons != NoElectrons    error("stop a")   end     end
    
    # Generate a full set of relativistic CSF from the given configurations and collect the associated level symmetries
    relconfList = ConfigurationR[]
    for  conf in configs
        ##x wa = Basics.generate("configuration list: relativistic", conf)
        wa = Basics.generateConfigurationRs(conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println(">>> include ", relconfList[i])    end   end
    ##x subshellList = Basics.generate("subshells: ordered list for relativistic configurations", relconfList)
    subshells = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshells; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        ##x newCsfs = Basics.generate("CSF list: from single ConfigurationR", relconf, subshells)
        newCsfs = Basics.generateCsfRs(relconf, subshells)
        append!( csfList, newCsfs)
    end
    
    # Select CSF with requested symmetry if needed
    if  length(levelSymmetries) == 0
        csfs = csfList          # Take all relativistic CSF into account
    else
        for  csf in CsfList
            if  LevelSymmetry(csf.J, csf.parity)  in  levelSymmetries   push!(csfs, csf)    end
        end
    end

    # Determine the number of electrons and the list of coreSubshells
    for  (k,sh)  in  enumerate(subshells)
        mocc = Basics.subshell_2j(sh) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc     is_filled = false;           break   end
        end
        if   is_filled    push!( coreSubshells, subshells[k])      else    break   end
    end
        
    # Initialize the orbitals
    if  typeof(settings.startScfFrom) == StartFromHydrogenic
        if  printout   println("> Start SCF process with hydrogenic orbitals.")   end
        # Generate start orbitals for the SCF field by using B-splines
        orbitals  = BsplinesN.generateOrbitalsHydrogenic(subshells, nuclearModel, primitives; printout=printout)
    elseif  typeof(settings.startScfFrom) == StartFromPrevious
        if  printout   println("> Start SCF process from given list of orbitals.energy")    end
        # Taking starting orbitals for the given dictionary; non-relativistic orbitals with a proper nuclear charge
        # are adapted if no orbital is found
        orbitals = Dict{Subshell, Orbital}()
        for  sh in subshells
            if  haskey(settings.startScfFrom.orbitals, sh)  
                orbitals[sh] = settings.startScfFrom.orbitals[sh]
            else
                println("Start orbitals do not contain an Orbital for subshell $sh ")
                orb          = HydrogenicIon.radialOrbital(subsh, nuclearModel.Z, grid)
                orb          = Orbital(orb.subshell, orb.isBound, true, orb.energy, orb.P, orb.Q, orb.Pprime, orb.Qprime, Radial.Grid())
                orbitals[sh] = orb 
            end
        end
    else  error("stop b")
    end
    
    basis = Basis(true, NoElectrons, subshells, csfs, coreSubshells, orbitals)
    return( basis )
end


"""
`SelfConsistent.performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Performs a SCF computation for the given list of configurations, the nuclear model as well as ASF settings.
        If explicit levelSymmetries are given, only these symmetries are considered. Internally, a proper set of primitives::Primitives 
        is initialized and used in the computations. The generated SCF field is controlled by the settings::AsfSettings.  
        A multiplet::Multiplet is returned.
"""
function performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid, 
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)
    
    # Generate primitives and initialize the many-electron basis
    primitives = BsplinesN.generatePrimitives(grid)    
    basis      = SelfConsistent.initializeBasis(configs, nm, primitives, settings; levelSymmetries, printout)
    
    # Solve a self-consistent field for this basis
    @show settings.scField
    if   typeof(settings.scField)  in  [Basics.DFSField, Basics.DFSwCPField, Basics.HSField]
        basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout) 
    elseif   settings.scField in [Basics.NuclearField()]  && settings.startScfFrom == StartFromHydrogenic() 
        # Return the basis as already generated.
    elseif   settings.scField in [Basics.ALField()] 
        multiplet = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
        basis     = multiplet.levels[1].basis
    elseif   typeof(settings.scField) == Basics.EOLField 
        multiplet = SelfConsistent.solveOptimizedLevelField(basis, nm, primitives, settings; printout=printout)
        basis     = multiplet.levels[1].basis
    else  error("stop a")
    end    
    
    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    mp = Hamiltonian.performCI(basis, nm, grid, settings, printout=printout)

    return( mp )
end


"""
`SelfConsistent.performSCF(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Performs a SCF computation for the given list of configurations, the nuclear model as well as ASF settings.
        If explicit levelSymmetries are given, only these symmetries are considered. Internally, a proper set of primitives::Primitives 
        is initialized and used in the computations. The generated SCF field is controlled by the settings::AsfSettings.  
        A multiplet::Multiplet is returned.
"""
function performSCF(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, 
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)
    
    # Generate primitives
    primitives = BsplinesN.generatePrimitives(grid)    
    
    # Solve a self-consistent field for this basis
    if   typeof(settings.scField)  in  [Basics.DFSField, Basics.DFSwCPField, Basics.HSField]
        basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout) 
    elseif   settings.scField in [Basics.NuclearField()]  && settings.startScfFrom == StartFromHydrogenic() 
        # Return the basis as already generated.
    elseif   settings.scField in [Basics.ALField()] 
        multiplet = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
        basis     = multiplet.levels[1].basis
    elseif   typeof(settings.scField) == Basics.EOLField 
        multiplet = SelfConsistent.solveOptimizedLevelField(basis, nm, primitives, settings; printout=printout)
        basis     = multiplet.levels[1].basis
    else  error("stop a")
    end    
    
    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    mp = Hamiltonian.performCI(basis, nm, settings, grid, printout=printout)

    return( mp )
end


"""
`SelfConsistent.solveAverageAtomField(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, scField::Basics.AbstractScField, 
                                      temp::Float64, radiusWS::Float64, primitives::BsplinesN.Primitives; printout::Bool=true)` 
    ... solves the self-consistent field for a given local average-atom potential as specified by scField 
        A (new) set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function solveAverageAtomField(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, scField::Basics.AbstractScField, 
                               temp::Float64, radiusWS::Float64, primitives::BsplinesN.Primitives; printout::Bool=true)
    # Determine the chemical potential
    chemMu    = JAC.Plasma.determineChemicalPotential(orbitals, temp, radiusWS, nuclearModel, primitives.grid);        @show chemMu
    
    # Extract the kappa's from orbitals
    kappas = Int64[];     for (k,v)  in  orbitals     push!(kappas, k.kappa)    end;    kappas = unique(kappas);   @show kappas

    ## Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
    # Define the storage for the calculations of matrices
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    # Set-up the overlap matrix; compute or fetch the diagonal 'overlap' blocks
    nsL = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = BsplinesN.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = BsplinesN.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # Determine the symmetry block of this basis and define storage for the kappa blocks and orbitals from the last iteration
    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end
    # Determine te nuclear potential once at the beginning
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, grid)
            
    # Start the SCF procedure for all symmetries
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false 
        if  NoIteration >  32
                println(">> Maximum number of SCF iterations = 32 is reached at accuracy " * 
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        #
        for kappa in kappas
            # (1) Re-compute the local potential
            if       scField == Basics.AaHSField()     wp = Basics.computePotentialAtomicAverageHS( grid, previousOrbitals, chemMu, temp)
            elseif   scField == Basics.AaDFSField()    wp = Basics.computePotentialAtomicAverageDFS(grid, previousOrbitals, chemMu, temp)
            else     error("stop potential")
            end
            pot = Basics.add(nuclearPotential, wp)
            
            # (2) Set-up the diagonal part of the Hamiltonian matrix
            wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
            # (3) Solve the generalized eigenvalue problem
            wc = Basics.diagonalize("generalized eigenvalues: LinearAlgebra", wa, wb)
            
            # (4) Analyse and print information about the convergence of the symmetry blocks and the occupied orbitals
            wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
            if  wcBlock > 1.0e-6   go_on = true   end     ## accuracyScf
            for  (k,v)  in  orbitals
                if      k.kappa == kappa
                    newOrbital = generateOrbitalFromPrimitives(k, wc, primitives)
                    wcOrbital  = Basics.analyzeConvergence(previousOrbitals[k], newOrbital)
                    if  wcOrbital > 1.0e-6   accuracyScf = wcOrbital;   go_on = true   end     ## accuracyScf
                        sa = "  $k::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "  
                        sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                        sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                        if  printout    println(sa)    end
                    ## println("  $sh  en [a.u.] = $(newOrbital.energy)   self-consistency = $(wcOrbital), $(wcBlock) [kappa=$kappa] ") 
                    previousOrbitals[k] = newOrbital
                end
            end
            # (5) Re-define the bsplineBlock
            bsplineBlock[kappa] = wc
        end
        chemMu              = JAC.Plasma.determineChemicalPotential(previousOrbitals, temp, radiusWS, nuclearModel, primitives.grid)
        if  go_on   nothing   else   break   end
    end
    
    analyzedOrbitals = Basics.analyze(previousOrbitals, printout=true)    
    return( analyzedOrbitals )
end


"""
`SelfConsistent.solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::BsplinesN.Primitives, 
                                       settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        for the average-level (AL) functional. In addition, the settings::AsfSettings are taken into account.
        A (new) basis::Basis is returned.
"""
function solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::BsplinesN.Primitives, 
                                settings::AsfSettings; printout::Bool=true)
    error("Not yet implemented; cf. Bsplines (old)")
    
    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, previousOrbitals)
    return( newBasis )
end


"""
`SelfConsistent.solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::BsplinesN.Primitives, 
                                         settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        for the average-level (AL) functional. In addition, the settings::AsfSettings are taken into account.
        A (new) multiplet::Multiplet is returned.
"""
function solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::BsplinesN.Primitives, 
                                  settings::AsfSettings; printout::Bool=true)
    error("Not yet implemented; this was never done so far.")
    
    multiplet = Multiplet()
    return( multiplet )
end


"""
`SelfConsistent.solveMeanFieldBasis(basis::Basis, nuclearModel::Nuclear.Model, primitives::BsplinesN.Primitives, 
                                    settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        the (local) mean-field potential as specified by the settings::AsfSettings. A (new) basis::Basis is returned.
"""
function solveMeanFieldBasis(basis::Basis, nuclearModel::Nuclear.Model, primitives::BsplinesN.Primitives, 
                             settings::AsfSettings; printout::Bool=true)
    ## Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
    # Define the storage for the calculations of matrices
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    # Set-up the overlap matrix; compute or fetch the diagonal 'overlap' blocks
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = BsplinesN.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = BsplinesN.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # Determine te nuclear potential once at the beginning
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    
    # Determine the symmetry block of this basis and define storage for the kappa blocks and orbitals from the last iteration
    kappas   = Int64[];   for sh in basis.subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)

    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(basis.orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end
    
    # Start the SCF procedure for all symmetries
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false 
        if  NoIteration >  settings.maxIterationsScf
                println(">> Maximum number of SCF iterations = $(settings.maxIterationsScf) is reached at accuracy " * 
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        #
        for kappa in kappas
            # (1) First re-define an (arbitrary) 'level' that represents the mean occupation for the local potential
            wBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, previousOrbitals)
            NoCsf  = length(wBasis.csfs)
            wmc    = zeros( NoCsf );   wN = 0.
            for i = 1:NoCsf   wmc[i] = Basics.twice(wBasis.csfs[i].J) + 1.0;   wN = wN + abs(wmc[i])^2    end
            for i = 1:NoCsf   wmc[i] = wmc[i] / sqrt(wN)   end
            wLevel = Level( AngularJ64(0), AngularM64(0), Basics.plus, 0, -1., 0., true, wBasis, wmc)
            # (2) Re-compute the local potential
            wp  = Basics.computePotential(settings.scField, primitives.grid, wLevel)
            pot = Basics.add(nuclearPotential, wp)
            # (3) Set-up the diagonal part of the Hamiltonian matrix
            wa = BsplinesN.setupLocalMatrix(kappa, primitives, pot, storage)
            # (4) Solve the generalized eigenvalue problem
            wc = Basics.diagonalize("generalized eigenvalues: LinearAlgebra", wa, wb)
            # (5) Analyse and print information about the convergence of the symmetry blocks and the occupied orbitals
            wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
            if  wcBlock > 1.000 * settings.accuracyScf   go_on = true   end
            for  sh in basis.subshells
                if      sh in settings.frozenSubshells   ## do nothing
                elseif  sh.kappa == kappa
                    newOrbital = BsplinesN.generateOrbitalFromPrimitives(sh, wc, primitives)
                    wcOrbital  = Basics.analyzeConvergence(previousOrbitals[sh], newOrbital)
                    if  wcOrbital > settings.accuracyScf   accuracyScf = wcOrbital;   go_on = true   end
                        sa = "  $sh::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "  
                        sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                        sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                        if  printout    println(sa)    end
                    ##x println("  $sh  en [a.u.] = $(newOrbital.energy)   self-consistency = $(wcOrbital), $(wcBlock) [kappa=$kappa] ") 
                    previousOrbitals[sh] = newOrbital
                end
            end
            # (6) Re-define the bsplineBlock
            bsplineBlock[kappa] = wc
        end
        if  go_on   nothing   else   break   end
    end
    
    analyzedOrbitals = Basics.analyze(previousOrbitals, printout=true)
    
    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, analyzedOrbitals)
    return( newBasis )
end

end # module
    
