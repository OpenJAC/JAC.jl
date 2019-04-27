using   JAC.BasicTypes
export  generate

"""
`JAC.Basics.generate()`  ... generates various (complex) lists and elaborate data from simple rather simple input; the return typically depends
                      on the given keystring and the input data. 

  + `("condensed multiplet: by single weight", multiplet::Multiplet)`  ... to condense/reduce the number of CSF in the basis of the
                      given multiplet due to a single 'weight'; a multiplet::Multiplet is returned.  **Not yet implemented !**
"""
function generate(sa::String, multiplet::Multiplet)
    !(sa == "condensed multiplet: by single weight")   &&   error("Unsupported keystring = $sa")
    error("Not yet implemented !")
end



"""
  + `("configuration list: NR, from basis", basis::Basis)`  ... to (re-) generate the list of NR configurations from the given basis;
                      a confList::Array{Configuration,1} is returned.
"""
function generate(sa::String, basis::Basis)
    !(sa == "configuration list: NR, from basis")   &&   error("Unsupported keystring = $sa")
    confList    = Configuration[]
    NoElectrons = sum( basis.csfs[1].occupation )
    for  csf in basis.csfs
        shellList = Dict{Shell,Int64}[]
        for  i = 1:length(csf.occupation)
             n = csf.subshell[i].n;    l = subshell_l( csf.subshell[i] );    sh = Shell(n, l)
            if     haskey(shellList, sh)    shellList(sh) = shellList(sh) + csf.occupation[i]
            else   shellList = merge( shellList, Dict(sh => csf.occupation[i]) )
            end
        end
        confNew = Configuration( shellList, NoElectrons )
        push!( confList, confNew )
    end
    confListNew = ConfigurationsExcludeDoubles(confList)
    
    return( confListNew )
end


"""
  + `("configuration list: NR, single-configuration", refConf::Configuration, NoExcitations::Int64, fromShells::Array{Shell,1},`
                      toShells::Array{Shell,1})  ... to generate a non-relativistic configuration list, including the given reference
                      configuration (refConf) and with all configurations that differ by NoExcitations from the fromShells into the 
                      toShells; an Array{Configuration,1} is returned.
"""
function generate(sa::String, refConf::Configuration, NoExcitations::Int64, fromShells::Array{Shell,1}, toShells::Array{Shell,1})
    !(sa == "configuration list: NR, single-configuration")   &&   error("Unsupported keystring = $sa")
    confList = [refConf]
    # First prepare a proper reference configuration that also includes all fromShells and toShells with zero occupation
    shellDict = deepcopy(refConf.shells)
    allShells = deepcopy(fromShells)
    for  sh in toShells
        if     sh in allShells
        else   append!(allShells, [sh] )
        end
    end
    
    for  sh in allShells
       if  !haskey(shellDict, sh )    shellDict = merge( shellDict, Dict( sh => 0 ) )    end
    end

    rconf = Configuration(shellDict, refConf.NoElectrons)
    println("aa", rconf)

    # Add single excitations
    if      NoExcitations == 0
    elseif  NoExcitations == 1
        for  fromsh in fromShells
            for  tosh in toShells
                confnew = deepcopy(rconf)
                confnew.shells[fromsh] = confnew.shells[fromsh] - 1
                confnew.shells[tosh]   = confnew.shells[tosh]   + 1
                println("bb", confnew)

                # Check that this is valid configuration and append, if appropriate
                add = true
                for (k,v) in  confnew.shells
                    if  v <= -1  ||  Parity(rconf) != Parity(confnew)    add = false    end
                end
                if  add 
                    for (k,v) in confnew.shells   
                       if  v == 0   delete!(confnew.shells, k)   end
                    end
                    push!(confList, confnew);   println("fromsh = $fromsh  tosh = $tosh")   
                end
            end
        end

    # Add double excitations
    elseif  NoExcitations == 2
        for  fromsha in fromShells
        for  fromshb in fromShells
            for  tosha in toShells
            for  toshb in toShells
                confnew = deepcopy(rconf)
                confnew.shells[fromsha] = confnew.shells[fromsha] - 1
                confnew.shells[fromshb] = confnew.shells[fromshb] - 1
                confnew.shells[tosha]   = confnew.shells[tosha]   + 1
                confnew.shells[toshb]   = confnew.shells[toshb]   + 1
                println("bb", confnew)

                # Check that this is valid configuration and append, if appropriate
                add = true
                for (k,v) in  confnew.shells
                    if  v <= -1  ||  Parity(rconf) != Parity(confnew)    add = false    end
                end
                if  add    
                    for (k,v) in confnew.shells   
                       if  v == 0   delete!(confnew.shells, k)   end
                    end
                    push!(confList, confnew)    
                end
            end
            end
        end
        end

    else
        error("stop a")
    end

    # Exclude all configurations that appear twice (and more) in the list
    confList = JAC.ConfigurationsExcludeDoubles(confList)

    return( confList )
end


"""
  + `("configuration list: relativistic", conf::Configuration)`  ... to split/decompose a non-relativistic configuration into an list
                      of relativistic ConfigurationR[]. The proper occupuation of the relativistic subshells is taken into account.
"""
function generate(sa::String, conf::Configuration)
    subshellList = Subshell[]
    confList     = ConfigurationR[]

    !(sa == "configuration list: relativistic")  &&   error("Unsupported keystring = $sa")

    initialize = true;    NoElectrons = 0
    for (k, occ)  in conf.shells
        NoElectrons     = NoElectrons + occ
        wa              = JAC.BasicTypes.shellSplitOccupation(k, occ)
        subshellListNew = Dict{JAC.Subshell,Int64}[]

        if  initialize
            subshellListNew = wa;    initialize = false
        else
            for  s in 1:length(subshellList)
                for  a in 1:length(wa)
                    wb = Base.merge( subshellList[s], wa[a] )
                    push!(subshellListNew, wb)
                end
            end
        end
        subshellList = deepcopy(subshellListNew)
    end
    
    for subsh in subshellList
       wa = ConfigurationR(subsh, NoElectrons)
       push!(confList, wa)
    end

    return( confList )
end


"""
  + `("CSF list: from single ConfigurationR", conf::ConfigurationR, subshellList::Array{Subshell,1}) ... to construct from a given 
                 (relativistic) configuration all possible CSF with regard to the subshell order as specified by subshellList; 
                 a list::Array{CsfR,1} is returned.
"""
function generate(sa::String, conf::ConfigurationR, subshellList::Array{Subshell,1})
    parity  = Basics.determineParity(conf)
    csfList = CsfR[];   useStandardSubshells = true;    subhshellList = Subshell[];   first = true;    previousCsfs = CsfR[]
    # 
    for  subsh in subshellList
        if   subsh in keys(conf.subshells)    occ = conf.subshells[subsh]    else    occ = 0    end
        if   first
            stateList   = provide("subshell states: antisymmetric, seniority", subsh, occ)
            currentCsfs = CsfR[]
            for  state in stateList
                push!( currentCsfs, CsfR( true, AngularJ64(state.Jsub2//2), parity, [state.occ], [state.nu],
                                               [AngularJ64(state.Jsub2//2)], [AngularJ64(state.Jsub2//2)], Subshell[]) )
            end
            previousCsfs = copy(currentCsfs)
            first        = false
        else
            # Now support also all couplings of the subshell states with the CSFs that were built-up so far
            stateList   = provide("subshell states: antisymmetric, seniority", subsh, occ)
            currentCsfs = CsfR[]
            ##x for  i = 1:length(previousCsfs)    println("generate-aa: ", previousCsfs[i])    end
            for  csf in  previousCsfs
                for  state in stateList
                    occupation = deepcopy(csf.occupation);    seniority = deepcopy(csf.seniority);    
                    subshellJ  = deepcopy(csf.subshellJ);     subshells = deepcopy(csf.subshells)
                    push!(occupation, state.occ);   push!(seniority, state.nu);   push!(subshellJ, AngularJ64(state.Jsub2//2) ) 
                    push!(subshells, subsh)
                    newXList = oplus( csf.subshellX[end], AngularJ64(state.Jsub2//2) )
                    for  newX in newXList
                        subshellX = deepcopy(csf.subshellX);   push!(subshellX, newX) 
                        push!( currentCsfs, CsfR( true, subshellX[end], parity, occupation, seniority, subshellJ, subshellX, Subshell[]) ) 
                    end
                end
            end
            previousCsfs = copy(currentCsfs)
        end
    end
    
    return( previousCsfs )
end

  
"""
  + `("shells: ordered list for NR configurations", confs::Array{Configuration,1})`  ... to generate for confs, i.e. all the given
               (non-relativistic) configurations, a common and ordered shell list; a list::Array{Shell,1} is returned.
"""
function generate(sa::String, confs::Array{Configuration,1})
    shells = Shell[]   

    !(sa == "shells: ordered list for NR configurations")  &&   error("Unsupported keystring = $sa")

    wa = Constants.give("ordered shell list: non-relativistic", 7)
    for  a in wa
        for  cf in 1:length(confs)
            ks = keys(confs[cf].shells)
            if  a in ks   push!(shells, a);    break    end
        end 
    end

    return( shells )
end

  
"""
  + `("subshells: ordered list for relativistic configurations", confs::Array{ConfigurationR,1})`  ... to generate for confs, i.e. all the given
                      (relativistic) configurations, common and ordered subshell list; a list::Array{Subshell,1} is returned.
"""
function generate(sa::String, confs::Array{ConfigurationR,1})
    subshells = Subshell[]   

    !(sa == "subshells: ordered list for relativistic configurations")  &&   error("Unsupported keystring = $sa")

    wa = Constants.give("ordered subshell list: relativistic", 7)
    for  a in wa
        for  cf in 1:length(confs)
            ks = keys(confs[cf].subshells)
            ## Do include 'empty' subshells into the subshell list if they are specified by the given configurations
            ## if  a in ks   &&   confs[cf].subshells[a]  !=  0     push!(subshells, a);    break    end
            if  a in ks     push!(subshells, a);    break    end
        end 
    end
    ##x println("***subshells = $subshells")
    return( subshells )
end

  
"""
  + `("subshells: ordered list for two bases", basisA::Basis,  basisB::Basis)`  ... to generate common and ordered subshell 
                  list for the two basis A and B; a list::Array{Subshell,1} is returned.
"""
function generate(sa::String, basisA::Basis,  basisB::Basis)
    function areEqual(nn::Int64, sha::Array{Subshell,1}, shb::Array{Subshell,1})
        # Determines whether the first nn subshells are equal in sha and shb (true) or not (false)
        for  i = 1:nx   
            if    sha[i] != shb[i]    return( false )   end
        end
        return( true )
    end
        
    subshells = Subshell[]   

    !(sa == "subshells: ordered list for two bases")  &&   error("Unsupported keystring = $sa")

    nx = min(length(basisA.subshells), length(basisB.subshells))
    if  areEqual(nx, basisA.subshells, basisB.subshells)
        # If subshell order is equal, any subshell order is accepted for those subshells that occur in both basis
        for  i = 1:nx   
            if    basisA.subshells[i] != basisB.subshells[i]   error("Inconsistent subshells of two bases.")
            else  push!( subshells, basisA.subshells[i])
            end
        end
        #
        if       length(basisA.subshells) > nx   
            for  i = nx+1:length(basisA.subshells)    push!(subshells, basisA.subshells[i])    end
        elseif   length(basisB.subshells) > nx   
            for  i = nx+1:length(basisB.subshells)    push!(subshells, basisB.subshells[i])    end
        end
    else
        println("basisA.subshells = $(basisA.subshells)")
        println("basisB.subshells = $(basisB.subshells)")
        # If subshell order is NOT equal, all subshell in basisA and basisB must follow standard order
        standardList = Constants.give("ordered subshell list: relativistic", 7)
        na = 0;  nb = 0
        for  sh in standardList   
            if  sh in basisA.subshells  ||   sh in basisB.subshells   push!( subshells, sh)    end    
        end
        nn = 0;   
        for sh in basisA.subshells   
            if    !(sh in subshells) push!( subshells, sh)    
            else  wb = findall(x->x==sh, subshells);    if wb[1] <= nn   error("stop a")   else   nn = wb[1]  end    
            end
        end
        nn = 0;   
        for sh in basisB.subshells   
            if    !(sh in subshells) push!( subshells, sh)    
            else  wb = findall(x->x==sh, subshells);    if wb[1] <= nn   error("stop a")   else   nn = wb[1]  end    
            end
        end
        println("*** extended subshells from two basis = $subshells")
    end
    ##x println("subshells from two basis = $subshells")

    return( subshells )
end


  
"""
  + `("single-electron spectrum: STO", N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0,
                                       beta_0::Float64=1.1)` ... to generate a complete one-electron spectrum with N positive and N negative
                                       states, and by using even-tempered Slater-type orbitals (STO) with parameters 
                                       ``\alpha_i = \alpha_0 \beta_0^i``; a spectrum::SingleElecSpectrum is returned where just N0 positive 
                                       and N_0 negative are kept for later use.  **Not yet implemented !**

    `("single-electron spectrum: STO, positive", N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0,
                                      beta_0::Float64=1.1)`  ... to generate the same but to return only the N_0 positive states.  
                                      **Not yet implemented !**
"""
function generate(sa::String, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0, beta_0::Float64=1.1)
    !(sa == "single-electron spectrum: STO"  &&  sa == "single-electron spectrum: STO, positive")   &&   error("Unsupported keystring = $sa")
    error("Not yet implemented !")

    #= * define N normalized states in stoplus::Array{Vector{Float64},1} and stominus::Array{Vector{Float64},1} 
         on the given grid and with given parameters
       * calculate one-electron (NxN) Hamiltonian matrix with potential h_pq and overlap S_pq, p,q = 1..N
       * solve generalized eigenvalue problem h c_vec = epsilon S c_vec
       * set spectrum::SingleElecSpectrum with just N0 functions
       * check normalization and orthogonality of the functions   =#

    return( nothing )  
end


"""
`JAC.generateLevelWithExtraElectron(newOrbital, symt::LevelSymmetry, level::Level)`  ... generates a (new) level with one extra electron
                                    in subshell newOrbital.subshell and with overall symmetry symt. The function assumes that all CSF in the basis
                                    of level have the same symmetry as the level itself, that the new subshell is not yet part of the basis and
                                    that CSF with the total symmetry can be constructed; it terminates with an error if one of this assumptions
                                    is violated. A newLevel::Level with the same number of CSF and the same representation (eigenvector)
                                    as the given level is returned. From a physics viewpoint, newLevel refers to a antisymmetrized product state
                                    of [level x orbital(sh, energy)] J^P (symt) and with an total energy = level.energy + newOrbital.energy.
"""
function generateLevelWithExtraElectron(newOrbital, symt::LevelSymmetry, level::Level)
    basis = level.basis;    newSubshells = copy(basis.subshells);    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    
    push!(newSubshells, newOrbital.subshell)
    Constants.define("relativistic subshell list", newSubshells; printout=false)
    newOrbitals = Base.merge( copy(basis.orbitals), Dict( newOrbital.subshell => newOrbital ))
   
    for  i = 1:length(basis.csfs)
        stateList   = provide("subshell states: antisymmetric, seniority", newOrbital.subshell, 1)
        if  basis.csfs[i].J != J   ||   basis.csfs[i].parity != parity    error("Improper symmetry of CSF.")   end
        if  length(stateList) != 1                                        error("Improper number of subshell states.")   end
        if  JAC.AngularMomentum.triangularDelta(J, BasicTypes.subshell_j(newOrbital.subshell), symt.J) != 1    
                                                                         error("Improper coupling of (new) subshell & total J.")   end
        substate   = stateList[1]
        occupation = copy(basis.csfs[i].occupation);    push!(occupation, substate.occ)
        seniority  = copy(basis.csfs[i].seniority);     push!(seniority,  substate.nu)
        subshellJ  = copy(basis.csfs[i].subshellJ);     push!(subshellJ,  AngularJ64(substate.Jsub2//2) )
        subshellX  = copy(basis.csfs[i].subshellX);     push!(subshellX,  symt.J)
        push!(newCsfs, CsfR(true, symt.J, symt.parity, occupation, seniority, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons+1, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(symt.J, AngularM64(symt.J), symt.parity, level.index, level.energy+newOrbital.energy, level.relativeOcc, 
                      true, newBasis, level.mc)
    return( newLevel )
end


"""
`JAC.Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level)`  ... generates a (new) level with one extra subshell sh but with the same
                                    overall symmetry as before. The function assumes that the new subshell is not yet part of the basis; 
                                    it terminates with an error if one of this assumptions is violated. A newLevel::Level with the same 
                                    symmetry, number of CSF, energy and the same representation (eigenvector) as the given level is returned. 
"""
function generateLevelWithExtraSubshell(sh::Subshell, level::Level)
    basis = level.basis;    newSubshells = copy(basis.subshells);    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    newOrbitals = copy(basis.orbitals)

    push!(newSubshells, sh)
    Constants.define("relativistic subshell list", newSubshells; printout=false)
   
    for  i = 1:length(basis.csfs)
        stateList   = provide("subshell states: antisymmetric, seniority", sh, 0)
        if  length(stateList) != 1                                        error("Improper number of subshell states.")   end
        substate   = stateList[1]
        occupation = copy(basis.csfs[i].occupation);    push!(occupation, substate.occ)
        seniority  = copy(basis.csfs[i].seniority);     push!(seniority,  substate.nu)
        subshellJ  = copy(basis.csfs[i].subshellJ);     push!(subshellJ,  AngularJ64(substate.Jsub2//2) )
        subshellX  = copy(basis.csfs[i].subshellX);     push!(subshellX,  subshellX[end])
        push!(newCsfs, CsfR(true, basis.csfs[i].J, basis.csfs[i].parity, occupation, seniority, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(level.J, level.M, level.parity, level.index, level.energy, level.relativeOcc, true, newBasis, level.mc)
    return( newLevel )
end


"""
`JAC.Basics.generateLevelWithSymmetryReducedBasis(level::Level)`  ... generates a level with a new basis and representation that only includes the
                                                                CSF with the same symmetry as the level itself. Both, the underlying CSF basis
                                                                and the eigenvectors are adopted accordingly. A (new) level::levelR is returned.
"""
function generateLevelWithSymmetryReducedBasis(level::Level)
    basis = level.basis;    newCsfs = CsfR[];    newMc = Float64[]
    
    for  i = 1:length(basis.csfs)
        if  basis.csfs[i].J == level.J  &&   basis.csfs[i].parity == level.parity 
            push!(newCsfs, basis.csfs[i]);    push!(newMc, level.mc[i])
        end
    end
    newBasis = Basis(true, basis.NoElectrons, basis.subshells, newCsfs, basis.coreSubshells, basis.orbitals)
    newLevel = Level(level.J, level.M, level.parity, level.index, level.energy, level.relativeOcc, true, newBasis, newMc)
    return( newLevel )
end



