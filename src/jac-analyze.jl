
"""
`Basics.analyze()`  ... analyzes various data structures with regard to requested details.

  + `("level decomposition: % of NR configurations", level::Level)`  
    ... to analyze and list the non-relativistic configurations with a weight larger than 5%. A list of NR configurations with their 
        corresponding weights are printed, but nothing is returned.
"""
function Basics.analyze(sa::String, level::Level)
    !(sa == "level decomposition: % of NR configurations")   &&   error("Unsupported keystring = $sa")
    !(level.hasBasis)                                        &&   error("Levels without a basis cannot be analyzed.")
    
    confList   = generate("configuration list: NR, from basis", level.basis)
    percentage = zeros( length(confList) )
    
    # Now generate for each CSF in this basis a Configuration, compare with confList and add contribution to the corresponding percentage
    for  k = 1:length(level.basis.csfs)
        csf = level.basis.csfs[k]
        shellList = Dict{Shell,Int64}[]
        for  i = 1:length(csf.occupation)
            n = csf.subshell[i].n;    l = Basics.subshell_l( csf.subshell[i] );    sh = Shell(n, l)
            if     haskey(shellList, sh)    shellList(sh) = shellList(sh) + csf.occupation[i]
            else   shellList = merge( shellList, Dict(sh => csf.occupation[i]) )
            end
        end
        confNew = Configuration( shellList, NoElectrons )
        # Compare and add the weight to the right configuration
        for  i = 1:length(confList)
            if   confList[i] == confNew    percentage[i] = percentage[i] + level.mc[k]^2;    break    end 
        end
    end
  
    # Check that all weights sum up properly
    wa = sum( percentage );   abs(1. - wa) > 1.0e-3  &&    error("Total percentage $wa must add to 1.")
    # Sort the percentage and print both, configurations and percentage
    wb = sortperm( percentage )
    for i in wb
        println("Configuration: $(confList[i])  ... with  $(percentage[i]) %")
    end

    return( nothing )
end



"""
  + `("level decomposition: % of jj-coupled CSF", level::Level, N::Int64)`  
    ... to anaylze and list (up to) N relativistic CSF, together with their weight |c_n|^2 in the expansion of the given level. 
        A list of CSF and their corresponding weights are printed, but nothing is returned otherwise.  **Not yet implemented !**
"""
function Basics.analyze(sa::String, level::Level, N::Int64)
    !(sa == "level decomposition: % of jj-coupled CSF")   &&   error("Unsupported keystring = $sa")
    error("Not yet implemented !")

    return( nothing )
end


    
"""
`Basics.analyzeConvergence()`  
    ... analyzes the convergence for two (subsequently generated) orbitals or related quantities; a value::Float64 is
        returned if not indicated otherwise.

  + `(blocka::JAC.Eigen, blockb::JAC.Eigen)`  
    ... to analyze the convergence for two blocks of eigenvalues und eigenvectors; the value > 0 indicates the relative maximum 
        difference abs( a-b/a+b ) between two eigenvalues. The value is set arbitrarely set to 1.0e+2 if the number of eigenvalues
        differ in the two blocks.
"""
function Basics.analyzeConvergence(blocka::Basics.Eigen, blockb::Basics.Eigen)
    va = blocka.values;   vb = blockb.values
    if       length(va) !=  length(vb)    wa = 1.0e+2   
    else     wa = 0.;       for  i = 1:length(va)   wa = max(wa, abs( (va[i]-vb[i]) / (va[i]+vb[i]) ))   end
    end
    return( wa )
end


    
"""
  + `(a::JAC.Radial.Orbital, b::JAC.Radial.Orbital)`  
    ... to analyze the convergence for two radial orbitals; the value > 0 indicates the relative maximum difference abs( a-b/a+b ) 
        between the energies or (large or small) components of the orbitals. The function assumes that the orbitals are represented 
        on the same grid and that the components are zero if only one of the orbitals is defined beyond some maximum grid point.
"""
function Basics.analyzeConvergence(a::JAC.Radial.Orbital, b::JAC.Radial.Orbital)
    if       a.subshell !=  b.subshell                                             error("stop a")   
    elseif   a.useStandardGrid !=  b.useStandardGrid   ||  !(a.useStandardGrid)    error("stop b")   
    else     wa = 0.;   wa = max(wa, abs( (a.energy-b.energy) / (a.energy+b.energy) ));   nx = min(length(a.P), length(b.P))
             #=== for  i = 1:nx 
                 wp = abs(a.P[i])+abs(b.P[i]);    if  wp > 1.0e-10   wa = max(wa, abs( (a.P[i]-b.P[i]) / wp ))   end 
                 ##x if  abs( (a.P[i]-b.P[i]) / wp ) > 1.01    Defaults.warn(AddWarning, "analyzeConv: P, i = $i Pa=$(a.P[i]) Pb=$(b.P[i])")   end
                 wp = abs(a.Q[i])+abs(b.Q[i]);    if  wp > 1.0e-10   wa = max(wa, abs( (a.Q[i]-b.Q[i]) / wp ))   end   
                 ##x if  abs( (a.Q[i]-b.Q[i]) / wp ) > 1.01    Defaults.warn(AddWarning, "analyzeConv: Q, i = $i Qa=$(a.Q[i]) Qb=$(b.Q[i]) ")   end
             end ===#
    end
    return( wa )
end

