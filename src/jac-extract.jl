
"""
`Basics.extractNoOpenShells(conf::Configuration)`  
    ... determine the number of open (nonrelativistic) shells in conf; a singleton of type AbstractOpenShell is returned.
"""
function Basics.extractNoOpenShells(conf::Configuration)
    ##x println("conf = $conf")
    ns = 0
    for  (shell, occ)  in conf.shells
        if      occ == 0  ||  occ == 2*(2*shell.l + 1)
        else    ns = ns + 1
        end
    end
    ##x println("conf = $conf  ns = $ns")
    
    if      ns == 0   return ( LSjj.ZeroOpenShell() )
    elseif  ns == 1   return ( LSjj.OneOpenShell() )
    elseif  ns == 2   return ( LSjj.TwoOpenShells() )
    elseif  ns == 3   return ( LSjj.ThreeOpenShells() )
    else    error("stop a")
    end
    
end


"""
`Basics.extractNonrelativisticShellList(subshells::Array{Subshell,1})  
    ... extract a nonrelativistic shell list from the given subshell list; the shells appear in the same order
        as the corresponding subshells. A shellList::Array{Subshell,1} is returned.
"""
function Basics.extractNonrelativisticShellList(subshells::Array{Subshell,1}) 
    shellList = Shell[]

    for  subsh  in subshells
        l = Basics.subshell_l(subsh);   shell = Shell( subsh.n, l)
        if  shell in shellList      else  push!(shellList, shell)   end
    end 
    ##x println("Basics.extractNonrelativisticShellList()::  shellList = $shellList")
    
    return( shellList )
end


"""
`Basics.extractNonrelativisticConfigurations(basis::Basis)  
    ... extract all nonrelativistic configurations that contribute to the given set of CSF in basis.csfs. 
        A confList::Array{Configuration,1} is returned.
"""
function Basics.extractNonrelativisticConfigurations(basis::Basis)
    confList  = Configuration[]
    subshells = basis.subshells
    
    for  csf  in basis.csfs
        newShells = Dict{Shell,Int64}();    NoElectrons = 0
        for  i = 1:length(subshells)
            n = subshells[i].n;    l = Basics.subshell_l(subshells[i]);    occ = csf.occupation[i]
            shell = Shell(n,l)
            if    haskey(newShells, shell)    newShells[shell] = newShells[shell] + occ;   NoElectrons = NoElectrons + occ
            else  
                if   occ > 0  newShells = Base.merge( newShells, Dict( shell => occ));     NoElectrons = NoElectrons + occ   end
            end
        end
        if  basis.NoElectrons != NoElectrons    error("stop a")
        else
            conf = Configuration(newShells, NoElectrons)
            if  conf in confList    continue;    else    push!( confList,  conf)    end
        end
    end 
    ##x println("Basics.extractNonrelativisticConfigurations()::  confList = $confList")
    
    return( confList )
end



"""
`Basics.extractNonrelativisticConfigurationFromCsfR(csf::CsfR,  basis::Basis)`  
    ... extract the nonrelativistic configuration from the given CSF that is defined in basis.csfs. 
        A conf::Configuration is returned.
"""
function  Basics.extractNonrelativisticConfigurationFromCsfR(csf::CsfR,  basis::Basis)
    subshells = basis.subshells
    
    newShells = Dict{Shell,Int64}();    NoElectrons = 0
    for  i = 1:length(subshells)
        n = subshells[i].n;    l = Basics.subshell_l(subshells[i]);    occ = csf.occupation[i]
        shell = Shell(n,l)
        if    haskey(newShells, shell)    newShells[shell] = newShells[shell] + occ;   NoElectrons = NoElectrons + occ
        else  
            if   occ > 0  newShells = Base.merge( newShells, Dict( shell => occ));     NoElectrons = NoElectrons + occ   end
        end
    end
    
    if      basis.NoElectrons != NoElectrons    error("stop a")
    else    conf = Configuration(newShells, NoElectrons)
    end 
    ##x println("Basics.extractNonrelativisticConfigurationFromCsfR()::  conf = $conf")
    
    return( conf )
end


"""
`Basics.extractOpenShellQNfromCsfR(csfR::CsfR, basis::Basis)`  
    ... extracts the quantum numbers of the nonrelativistic open shells from the given csfr. Here, we assume that csfR is 
        defined in the given basis, and it is checked that the subshells of each nonrelativistic shell occur subsequently
        in basis and in the right sequence j = l - 1/2, j = l + 1/2, respectively. 
        An rQN::Dict{Shell, Tuple{Tuple{Int64,Int64,Int64,Int64,Int64}, Tuple{Int64,Int64,Int64,Int64,Int64}}}
                = Dict(shell => (rQNm, rQNp) ) is returned with rQNm = (idxm, Nm, Qm, Jm, Jmx) and 
        rQNp = (idxp, Np, Qp, Jp, Jpx). The two tuples rQNm, rQNp contains all the subshell quantum numbers and their
        coupling in csfR.
"""
function Basics.extractOpenShellQNfromCsfR(csfR::CsfR, basis::Basis)
    wa = Dict{Shell,Tuple{NTuple{5,Int64},NTuple{5,Int64}}}()
    for  s = 1:length(basis.subshells)
        sx = s + 1;     subsh = basis.subshells[s];     n = subsh.n;        occ   = csfR.occupation[s]
        l  = Basics.subshell_l(subsh);    j2 = Basics.subshell_2j(subsh);   shell = Shell(n,l)
        # Get l,j from next subshell ... if still defined
        if  sx <= length(basis.subshells)       subshx = basis.subshells[sx]
            nx = subshx.n;    lx = Basics.subshell_l(subshx);   j2x = Basics.subshell_2j(subshx);   occx   = csfR.occupation[sx]
        end
        # Now collect the quantum numbers into tuples
        if  l == 0   &&   0 < occ < 2
            ## wa = Base.merge( wa, Dict( shell => ( (s, csfR.occupation[s], csfR.seniority[s], 
            ##                                        Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ), 
            ##                                        (-9, 0, -9, 0, Basics.twice(csfR.subshellX[s]) ) ) ) )
            wa = Base.merge( wa, Dict( shell => ( (-9, 0, -9, 0, 0 ),
                                                  (s, csfR.occupation[s], csfR.seniority[s], 
                                                   Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ) ) ) )
        elseif  2l + 1 == j2    continue
        elseif  2l - 1 == j2    &&    (n != nx  ||  l != lx  ||  j2x != j2 + 2)
                error("stop a")
        elseif  0 < occ + occx < 2* (2l + 1)
            wa = Base.merge( wa, Dict(shell => ( (s, csfR.occupation[s], csfR.seniority[s], 
                                                  Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ), 
                                                 (sx, csfR.occupation[sx], csfR.seniority[sx], 
                                                  Basics.twice(csfR.subshellJ[sx]), Basics.twice(csfR.subshellX[sx]) ) ) ) )
        end
    end
    
    ##x println("Basics.extractOpenShellQNfromCsfR()::  csfR = $csfR  \n  wa = $wa")
    
    return( wa )
end


"""
`Basics.extractShellOccupationFromCsfR(csfR::CsfR, basis::Basis)`  
    ... extracts the nonrelativistic shell list as well as the (nonrelativistic) occupation of the given CsfR. 
        Here, we assume that csfR is defined in the given basis and that all subshells occur in standard sequence
        for a jjJ-LSJ transformation. A tuple tp::Tuple{Array{Shell,1},Array{Int64,1}}  of two arrays of equal 
        length is returned that contains the nonrelativistic shells and their occupation in the specification of 
        csfR. 
"""
function Basics.extractShellOccupationFromCsfR(csfR::CsfR, basis::Basis)
    shellList = Basics.extractNonrelativisticShellList(basis.subshells) 
    occList   = Int64[]
    
    for  s = 1:length(basis.subshells)
        l = Basics.subshell_l(basis.subshells[s]);   j2 = Basics.subshell_2j(basis.subshells[s])
        if      l == 0          push!( occList, csfR.occupation[s])
        elseif  2*l + 1 == j2   continue    # contributions from j = l + 1/2 are collected together with j = l - 1/2 below
        else    push!( occList, csfR.occupation[s] + csfR.occupation[s+1])
        end
    end 
    ##x println("Basics.extractShellOccupationFromCsfR()::  shellList = $shellList,  occList = $occList")
    
    return( (shellList, occList) )
end
