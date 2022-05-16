    
    export recastAG

    """
    `Basics.add(pota::Radial.Potential, potb::Radial.Potential)`  
        ... to add two radial potentials together that are defined on the same grid. A potential::RadialPotential is returned 
            that inherits its radial size from the potential that is defined in a larger range of r-values.
    """
    function  Basics.add(pota::Radial.Potential, potb::Radial.Potential)
        if  pota.grid.NoPoints != potb.grid.NoPoints  ||   pota.grid.rnt != potb.grid.rnt  ||  pota.grid.h != potb.grid.h  ||  
            pota.grid.hp != potb.grid.hp                   error("stop a")   end
        name = pota.name * "+" * potb.name;   nmax = max(length(pota.Zr), length(potb.Zr));   nmin = min(length(pota.Zr), length(potb.Zr))
        Zr   = zeros(nmax);   
        nx = length(pota.Zr);    Zr[1:nx] = Zr[1:nx] + pota.Zr[1:nx] 
        nx = length(potb.Zr);    Zr[1:nx] = Zr[1:nx] + potb.Zr[1:nx]
        Zr[nmin+1:nmax] .= Zr[nmin] 

        potential = Radial.Potential(name, Zr, pota.grid)
        return( potential )
    end



    """
    `Basics.addZerosToCsfR(nz::Int64, csf::CsfR)`  
        ... 'adds' a number of zeros to the occupation, seniorityNr, etc. of a given CsfR.
    """
    function  Basics.addZerosToCsfR(nz::Int64, csf::CsfR)
        !csf.useStandardSubshells   &&   error("Zeros can only be added to CSF with standard subshell order.")
        
        occupation = csf.occupation;   seniorityNr = csf.seniorityNr;    subshellJ = csf.subshellJ;    subshellX = csf.subshellX 
        
        for  i = 1:nz
            push!(occupation, 0);    push!(seniorityNr, 0);    push!(subshellJ, AngularJ64(0));    push!(subshellX, subshellJ[end])
        end
        
        newCsf = CsfR( csf.useStandardSubshells, csf.J, csf.parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[])
        return(newCsf)
    end

    
    
    """
    `Basics.analyze("level decomposition: % of NR configurations", level::Level)`  
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
    `Basics.analyze("level decomposition: % of jj-coupled CSF", level::Level, N::Int64)`  
        ... to anaylze and list (up to) N relativistic CSF, together with their weight |c_n|^2 in the expansion of the given level. 
            A list of CSF and their corresponding weights are printed, but nothing is returned otherwise.  **Not yet implemented !**
    """
    function Basics.analyze(sa::String, level::Level, N::Int64)
        !(sa == "level decomposition: % of jj-coupled CSF")   &&   error("Unsupported keystring = $sa")
        error("Not yet implemented !")

        return( nothing )
    end



    """
    `Basics.analyze(orbitals::Dict{Subshell, Orbital}; printout::Bool=false)`  
        ... to anaylze the phase of the orbitals and to ensure that all orbitals start `positve' near r --> 0; 
            a set of newOrbitals::Dict{Subshell, Orbital} is returned.
    """
    function Basics.analyze(orbitals::Dict{Subshell, Orbital}; printout::Bool=false)
        newOrbitals = Dict{Subshell, Orbital}();    mOrbital = Orbital(Subshell("1s_1/2"), -1.0)
        # Analyze and modify the orbital if necessary
        for (k,v) in orbitals
            found = false
            for  P in v.P   
                if      abs(P) > 1.0e-5  &&  abs(P) == P
                    mOrbital = v;   found = true;   break
                elseif  abs(P) > 1.0e-5  &&  abs(P) == -P
                    mOrbital = Orbital( v.subshell, v.isBound, v.useStandardGrid, v.energy, -1.0 .* v.P,  -1.0 .* v.Q,  
                                        -1.0 .* v.Pprime,  -1.0 .* v.Qprime, v.grid)
                    if  printout    println(">>> Sign changed for orbital $k")   end                    
                    found = true;   break
                end
            end
            if  !found  error("stop a")     end
            newOrbitals[k] = mOrbital
        end

        return( newOrbitals )
    end


        
    """
    `Basics.analyzeConvergence(blocka::Basics.Eigen, blockb::Basics.Eigen)`  
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
    `Basics.analyzeConvergence(a::Radial.Orbital, b::Radial.Orbital)`  
        ... to analyze the convergence for two radial orbitals; the value > 0 indicates the relative maximum difference abs( a-b/a+b ) 
            between the energies or (large or small) components of the orbitals. The function assumes that the orbitals are represented 
            on the same grid and that the components are zero if only one of the orbitals is defined beyond some maximum grid point.
    """
    function Basics.analyzeConvergence(a::Radial.Orbital, b::Radial.Orbital)
        if       a.subshell !=  b.subshell                                             error("stop a")   
        elseif   a.useStandardGrid !=  b.useStandardGrid   ||  !(a.useStandardGrid)    @show a.useStandardGrid, b.useStandardGrid;   error("stop b")   
        else     wa = 0.;   wa = max(wa, abs( (a.energy-b.energy) / (a.energy+b.energy) ));   nx = min(length(a.P), length(b.P))
                #=== for  i = 1:nx 
                    wp = abs(a.P[i])+abs(b.P[i]);    if  wp > 1.0e-10   wa = max(wa, abs( (a.P[i]-b.P[i]) / wp ))   end 
                    wp = abs(a.Q[i])+abs(b.Q[i]);    if  wp > 1.0e-10   wa = max(wa, abs( (a.Q[i]-b.Q[i]) / wp ))   end   
                end ===#
        end
        return( wa )
    end

        

    """
    `Basics.determineEnergySharings(energy::Float64, NoEnergySharings::Int64)`  
        ... to determine the NoEnergySharings sharings by using the Gauss-Legendre integration points in the interval
            [0., energy]; an Array{Tuple{Float64,Float64,Float64},1} is returned. This methods assumes that the package 
            GaussQuadrature is 'used'.
    """
    function Basics.determineEnergySharings(energy::Float64, NoEnergySharings::Int64)
        xx, ww = GaussQuadrature.legendre(NoEnergySharings)
        sharings = Tuple{Float64,Float64,Float64}[]
        for  (ix, x) in enumerate(xx)
            xrel = (x+1.) / 2
            push!(sharings, (xrel*energy, energy - xrel*energy, ww[ix]*energy) )
        end
        return( sharings )
    end



    """
    `Basics.determineHoleShells(conf::Configuration)`  
        ... to determine the hole-shells of a given non-relativistic configuration; a shellList::Array{Shell,1} is returned.
    """
    function Basics.determineHoleShells(conf::Configuration)

        shellList = Shell[]
        for  (k,v) in conf.shells
        if   conf.shells[k] < 2*(2k.l+1)   push!(shellList, k)    end 
        end

        return( shellList )
    end



    """
    `Basics.determineMeanEnergy(conf::Configuration, orbitals::Dict{Subshell, Orbital}, nm::Nuclear.Model, grid::Radial.Grid)`  
        ... to determine the mean energy of a given non-relativistic configuration; this 'mean' energy is calculated as the mean
            energy of all single-CSF levels and by using the the given set of orbitals. A value::Float64 is returned.
    """
    function Basics.determineMeanEnergy(conf::Configuration, orbitals::Dict{Subshell, Orbital}, nm::Nuclear.Model, grid::Radial.Grid)
        asfSettings = AsfSettings();    meanEnergy = 0.;    nlev = 0
        multiplet   = Basics.performCI([conf], orbitals, nm, grid, asfSettings, printout=false)
        for  level in multiplet.levels
            nl = Basics.twice(level.J) + 1;     nlev = nlev + nl;     meanEnergy = meanEnergy + nl * level.energy
        end  
        meanEnergy = meanEnergy / nlev
        
        if true  println(">>> mean energy of $conf  is:  $meanEnergy  [a.u.]")   end
        
        return( meanEnergy )
    end



    """
    `Basics.determineMeanEnergy(multiplet::Multiplet)`  
        ... to determine the mean energy of a given multiplet. A value::Float64 is returned.
    """
    function Basics.determineMeanEnergy(multiplet::Multiplet)
        meanEnergy = 0.;    nlev = 0
        for  level in multiplet.levels
            nl = Basics.twice(level.J) + 1;     nlev = nlev + nl;     meanEnergy = meanEnergy + nl * level.energy
        end  
        meanEnergy = meanEnergy / nlev
        
        if false  println(">>> mean energy of $(nlev)-level multiplet is:  $meanEnergy  [a.u.]")   end
        
        return( meanEnergy )
    end


    """
    `Basics..determineNonorthogonalShellOverlap(finalSubshells::Array{Subshell,1}, initialSubshells::Array{Subshell,1}, 
                                                grid::Radial.Grid)`  
        ... determines the overlap of those shells that do not agree in both sets but have the same symmetry; a zero
            overlap and Subshell(99,-99) are returned in cases that no such pair can be found.
            A triple (shella, shellb, overlap) is returned.
    """
    function Basics.determineNonorthogonalShellOverlap(finalSubshells::Array{Subshell,1}, initialSubshells::Array{Subshell,1}, 
                                                       grid::Radial.Grid)
        neSubshells = setdiff(union(finalSubshells, initialSubshells),  intersect(finalSubshells, initialSubshells))
        ret = (Subshell(99,-99), Subshell(99,-99), 0.);     ne = length(neSubshells)
        ##x @show  neSubshells
        for  a = 1:ne
            for  b = a+1:ne
                if   neSubshells[a].kappa == neSubshells[b].kappa   
                    println(">>>> non-zero overlap integral < $(neSubshells[a]) | $(neSubshells[b]) > ")
                    ret = (neSubshells[a], neSubshells[b], 0.1)
                end
            end
        end
        
        if  ret[3] == 0.  println(">>>> zero overlap integral for $neSubshells ")  end
        
        return( ret )
    end



    """
    `Basics.determineParity(conf::Configuration)`  ... to determine the parity of a given non-relativistic configuration.
    """
    function Basics.determineParity(conf::Configuration)

        par = 1
        for  (k,v) in conf.shells
        if   iseven(k.l)    p = 1   else   p = -1    end    
        par = par * (p^v)
        end

        if       par == 1    return( Basics.plus )  
        elseif   par == -1   return( Basics.minus )
        else     error("stop b")
        end  
    end


    
    """
    `Basics.determineParity(conf::ConfigurationR)`  ... to determine the parity of a given relativistic configuration.
    """
    function Basics.determineParity(conf::ConfigurationR)

        par = 1
        for  (k,v) in conf.subshells
        l = Basics.subshell_l(k)
        if   iseven(l)    p = 1   else   p = -1    end    
        par = par * (p^v)
        end

        if       par == 1    return( Basics.plus )  
        elseif   par == -1   return( Basics.minus )
        else     error("stop b")
        end  
    end


    """
    `Basics.determinePolarizationVector(q::Int64, kind::Basics.AbstractPolarization; star::Bool)`  
        ... determines the q-th component of the (complex and unit) polarization vector u for the given kind of polarization;
            for star = true, the complex-conjugate component is returned.  a comp::Float64 is returned.
    """
    #---OLD VERSIONS:---
    #function Basics.determinePolarizationVector(q::Int64, kind::Basics.LinearPolarization; star::Bool=false)
    #    if      q == 1   return( -1/sqrt(2.) )
    #    elseif  q == -1  return( 1/sqrt(2.) )
    #    else   return( 0 )
    #    end
    #end
    #
    #function Basics.determinePolarizationVector(q::Int64, kind::Basics.LeftCircular; star::Bool=false)
    #    if  q == 1   return( -1. )   else   return( 0. )                end
    #end
    
    #function Basics.determinePolarizationVector(q::Int64, kind::Basics.RightCircular; star::Bool=false)
    #    if  q == -1  return( 1. )   else   return( 0. )                end
    #end
    
    #---NEW VERSIONS:---
    function Basics.determinePolarizationVector(q::Int64, kind::Basics.LinearPolarization; star::Bool)
        if  star == false
            if      q == 1   return( -1/sqrt(2.) )
            elseif  q == -1  return( 1/sqrt(2.) )
            else   return( 0 )
            end
        else
            if      q == -1   return( 1/sqrt(2.) )
            elseif  q == 1  return( -1/sqrt(2.) )
            else   return( 0 )
            end
        end
    end
    #
    function Basics.determinePolarizationVector(q::Int64, kind::Basics.LeftCircular; star::Bool)
        if  star == false
            if  q == 1   return( -1. )   else   return( 0. )                end
        else
            if  q == -1   return( 1. )   else   return( 0. )                end
        end
    end
    #
    function Basics.determinePolarizationVector(q::Int64, kind::Basics.RightCircular; star::Bool)
        if  star == false
            if  q == -1  return( 1. )   else   return( 0. )                end
        else
            if  q == 1  return( -1. )   else   return( 0. )                end
        end
    end
    #
    function Basics.determinePolarizationVector(q::Int64, kind::Basics.LeftElliptical; star::Bool)
        epsilon = kind.ellipticity
        if  star == false
            if  q == 1   return( -1/sqrt(2*(1. + epsilon^2)) * (1. - epsilon) )   else   return( 1/sqrt(2*(1. + epsilon^2)) * (1. + epsilon) )                end
        else
            if  q == -1   return( 1/sqrt(2*(1. + epsilon^2)) * (1. - epsilon) )   else   return( -1/sqrt(2*(1. + epsilon^2)) * (1. + epsilon) )               end
        end
    end
    #
    function Basics.determinePolarizationVector(q::Int64, kind::Basics.RightElliptical; star::Bool)
        epsilon = kind.ellipticity
        if  star == false
            if  q == -1  return( 1/sqrt(2*(1. + epsilon^2)) * (1. - epsilon) )   else   return( -1/sqrt(2*(1. + epsilon^2)) * (1. + epsilon) )               end
        else
            if  q == 1  return( -1/sqrt(2*(1. + epsilon^2)) * (1. - epsilon) )   else   return( 1/sqrt(2*(1. + epsilon^2)) * (1. + epsilon) )                end
        end
    end

    """
    `Basics.determinePolarizationLambda(kind::Basics.AbstractPolarization)`  
        ... determines the (standard) lambda value for the given kind of polarization; a lambda::Int64 is returned.
    """
    function Basics.determinePolarizationLambda(kind::Basics.LeftCircular)   return(  1 )    end
    function Basics.determinePolarizationLambda(kind::Basics.RightCircular)  return( -1 )    end
    function Basics.determinePolarizationLambda(kind::Basics.LeftElliptical)   return(  1 )    end
    function Basics.determinePolarizationLambda(kind::Basics.RightElliptical)  return( -1 )    end


    """
    `Basics.determineSelectedLines(lineList::Array{Tuple{Int64,Int64},1}, initialMultiplet::Multiplet, finalMultiplet::Multiplet)`  
        ... to determine the specified lines as tupels of initial- and final levels. A level index is 0 in lineList always refers to
            all levels from the corresponding multiplet. A (unique) newList::Array{Tuple{Int64,Int64},1} is returned.
    """
    function Basics.determineSelectedLines(lineList::Array{Tuple{Int64,Int64},1}, initialMultiplet::Multiplet, finalMultiplet::Multiplet)
        newList = Tuple{Int64,Int64}[]
        for tupl in lineList
            if      tupl[1] >  length(initialMultiplet.levels)   error("Initial level index $(tupl[1]) must be <= No. initial levels.")
            elseif  tupl[2] >  length(finalMultiplet.levels)     error("Final level index $(tupl[2]) must be <= No. final levels.")
            elseif  tupl[1] == 0  &&   tupl[2] == 0    error("Tupel (0,0) is not allowed.")
            elseif  tupl[1] == 0       for i = 1:length(initialMultiplet.levels)   push!(newList, (i, tupl[2]))   end
            elseif  tupl[2] == 0       for i = 1:length(finalMultiplet.levels)     push!(newList, (tupl[1], i))   end
            else                                                                   push!(newList, (tupl[1], tupl[2]))
            end
        end
        newList = unique(newList)
        return( newList)
    end     


    """
    `Basics.determineSelectedPathways(pathwayList::Array{Tuple{Int64,Int64,Int64},1}, initialMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                      finalMultiplet::Multiplet)`  
        ... to determine the specified pathways as tupels of initial-, intermediate- and final levels. A level index 0 in pathwayList always 
            refers to all levels from the corresponding multiplet. A (unique) newList::Array{Tuple{Int64,Int64,Int64},1} is returned.
    """
    function Basics.determineSelectedPathways(pathwayList::Array{Tuple{Int64,Int64,Int64},1}, initialMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                    finalMultiplet::Multiplet)
        newList = Tuple{Int64,Int64,Int64}[]
        for tupl in pathwayList
            if      tupl[1] >  length(initialMultiplet.levels)       error("Initial level index $(tupl[1]) must be <= No. initial levels.")
            elseif  tupl[2] >  length(intermediateMultiplet.levels)  error("Intermediate level index $(tupl[2]) must be <= No. intermediate levels.")
            elseif  tupl[3] >  length(finalMultiplet.levels)         error("Final level index $(tupl[3]) must be <= No. final levels.")
            elseif  (tupl[1] + tupl[2] == 0  || tupl[1] + tupl[3] == 0  ||  tupl[2] + tupl[3] == 0)   &&  error("Only tupels with one 0-index are allowed.")
            elseif  tupl[1] == 0       for i = 1:length(initialMultiplet.levels)          push!(newList, (i, tupl[2], tupl[3]))   end
            elseif  tupl[2] == 0       for i = 1:length(intermediateMultiplet.levels)     push!(newList, (tupl[1], i, tupl[3]))   end
            elseif  tupl[3] == 0       for i = 1:length(finalMultiplet.levels)            push!(newList, (tupl[1], tupl[2], i))   end
            else                                                                          push!(newList, (tupl[1], tupl[2], tupl[3]))
            end
        end
        newList = unique(newList)
        return( newList)
    end     

    
    
    """
    `Basics.diagonalize("matrix: LinearAlgebra", matrixA::Array{Float64,2}; range=(0:0)::UnitRange{Int64})`  
        ... to apply the standard the standard LinearAlgebra.eigen() method from Julia for a symmetrc matrix; 
            only the upper-triangular parts of matrixA is used by an explicit symmetrization; ; an eigen::Basics.Eigen is returned.
    """
    function Basics.diagonalize(sa::String, matrixA::Array{Float64,2}; range=(0:0)::UnitRange{Int64})
        if       sa == "matrix: LinearAlgebra" 
            # Use the LinearAlgebra.eigen method from Julia for one symmetric matrix  
            mA = LinearAlgebra.Symmetric(matrixA)
            if  range == 0:0    
                wa = LinearAlgebra.eigen( mA )
                vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
                for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d])    end
            else                
                wa = LinearAlgebra.eigen( mA, range )
                vectors = Vector{Float64}[];    for  i = range   push!( vectors, wa.vectors[:,1] )   end
            end
            #
            wc = Basics.Eigen( wa.values, vectors )
            return( wc )
        #== elseif       sa == "matrix: Julia, eigfact" 
            # Use the standard eigfact() method from Julia for a quadratic, full matrix   
            wa = eigen( matrix )
            vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
            for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d])    end
            wc = Basics.Eigen( wa.values, vectors )
            return( wc )  ==#
        else     error("Unsupported keystring = $sa")
        end
    end


    """
    `Basics.diagonalize("generalized eigenvalues: LinearAlgebra", matrixA::Array{Float64,2}, matrixB::Array{Float64,2})`  
        ... to apply the standard LinearAlgebra.eigen() method from Julia for a generalized eigenvalue problem with two symmetric
            matrices; only the upper-triangular parts of matrixA and matrixB are used by an explicit symmetrization; 
            an eigen::Basics.Eigen is returned.
    """
    function Basics.diagonalize(sa::String, matrixA::Array{Float64,2}, matrixB::Array{Float64,2})
        if       sa == "generalized eigenvalues: LinearAlgebra" 
            # Use the LinearAlgebra.eigen method from Julia for two symmetric matrices   
            mA = LinearAlgebra.Symmetric(matrixA)
            mB = LinearAlgebra.Symmetric(matrixB)
            wa = LinearAlgebra.eigen(mA,mB)
            vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
            for   i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d] )    end
            wc      = Basics.Eigen( wa.values, vectors )
            return( wc )
        #== elseif       sa == "generalized eigenvalues: Julia, eigfact" 
            # Use the standard eigfact() method from Julia for two quadratic, full matrices   
            wa = eigen( matrixA, matrixB )
            vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
            # for  i = 0:d-1    push!(vectors, real(wb[i*d+1:i*d+d]) )    end
            for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d] )    end
            # wc = Basics.Eigen( real(wa.values), vectors )
            wc = Basics.Eigen( wa.values, vectors )
            return( wc )  ==#
        else     error("Unsupported keystring = $sa")
        end
    end


    """
    `Basics.diracDelta(x::Float64, dx::Float64)`  
        ... evaluates Dirac's function  delta(x) = 0  for abs(x) > dx/2   and   delta(x) = 1/dx    for abs(x) <= dx/2;
            a values::Float64 is returned.
    """
    function Basics.diracDelta(x::Float64, dx::Float64)    
        if  abs(x) > dx/2   return( 0. )    else    return( 1/dx )  end
    end
    

    """
    `Basics.display("constants")`  or  `("physical constants")`  
        ... to display (all) currently defined physical constants; nothing is returned if not indicated otherwise. 
            Cf. Defaults.setDefaults().

    `Basics.display("settings")`     ... to display (all) currently defined settings of the JAC module.
    """
    function Basics.display(sa::String)

        if        sa in ["constants", "physical constants"]
            println("Physical constants are defines as follows:  \n", 
                    "------------------------------------------  \n")
            sb = "  + Fine-structure constant:    " * string( Defaults.getDefaults("alpha") );                println(sb)        
            sb = "  + Electron mass [kg]:         " * string( Defaults.getDefaults("electron mass: kg") );    println(sb)        
            sb = "  + Electron mass [amu]:        " * string( Defaults.getDefaults("electron mass: amu") );   println(sb)        
            println()

        elseif   sa == "settings"
            println("Current settings of the JAC module:  \n", 
                    "-----------------------------------  \n")
            sb = "  + Framework:                              " * string( Defaults.getDefaults("framework") );                println(sb)        
            sb = "  + Energy unit:                            " * string( Defaults.getDefaults("unit: energy") );             println(sb)        
            sb = "  + Rate and transition probability unit:   " * string( Defaults.getDefaults("unit: rate") );               println(sb)        
            sb = "  + Cross section unit:                     " * string( Defaults.getDefaults("unit: cross section") );      println(sb) 
            sb = "  + Time unit:                              " * string( Defaults.getDefaults("unit: time") );               println(sb) 
            println()       
            
            if      Defaults.getDefaults("standard grid") != false    println("  + A standard grid has been defined; cf. Defaults.getDefaults()" )
            elseif  Defaults.getDefaults("standard grid") == false    println("  + No standard grid has yet been defined; cf. Defaults.setDefaults()" )
            end
            println()       
    
        else    error("Unsupported keystring:: $sa")
        end

        return( nothing )
    end

    
    
    """
    `Basics.display(stream::IO, configs::Array{Configuration,1})`  
        ... displays the generated list of configurations in a compact form; nothing is returned in this case.
    """
    function  Basics.display(stream::IO, configs::Array{Configuration,1})
        nx = 85
        println(" ")
        println("  Generated list of (non-relativistic) configurations that contribute to the representation:")
        ## println(" ")
        println("  ", TableStrings.hLine(nx))
        for  i = 1:length(configs)
            sa = string("(", i, ")" );    sa = TableStrings.flushright(9, sa, na=3);   sa = sa * string(configs[i])
            println(sa)
        end 
        println("  ", TableStrings.hLine(nx))
        
        return( nothing )
    end

    
    
    """
    `Basics.display(stream::IO, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid; longTable::Bool=false)`  
        ... displays the (generated)  orbitals in a neat and compact form; the default (longTable==false) just prints the subshell,
            isbound and energies. If longTable==true, it prints in addition the r_max, <r>, ... values.  nothing is returned.
    """
    function  Basics.display(stream::IO, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid; longTable::Bool=false)
        println(stream, " ")
        println(stream, "  Relativistic orbitals:")
        println(stream, " ")
        if  longTable
            nx = 108
            println(stream, "  ", TableStrings.hLine(nx))
            println(stream, "   Subshell   isBound   energy [a.u.]     energy " * TableStrings.inUnits("energy") * 
                            "   st-grid   r_max [a.u.]    <r> [a.u.]    <r^2> [a.u.]  ")
            println(stream, "  ", TableStrings.hLine(nx))
            #
            for  sh  in Defaults.GBL_STANDARD_SUBSHELL_LIST
                if haskey(orbitals, sh)
                    v = orbitals[sh]
                    sen_au = @sprintf("%.8e", v.energy)
                    sen    = @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", v.energy))
                    if !v.useStandardGrid   error("Non-standard grids not yet implemented.")    end
                    p2 = zeros( size(v.P, 1) );    p2 = p2 .+ (v.P .* v.P);    p2 = p2 .+ (v.Q .* v.Q)
                    imax    = findmax(p2)
                    rmax    = grid.r[imax[2]]
                    rmean   = RadialIntegrals.rkDiagonal(1, v, v, grid)
                    r2mean  = RadialIntegrals.rkDiagonal(2, v, v, grid)
                    srmax   = @sprintf("%.5e", rmax)
                    srmean  = @sprintf("%.5e", rmean)
                    sr2mean = @sprintf("%.5e", r2mean)
                    sa      = "     $sh    $(v.isBound)   " * sen_au * "   " * sen * "   $(v.useStandardGrid)     " 
                    sa      = sa * srmax * "    " * srmean * "    " * sr2mean
                    println(sa)
                end
            end 
            println("  ", TableStrings.hLine(nx))
        else
            nx = 62
            println(stream, "  ", TableStrings.hLine(nx))
            println(stream, "   Subshell   isBound   energy [a.u.]     energy " * TableStrings.inUnits("energy") * "   st-grid ")
            println(stream, "  ", TableStrings.hLine(nx))
            #
            for  sh  in Defaults.GBL_STANDARD_SUBSHELL_LIST
                if haskey(orbitals, sh)
                    v = orbitals[sh]
                    sen_au = @sprintf("%.8e", v.energy)
                    sen    = @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", v.energy))
                    sa     = "     $sh    $(v.isBound)   " * sen_au * "   " * sen * "   $(v.useStandardGrid)"
                    println(sa)
                end
            end 
            println("  ", TableStrings.hLine(nx))
        end
        
        return( nothing )
    end

    
    
    """
    `Basics.display(stream::IO, channels::Array{AtomicState.GreenChannel,1})`  
        ... displays the (generated Green function)  channels in a neat and compact form; nothing is returned in this case.
    """
    function  Basics.display(stream::IO, channels::Array{AtomicState.GreenChannel,1})
        println(stream, " ")
        println(stream, "  Green function channels:")
        println(stream, " ")
        for  channel in channels
            energies = Float64[]
            for  level in channel.gMultiplet.levels   push!(energies, level.energy)  end
            minenergy = minimum(energies);    maxenergy = maximum(energies)
            println(stream, "  Channel with $(channel.symmetry) symmetry, $(length(channel.gMultiplet.levels)) levels and energies [Hartree]:" *
                            " $minenergy ... $maxenergy")
        end
        
        return( nothing )
    end
   

    """
    `Basics.displayLevels(stream::IO, multiplets::Array{Multiplet,1}; N::Int64=10)`  
        ... display on stream the level energies of the N lowest levels in ascending order, and together with the configuration 
            of most leading term in the level expansion. A neat table is printed but nothing is returned otherwise.
    """
    function Basics.displayLevels(stream::IO, multiplets::Array{Multiplet,1}; N::Int64=10)
        allLevels = Level[];    nx = 104
        for  multiplet  in multiplets
            for  level in  multiplet.levels   push!(allLevels, level)      end
        end
        sortedLevels = Base.sort( allLevels, lt=Base.isless)
        energy0      = sortedLevels[1].energy
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "    J Parity   No. electrons    Energy " * TableStrings.inUnits("energy") * "        Leading configuration") 
        println(stream, "  ", TableStrings.hLine(nx))
        for  (n, level)  in  enumerate(sortedLevels)
            if  n > N   break   end
            mc2  = level.mc .* level.mc;   index = findmax(mc2)[2]
            conf = Basics.extractNonrelativisticConfigurationFromCsfR(level.basis.csfs[index],  level.basis)
            #
            sa   = TableStrings.flushright(12, string(LevelSymmetry(level.J, level.parity)))  * 
                   TableStrings.flushright(12, string(conf.NoElectrons))                      * 
                   TableStrings.flushright(20, @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", level.energy - energy0))) * "       " *
                   TableStrings.flushleft(64, string(conf)) 
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        
        return( nothing )
    end
   

    """
    `Basics.displayOrbitalOverlap(stream::IO, leftOrbitals::Dict{Subshell, Orbital}, rightOrbitals::Dict{Subshell, Orbital}, 
                                  grid::Radial.Grid; kappas::Array{Int64,1} = [-1,1,-2])`  
        ... display on stream the overlap of the orbitals from the set of leftOrbs and rightOrbs, and for the given
            kappa-symmetries. A neat table is printed but nothing is returned otherwise.
    """
    function Basics.displayOrbitalOverlap(stream::IO, leftOrbitals::Dict{Subshell, Orbital}, rightOrbitals::Dict{Subshell, Orbital}, 
                                          grid::Radial.Grid; kappas::Array{Int64,1} = [-1,1,-2])
        nx = 150
        println(stream, "\n  Orbital overlaps") 
        println(stream,   "  ================") 
        #
        for kappa in kappas
            sa = "\n  "
            for  (kl,vl) in leftOrbitals
                if   kl.kappa != kappa          continue    end
                for  (kr,vr) in rightOrbitals
                    if   kr.kappa != kappa      continue    end
                    overlap = RadialIntegrals.overlap(vl, vr, grid)
                    if     length(sa) < nx - 20   sa =    sa * "    |<$kl|$kr>|^2 = " * @sprintf("%.4e", overlap*overlap)
                    else   println(stream, sa);   sa =  "  " * "    |<$kl|$kr>|^2 = " * @sprintf("%.4e", overlap*overlap) 
                    end
                end
            end
            println(stream, sa)
        end
        ## println(stream, "  ", TableStrings.hLine(nx))
        
        return( nothing )
    end
   

    """
    `   + (stream::IO, channels::Array{AtomicState.GreenChannel,1}; N::Int64=10)`  
        ... display on stream the lowest N level energies of all levels from the (Green) channles in ascending order, and together with the 
            configuration of most leading term in the level expansion. A neat table is printed but nothing is returned otherwise.
    """
    function Basics.displayLevels(stream::IO, channels::Array{AtomicState.GreenChannel,1}; N::Int64=10)
        multiplets = Multiplet[]
        for  channel in channels    push!(multiplets, channel.gMultiplet)   end
        Basics.displayLevels(stream, multiplets)
        
        return( nothing )
    end


    #=====
    """
    `Basics.excludeDoubles(confList::Array{Configuration,1})`  
        ... to exlude from the (non-relativistic) confList all 'doubles', ie. configurations that occured before in the list; 
            a new confList::Array{Configuration,1} is returned.
    """
    function Basics.excludeDoubles(confList::Array{Configuration,1})
        confListNew = [ deepcopy(confList[1]) ]

        for  ic  in 2:length(confList)
            add = true
            for  id in 1:ic-1    if   confList[ic] == confList[id]    add = false;   break   end    end
            if  add   push!( confListNew, confList[ic])   end
        end    
        
        return( confListNew )
    end
        
        
    """
    `Basics.excludeDoubles(csfList::Array{CsfR,1})`  
        ... to exlude from the (relativistic) csfList all 'doubles', ie. CSF that occured already before in the list; 
            a new csfList::Array{CsfR,1} is returned.
    """
    function Basics.excludeDoubles(confList::Array{CsfR,1})
        error("Not yet implemented !")
        
        csfListNew = [ deepcopy(csfList[1]) ]
        
        for  ic  in 2:length(csfList)
            add = true
            for  id in 1:ic-1
                if   csfList[ic] == csfList[id]    add = false;   break   end
            end
            if  add   push!( csfListNew, csfList[ic])	end
        end    
        return( csfListNew )
    end
    =====#
    
        
    """
    `Basics.excludeConfigurations(confList::Array{Configuration,1},NoElectrons::Int64)`  
        ... excludes all configurations from confList with less than NoElectrons; 
            a newConfList::Array{Configuration,1} is returned. 
    """
    function Basics.excludeConfigurations(confList::Array{Configuration,1},NoElectrons::Int64)
        newConfList = Configuration[]
        for  conf in confList       
            if  conf.NoElectrons >= NoElectrons     push!(newConfList, conf)    end
        end
        
        return( newConfList )
    end
    
        
    """
    `Basics.extractLeadingConfiguration(level::ManyElectron.Level)`  
        ... extract the leading configuration of the given level; a conf::Configuration is returned.
    """
    function Basics.extractLeadingConfiguration(level::ManyElectron.Level)
        allConfs = Basics.extractNonrelativisticConfigurations(level.basis)
        weights  = zeros(length(allConfs))
        for  (ia, allConf) in  enumerate(allConfs)
            for (ic, csf) in enumerate(level.basis.csfs)
                if  allConf == Basics.extractNonrelativisticConfigurationFromCsfR(csf, level.basis)     
                    weights[ia] = weights[ia] + level.mc[ic]
                end
            end
        end
        # Determine index of maximum and return the corresponding configuration
        wx   = findmax(weights)
        conf = allConfs[ wx[2] ]
        
        return( conf )
    end
    
        
    """
    `Basics.extractLeadingConfiguration(cLevel::Cascade.Level)`  
        ... extract the leading configuration of the given level; a conf::Configuration is returned.
    """
    function Basics.extractLeadingConfiguration(cLevel::Cascade.Level)
        # First extract the right ManyElectron.Level that is to be analyzed
        level = Level()
        if length(cLevel.parents) > 0   
            parent = cLevel.parents[1]
            if      parent.process == Basics.Auger()       level = parent.lineSet.linesA[parent.index].finalLevel
            elseif  parent.process == Basics.Radiative()   level = parent.lineSet.linesR[parent.index].finalLevel
            else    error("stop a")
            end
        elseif length(cLevel.daugthers) > 0   
            daugth = cLevel.daugthers[1]
            if      daugth.process == Basics.Auger()       level = daugth.lineSet.linesA[daugth.index].initialLevel
            elseif  daugth.process == Basics.Radiative()   level = daugth.lineSet.linesR[daugth.index].initialLevel
            else    error("stop b")
            end
        else        error("stop c")
        end
    
        allConfs = Basics.extractNonrelativisticConfigurations(level.basis)
        weights  = zeros(length(allConfs))
        for  (ia, allConf) in  enumerate(allConfs)
            for (ic, csf) in enumerate(level.basis.csfs)
                if  allConf == Basics.extractNonrelativisticConfigurationFromCsfR(csf, level.basis)     
                    weights[ia] = weights[ia] + level.mc[ic]
                end
            end
        end
        # Determine index of maximum and return the corresponding configuration
        wx   = findmax(weights)
        conf = allConfs[ wx[2] ]
        
        return( conf )
    end
    
        
    """
    `Basics.extractNoOpenShells(conf::Configuration)`  
        ... determine the number of open (nonrelativistic) shells in conf; a singleton of type AbstractOpenShell is returned.
    """
    function Basics.extractNoOpenShells(conf::Configuration)
        ns = 0
        for  (shell, occ)  in conf.shells
            if      occ == 0  ||  occ == 2*(2*shell.l + 1)
            else    ns = ns + 1
            end
        end
        
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
        
        return( shellList )
    end


    """
    `Basics.extractNonrelativisticShellList(confs::Array{Configuration,1})  
        ... extract a nonrelativistic shell list from the given configuration list; the shells are not ordered. 
            A shellList::Array{Subshell,1} is returned.
    """
    function Basics.extractNonrelativisticShellList(confs::Array{Configuration,1}) 
        shellList = Shell[]
        
        for conf in confs
            for (k,v) in conf.shells
                if  k in shellList      else    push!(shellList, k)     end
            end
        end
        @show shellList
        
        return( shellList )
    end


    """
    `Basics.extractNonrelativisticShellList(multiplet::Multiplet)  
        ... extract a nonrelativistic shell list from the given multiplet; however, the shells are not ordered. 
            A shellList::Array{Shell,1} is returned.
    """
    function Basics.extractNonrelativisticShellList(multiplet::Multiplet) 
        shellList = Shell[]
        confs     = Basics.extractNonrelativisticConfigurations(multiplet.levels[1].basis)
        
        for conf in confs
            for (k,v) in conf.shells
                if  k in shellList      else    push!(shellList, k)     end
            end
        end
        @show shellList
        
        return( shellList )
    end


    """
    `Basics.extractNonrelativisticShellList(multiplets::Array{Multiplet,1})  
        ... extract a nonrelativistic shell list from the given list of multiplets; however, the shells are not ordered. 
            A shellList::Array{Shell,1} is returned.
    """
    function Basics.extractNonrelativisticShellList(multiplets::Array{Multiplet,1}) 
        shellList = Shell[]
        for  mp in multiplets
            shellList = Basics.merge(shellList, Basics.extractNonrelativisticShellList(mp))
        end
        @show shellList
        
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
                ## wa = Base.merge( wa, Dict( shell => ( (s, csfR.occupation[s], csfR.seniorityNr[s], 
                ##                                        Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ), 
                ##                                        (-9, 0, -9, 0, Basics.twice(csfR.subshellX[s]) ) ) ) )
                wa = Base.merge( wa, Dict( shell => ( (-9, 0, -9, 0, 0 ),
                                                    (s, csfR.occupation[s], csfR.seniorityNr[s], 
                                                    Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ) ) ) )
            elseif  2l + 1 == j2    continue
            elseif  2l - 1 == j2    &&    (n != nx  ||  l != lx  ||  j2x != j2 + 2)
                    error("stop a")
            elseif  0 < occ + occx < 2* (2l + 1)
                wa = Base.merge( wa, Dict(shell => ( (s, csfR.occupation[s], csfR.seniorityNr[s], 
                                                    Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ), 
                                                    (sx, csfR.occupation[sx], csfR.seniorityNr[sx], 
                                                    Basics.twice(csfR.subshellJ[sx]), Basics.twice(csfR.subshellX[sx]) ) ) ) )
            end
        end
        
        return( wa )
    end


    """
    `Basics.extractRelativisticSubshellList(confList::Array{Configuration,1})`  
        ... extract all (relativistic) subshells that (will) contribute to the given configuration list, and for which (start) 
            orbitals will be needed in course of the computation. A subshellList::Array{Subshell,1} is returned.
    """
    function Basics.extractRelativisticSubshellList(confList::Array{Configuration,1})
        # First extract all shells, then unique this list and, finally, convert in corresponding subshells
        shellList = Shell[]
        
        for  conf in confList
            for  (k,v) in conf.shells   push!( shellList, k)    end
        end
        
        shellList = unique(shellList)
        subshellList = Subshell[]
        
        for  sh in shellList    
            if  sh.l == 0    push!( subshellList, Subshell(sh.n, -1))
            else             push!( subshellList, Subshell(sh.n, sh.l));    push!( subshellList, Subshell(sh.n, -sh.l -1))
            end
        end
        subshellList = Base.sort(subshellList)
        println(">> Generated subshell list:  $subshellList ")
        return( subshellList )
    end


    """
    `Basics.extractRelativisticSubshellList(csf::CsfR, subshells::Array{Subshell,1})`  
        ... extract all (relativistic) subshells that are occupied in the given CSF with regard to the given subshells. 
            A subshellList::Array{Subshell,1} is returned.
    """
    function Basics.extractRelativisticSubshellList(csf::CsfR, subshells::Array{Subshell,1})
        reducedSubshells = Subshell[]
        
        for (i, subsh) in enumerate(subshells)
            if  csf.occupation[i]  != 0     push!(reducedSubshells, subsh)  end
        end
        
        return( reducedSubshells )
    end


    """
    `Basics.extractRelativisticSubshellList(level::Level)`  
        ... extract all (relativistic) subshells that are occupied in the given level in CSF with weigth > 0.1. 
            A subshellList::Array{Subshell,1} is returned.
    """
    function Basics.extractRelativisticSubshellList(level::Level)
        # Only consider CSF in the basis with non-zero mixing coefficients
        reducedSubshells = Subshell[]
        ##x wa = findmax(abs.(level.mc));    wb = Basics.extractRelativisticSubshellList(level.basis.csfs[ wa[2] ], level.basis.subshells)
        ##x @show wa, wa[2], wb

        for (i, coeff) in enumerate(level.mc)
            if  coeff*coeff  >  0.1   
                redSubshells     = Basics.extractRelativisticSubshellList(level.basis.csfs[i], level.basis.subshells)
                reducedSubshells = append!(reducedSubshells, redSubshells)
            end
        end
        reducedSubshells     = unique(reducedSubshells)
        
        return( reducedSubshells )
    end


    """
    `Basics.extractRelativisticSubshellList(rep::Representation)`  
        ... extract all (relativistic) subshells that (will) contribute to the given configuration list, and for which (start) 
            orbitals will be needed in course of the computation. A subshellList::Array{Subshell,1} is returned.
    """
    function Basics.extractRelativisticSubshellList(rep::Representation)
        # First extract all shells, then unique this list and, finally, convert in corresponding subshells
        shellList = Shell[]
        
        for  conf in rep.refConfigs
            for  (k,v) in conf.shells   push!( shellList, k)    end
        end
        
        for step in rep.repType.steps
            append!( shellList, step.seFrom);    append!( shellList, step.seTo)    
            append!( shellList, step.deFrom);    append!( shellList, step.deTo)   
            append!( shellList, step.teFrom);    append!( shellList, step.teTo)    
            append!( shellList, step.qeFrom);    append!( shellList, step.qeTo)   
        end
        
        shellList = unique(shellList)
        subshellList = Subshell[]
        
        for  sh in shellList    
            if  sh.l == 0    push!( subshellList, Subshell(sh.n, -1))
            else             push!( subshellList, Subshell(sh.n, sh.l));    push!( subshellList, Subshell(sh.n, -sh.l -1))
            end
        end
        
        return( subshellList )
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
        
        return( (shellList, occList) )
    end
