
"""
`Basics.determineEnergySharings(energy::Float64, NoEnergySharings::Int64)`  
    ... to determine the NoEnergySharings sharings by using the Gauss-Legendre integration points in the interval
        [0., energy]; an Array{Tuple{Float64,Float64},1} is returned. This methods assumes that the package 
        GaussQuadrature is 'used'.
"""
function Basics.determineEnergySharings(energy::Float64, NoEnergySharings::Int64)
    xx, ww = GaussQuadrature.legendre(NoEnergySharings)
    sharings = Tuple{Float64,Float64}[]
    for  x in xx
        xrel = (x+1.) / 2
        push!(sharings, (xrel*energy, energy - xrel*energy) )
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
`Basics.determineSelectedLines(lineList::Array{Tuple{Int64,Int64},1}, initialMultiplet::Multiplet, finalMultiplet::Multiplet)`  ... to determine 
                            the specified lines as tupels of initial- and final levels. A level index is 0 in lineList always refers to
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
                               finalMultiplet::Multiplet)`  ... to determine the specified pathways as tupels of initial-, intermediate-
                               and final levels. A level index 0 in pathwayList always refers to all levels from the corresponding multiplet. 
                               A (unique) newList::Array{Tuple{Int64,Int64,Int64},1} is returned.
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



