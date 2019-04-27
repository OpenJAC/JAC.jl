
## export  exclude

"""
`Basics.excludeDoubles()  ... exclude 'double entries' from various list; see below:

  + `(confList::Array{Configuration,1})`  ... to exlude from the (non-relativistic) confList all 'doubles', ie. configurations 
                                              that occured before in the list; a new confList::Array{Configuration,1} is returned.

  + `(csfList::Array{CsfR,1})`  ... to exlude from the (relativistic) csfList all 'doubles', ie. CSF that occured already before in the list; 
                                    a new csfList::Array{CsfR,1} is returned.
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

