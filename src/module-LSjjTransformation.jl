
"""
`module JAC.LSjjTransformation`  ... a submodel of JAC that contains methods and (numerical) values for performing the 
                                     LS-jj transformation; this transformation is mainly based on global data lists which are only
                                     accessible within this module-
"""
module  LSjjTransformation

    using  JAC
    
    
    """
    `struct  CsfNR`  ... defines a type for a relativistic configuration state function (CSF) in terms of its orbitals (sequence of orbitals), 
    		 their occupation as well as the angular momenta, seniority and the coupling of the corresponding antisymmetric subshell states.
    
        + useStandardSubshells ::Bool                 ... Determines whether the subshell list has to be obtained from some outer structure
                                                          given by the subshells field below. 
        + J                    ::AngularJ64           ... Total angular momentum J.
        + parity               ::Parity               ... Total parity.
        + occupation           ::Array{Int64,1}       ... occupation of the orbitals with regard to the specified orbital list.
        + seniority            ::Array{Int64,1}       ... list of seniority values of the antisymmetric subshell states
        + subshellJ            ::Array{AngularJ64,1}  ... list of J-values of the antisymmetric subshell states.
        + subshellX            ::Array{AngularJ64,1}  ... intermediate X-values in the coupling of the antisymmetric subshell states.
        + subshells            ::Array{Subshell,1}    ... Explicitly given subshell list if useStandardSubshells == false.
    """
    struct  CsfNR
        useStandardSubshells   ::Bool 
        J		               ::AngularJ64
        parity  	           ::Parity
        occupation	           ::Array{Int64,1}
        seniority	           ::Array{Int64,1}
        subshellJ	           ::Array{AngularJ64,1}
        subshellX	           ::Array{AngularJ64,1}
        subshells	           ::Array{Subshell,1}
    end
    
    
    # `Base.string(csf::CsfR)`  ... provides a String notation for csf::CsfR.
    function Base.string(csf::CsfR) 
    	  sa = "\n   CSF: "
        if    csf.useStandardSubshells

    	      for  i in 1:length(csf.occupation)
    		    sa = sa * JAC.subshellStateString(string(JAC.JAC_STANDARD_SUBSHELL_LIST[i]), csf.occupation[i], csf.seniority[i], 
                                                  csf.subshellJ[i], csf.subshellX[i])
    		    sa = sa * ", "
    	      end
    	      sa = sa * ": J=" * string(csf.J) * string(csf.parity)

            else
    	      for  i in 1:length(csf.subshells)
    		    sa = sa * JAC.subshellStateString(string(csf.subshells[i]), csf.occupation[i], csf.seniority[i], 
                                                  csf.subshellJ[i], csf.subshellX[i])
    		    sa = sa * ", "
    	      end
    	      sa = sa * ": J=" * string(csf.J) * string(csf.parity)
    	  end
    	  return( sa )
    end
    

    """
    `Base.:(==)(csfa::CsfR, csfb::CsfR)`  ... compares (recursively) two relativistic CSFs and return true if all subfields are equal, 
    						          and false otherwise
    """
    function  Base.:(==)(csfa::CsfR, csfb::CsfR)
        # Stop the comparison with an error message if the CSF are not defined w.r.t a common subshell list.
        if   !(csfa.useStandardSubshells)   ||   !(csfb.useStandardSubshells)	## ||   length(csfa.occupation) != length(csfb.occupation)
            error("To determine the equivalence of two CSF requires that both are defined on a common subshell list.")
        end
    
        if  length(csfa.occupation) != length(csfb.occupation)      return( false )    end
        if  csfa.J  !=  csfb.J  ||  csfa.parity  !=  csfb.parity	return( false )    end
    	if  csfa.occupation      !=  csfb.occupation			    return( false )    end
        if  csfa.seniority       !=  csfb.seniority			        return( false )    end
        if  csfa.subshellJ       !=  csfb.subshellJ			        return( false )    end
        if  csfa.subshellX       !=  csfb.subshellX			        return( false )    end
    	
    	  return( true )
    end


"""
  + `("CSF list: from single ConfigurationR", conf::ConfigurationR, subshellList::Array{Subshell,1}) ... to construct from a given 
                 (relativistic) configuration all possible CSF with regard to the subshell order as specified by subshellList; 
                 a list::Array{CsfR,1} is returned.
"""
function generate(sa::String, conf::ConfigurationR, subshellList::Array{Subshell,1})
    parity  = JAC.determineParity(conf)
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



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


    """
    `struct  SubshellTermLS`  ... defines a struct for a subshell term (state) in LS-coupling

        + l           ::Int32    ... angular momentum l
        + w           ::Int32    ... occupation
        + QQ          ::Int32    ... subshell total quasispin 2*Q
        + LL          ::Int32    ... subshell total angular momentum 2*L
        + SS          ::Int32    ... subshell total angular momentum 2*S
    """  
    struct  SubshellTermLS
        l             ::Int32
        w             ::Int32
        QQ            ::Int32
        LL            ::Int32
        SS            ::Int32
    end


    """
    `struct  LS_jj_qn`  ... defines a struct for the generalized quantum numbers (qn) of a LS_jj matrix element.

        + w           ::Int32    ... occupation
        + QQ          ::Int32    ... subshell total quasispin 2*Q
        + LL          ::Int32    ... subshell total angular momentum 2*L
        + SS          ::Int32    ... subshell total angular momentum 2*S
        + JJ          ::Int32    ... subshell total angular momentum 2*J
        + Nm          ::Int32    ... subshell Nm
        + Qm          ::Int32    ... subshell quantum number 2*Qm
        + Jm          ::Int32    ... subshell quantum number 2*Jm
        + Qp          ::Int32    ... subshell quantum number 2*Qp
        + Jp          ::Int32    ... subshell quantum number 2*Jp 
    """  
    struct  LS_jj_qn
        w             ::Int32
        QQ            ::Int32
        LL            ::Int32
        SS            ::Int32
        JJ            ::Int32
        Nm            ::Int32
        Qm            ::Int32
        Jm            ::Int32
        Qp            ::Int32
        Jp            ::Int32
    end


    """
    `struct  LS_jj_me`  ... defines a struct for the generalized quantum numbers (qn) of a LS_jj matrix element.

        + qn          ::LS_jj_qn    ... list of quantum numbers of this LS_jj matrix element
        + factor      ::Int32       ... integer factor
        + nom         ::Int64       ... nominator
        + denom       ::Int64       ... denominator
    """  
    struct  LS_jj_me
        qn            ::LS_jj_qn
        factor        ::Int32 
        nom           ::Int64
        denom         ::Int64
    end
    
    const   termLS_s = [ SubshellTermLS(0, 0, 0, 0, 1),    SubshellTermLS(0, 0, 1, 0, 0) ]
    const   termLS_p = [ SubshellTermLS(1, 0, 0, 0, 3),    SubshellTermLS(1, 0, 2, 2, 1),    SubshellTermLS(1, 0, 0, 4, 1), 
                         SubshellTermLS(1, 0, 3, 0, 0),    SubshellTermLS(1, 0, 1, 2, 2),    SubshellTermLS(1, 0, 1, 4, 0)  ]
    const   termLS_d = [ SubshellTermLS(2, 0, 0, 0, 5),    SubshellTermLS(2, 0, 0, 0, 1),    SubshellTermLS(2, 0, 2, 2, 3), 
                         SubshellTermLS(2, 0, 2, 2, 1),    SubshellTermLS(2, 0, 4, 4, 1),    SubshellTermLS(2, 0, 2, 4, 1),  
                         SubshellTermLS(2, 0, 0, 4, 3),    SubshellTermLS(2, 0, 0, 4, 1),    SubshellTermLS(2, 0, 2, 6, 3), 
                         SubshellTermLS(2, 0, 2, 6, 1),    SubshellTermLS(2, 0, 0, 6, 1),    SubshellTermLS(2, 0, 2, 8, 1),   
                         SubshellTermLS(2, 0, 0, 8, 3),    SubshellTermLS(2, 0, 0, 8, 1),    SubshellTermLS(2, 0, 2,10, 1), 
                         SubshellTermLS(2, 0, 0,12, 1),    SubshellTermLS(2, 0, 5, 0, 0),    SubshellTermLS(2, 0, 1, 0, 0),  
                         SubshellTermLS(2, 0, 3, 2, 2),    SubshellTermLS(2, 0, 1, 2, 2),    SubshellTermLS(2, 0, 1, 4, 4),   
                         SubshellTermLS(2, 0, 1, 4, 2),    SubshellTermLS(2, 0, 3, 4, 0),    SubshellTermLS(2, 0, 1, 4, 0), 
                         SubshellTermLS(2, 0, 3, 6, 2),    SubshellTermLS(2, 0, 1, 6, 2),    SubshellTermLS(2, 0, 1, 6, 0),   
                         SubshellTermLS(2, 0, 1, 8, 2),    SubshellTermLS(2, 0, 3, 8, 0),    SubshellTermLS(2, 0, 1, 8, 0),     
                         SubshellTermLS(2, 0, 1,10, 2),    SubshellTermLS(2, 0, 1,12, 0) ]


    """
    `JAC.LSjjTransformation.getLSjjCoefficient(l::Int64, N::Int64, qn::LS_jj_qn)`  
        ... to return the value of the LS-jj transformation matrix element for a given set of quantum numbers.
            Note that all (generalized) angular momentum quantum numbers except of l must be given twice the original numbers, 
            i.e. for the quantum numbers Q, L, S, J, Qm, Jm, Qp, Jp.
    """
    function getLSjjCoefficient(l::Int64, N::Int64, qn::LS_jj_qn)
        global  LS_jj_p_3,  LS_jj_p_4,  LS_jj_p_5,  LS_jj_p_6,  LS_jj_d_3,  LS_jj_d_4,  LS_jj_d_5,  LS_jj_d_6,  LS_jj_d_7,
                LS_jj_d_8,  LS_jj_d_9,  LS_jj_d_10
        #
        wa = 0.
        if      lshell == 0
        elseif  lshell == 1
            if      N == 3   ## Use data from the array LS_jj_p_3
                for  me  in  LS_jj_p_3    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 4   ## Use data from the array LS_jj_p_4
                for  me  in  LS_jj_p_4    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 5   ## Use data from the array LS_jj_p_5
                for  me  in  LS_jj_p_5    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 6   ## Use data from the array LS_jj_p_6
                for  me  in  LS_jj_p_6    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            else    error("stop a")
            end
        elseif  lshell == 2
            if      N == 3   ## Use data from the array LS_jj_d_3
                for  me  in  LS_jj_d_3    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 4   ## Use data from the array LS_jj_d_4
                for  me  in  LS_jj_d_4    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 5   ## Use data from the array LS_jj_d_5
                for  me  in  LS_jj_d_5    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 6   ## Use data from the array LS_jj_d_6
                for  me  in  LS_jj_d_6    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 7   ## Use data from the array LS_jj_d_7
                for  me  in  LS_jj_d_7    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 8   ## Use data from the array LS_jj_d_8
                for  me  in  LS_jj_d_8    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 9   ## Use data from the array LS_jj_d_9
                for  me  in  LS_jj_d_9    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 10  ## Use data from the array LS_jj_d_10
                for  me  in  LS_jj_d_10   if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            else    error("stop b")
            end
        else        error("stop e")
        end
        
        return(wa)
    end

    ## p^3  shell;  11 ME in total
    const   LS_jj_p_3 = [ LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 1, 0, 1, 2, 0),  1, 1, 1),
                          LS_jj_me( LS_jj_qn(0, 0, 0, 3, 3, 0, 1, 0, 1, 3), -1, 2, 9),
                          LS_jj_me( LS_jj_qn(0, 0, 0, 3, 3, 1, 0, 1, 0, 4),  1, 5, 9),
                          LS_jj_me( LS_jj_qn(0, 0, 0, 3, 3, 2, 1, 0, 1, 3),  1, 2, 9),
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 0, 1, 0, 1, 3),  1, 1, 2),
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 1, 0, 1, 3),  1, 1, 2),
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 0, 1, 0, 1, 3),  1, 5,18),
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 1, 0, 1, 0, 4),  1, 4, 9),
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 2, 1, 0, 1, 3), -1, 5,18),
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 1, 0, 1, 0, 4),  1, 1, 1)  ]
    ## p^4  shell;  9 ME in total
    const   LS_jj_p_4 = [ LS_jj_me( LS_jj_qn(0, 3, 0, 0, 0, 0, 1, 0, 2, 0),  1, 1, 3),
                          LS_jj_me( LS_jj_qn(0, 3, 0, 0, 0, 2, 1, 0, 2, 0),  1, 2, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 0, 1, 0, 2, 0), -1, 2, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 1, 0, 2, 0),  1, 1, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 1, 0, 1, 1, 3),  1, 1, 1),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 1, 0, 1, 1, 3),  1, 1, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 1, 0, 0, 4),  1, 2, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 1, 0, 1, 1, 3), -1, 2, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 1, 0, 0, 4),  1, 1, 3)  ]
    ## p^5  shell;  2 ME in total
    const   LS_jj_p_5 = [ LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 1, 0, 1, 2, 0),  1, 1, 1),
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 1, 0, 1, 3),  1, 1, 1)  ]
    ## p^6  shell;  1 ME in total
    const   LS_jj_p_6 = [ LS_jj_me( LS_jj_qn(0, 3, 0, 0, 0, 2, 1, 0, 2, 0),  1, 1, 1)  ]

    
    ## d^3  shell;  73 ME in total
    const   LS_jj_d_3 = [ LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 1, 1, 3, 1, 4), -1,  7, 15), ## 1
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 2, 0, 4, 2, 5), -1,  8, 15), ## 2
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 1, 1, 3, 1, 4), -1,  8, 15), ## 3
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 2, 0, 4, 2, 5),  1,  7, 15), ## 4
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 0, 2, 0, 0, 3), -1, 21,125), ## 5
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 3, 0),  1,  8,375), ## 6
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 1, 4), -1, 28, 75), ## 7
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 2, 0, 4, 2, 5), -1, 28, 75), ## 8
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 3, 0), -1,  8,125), ## 9
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 0, 2, 0, 0, 3), -1, 12, 25), ## 10
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 3, 0), -1,  7,150), ## 11
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 1, 4), -1,  1, 15), ## 12
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 0, 4, 2, 5),  1,  4, 15), ## 13
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 3, 0),  1,  7, 50), ## 14
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 1, 1, 3, 3, 0),  1,  3,  4), ## 16
        		          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 3, 1, 3, 3, 0),  1,  1,  4), ## 19
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 0, 2, 0, 0, 3), -1,  8, 25), ## 20
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 3, 0),  1,  7,100), ## 21
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 1, 4),  1,  2,  5), ## 22
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 3, 0), -1, 21,100), ## 24
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 0, 2, 0, 0, 3),  1,  4,125), ## 25
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 3, 0),  1, 14,125), ## 26
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 1, 4), -1,  4, 25), ## 27
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 2, 0, 4, 2, 5),  1,  9, 25), ## 28
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 3, 0), -1, 42,125), ## 29
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 0, 2, 0, 2, 5),  1, 24,125), ## 30
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 1, 1, 3, 1, 4), -1,  1, 25), ## 31
        		          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 1, 1, 3, 1, 8), -1, 72,125), ## 32
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 2, 2, 0, 2, 5), -1, 24,125), ## 33
                		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 0, 2, 0, 2, 5),  1,  1,  2), ## 35
        		          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 2, 2, 0, 2, 5),  1,  1,  2), ## 38
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 0, 2, 0, 2, 5), -1,  7,150), ## 40
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 1, 1, 3, 1, 4),  1, 16, 35), ## 41
        		          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 1, 1, 3, 1, 8), -1, 32,175), ## 42
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 2, 0, 2, 5),  1,  7,150), ## 43
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 0, 4, 2, 5),  1,  4, 15), ## 44
        		          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 0, 2, 0, 2, 5),  1, 28,375), ## 45
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 1, 1, 3, 1, 4), -1,  8,175), ## 46
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 1, 1, 3, 1, 8),  1,121,875), ## 47
        		          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 2, 0, 2, 5), -1, 28,375), ## 48
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 0, 4, 2, 5),  1,  2,  3), ## 49
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 0, 2, 0, 2, 5),  1, 14, 75), ## 50
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 1, 1, 3, 1, 4),  1, 16, 35), ## 51
        		          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 1, 1, 3, 1, 8),  1, 18,175), ## 52
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 2, 0, 2, 5), -1, 14, 75), ## 53
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 0, 4, 2, 5), -1,  1, 15), ## 54
        		          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 1, 1, 3, 1, 4), -1,  3, 35), ## 55
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 1, 1, 3, 1, 8),  1,  5,  7), ## 56
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 2, 0, 4, 2, 5),  1,  1,  5), ## 57
        		          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 1, 1, 3, 1, 4),  1,121,280), ## 58
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 1, 1, 3, 1, 8),  1, 15, 56), ## 59
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 2, 0, 4, 2, 5), -1,  3, 10), ## 60
        		          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 1, 1, 3, 1, 4), -1, 27, 56), ## 61
                		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 1, 1, 3, 1, 8),  1,  1, 56), ## 62
        	        	  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 2, 0, 4, 2, 5), -1,  1,  2), ## 63
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 0, 2, 0, 0, 9), -1, 12, 25), ## 64
        	        	  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 1, 1, 3, 1, 8),  1, 11, 25), ## 65
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 2, 0, 4, 2, 5),  1,  2, 25), ## 66
                		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 0, 2, 0, 0, 9), -1, 54,125), ## 67
        		          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 1, 1, 3, 1, 8), -1, 22,125), ## 68
                		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 2, 0, 4, 2, 5), -1, 49,125), ## 69
                		  LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 0, 2, 0, 0, 9), -1, 11,125), ## 70
        		          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 1, 1, 3, 1, 8), -1, 48,125), ## 71
                		  LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 2, 0, 4, 2, 5),  1, 66,125), ## 72
        	        	  LS_jj_me( LS_jj_qn(0, 2,10, 1,11, 1, 1, 3, 1, 8), -1,  1,  1) ]## 73
    ## d^4  shell;  198 ME in total 					    	       
    const   LS_jj_d_4 = [ LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 0, 2, 0, 3, 0),  1,  3,  10), ## 1
        		  LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 2, 2, 0, 3, 0),  1,  3,   5), ## 3
        		  LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 4, 2, 0, 3, 0),  1,  1,  10), ## 5
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 0, 2, 0, 3, 0), -1,  7, 250), ## 6
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 1, 1, 3, 0, 3), -1, 64, 125), ## 7
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 2, 2, 0, 3, 0),  1,  7, 125), ## 8
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 2, 0, 4, 1, 4),  1,  8,  25), ## 9
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 4, 2, 0, 3, 0), -1, 21, 250), ## 10
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 0, 2, 0, 3, 0), -1,  8,  15), ## 11
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 2, 2, 0, 3, 0),  1,  1,  15), ## 13
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 4, 2, 0, 3, 0),  1,  2,   5), ## 15
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 0, 2, 0, 3, 0), -1, 28, 375), ## 16
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 1, 1, 3, 0, 3),  1, 54, 125), ## 17
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 2, 0, 3, 0),  1, 56, 375), ## 18
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 0, 4, 1, 4),  1,  3,  25), ## 19
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 4, 2, 0, 3, 0), -1, 28, 125), ## 20
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 0, 2, 0, 3, 0), -1,  8, 125), ## 21
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 1, 1, 3, 0, 3), -1,  7, 125), ## 22
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 2, 2, 0, 3, 0),  1, 16, 125), ## 23
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 2, 0, 4, 1, 4), -1, 14,  25), ## 24
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 4, 2, 0, 3, 0), -1, 24, 125), ## 25
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 1, 1, 3, 2, 5),  1,  2,   3), ## 26
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 3, 1, 3, 2, 5),  1,  1,   3), ## 29
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 1, 1, 3, 2, 5),  1,  7,  75), ## 30
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 1, 1, 3, 0, 3),  1,  9,  25), ## 31
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 2, 0, 4, 1, 4),  1,  9,  25), ## 32
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 3, 1, 3, 2, 5), -1, 14,  75), ## 33
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 1, 1, 3, 2, 5),  1,  6,  25), ## 34
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 1, 1, 3, 0, 3), -1,  7,  50), ## 35
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 2, 0, 4, 1, 4), -1,  7,  50), ## 36
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 3, 1, 3, 2, 5), -1, 12,  25), ## 37
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 2, 1, 1, 3, 0, 3),  1,  1,   2), ## 39
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 2, 2, 0, 4, 1, 4), -1,  1,   2), ## 40
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 0, 2, 0, 1, 4),  1, 56, 375), ## 42
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 1, 1, 3, 2, 5),  1, 42, 125), ## 43
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 2, 0, 1, 4),  1,112, 375), ## 45
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 0, 4, 3, 0), -1,  6, 125), ## 46
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 3, 1, 3, 2, 5),  1, 21, 125), ## 49
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 0, 2, 0, 1, 4), -1, 64, 375), ## 50
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 1, 1, 3, 2, 5), -1,  3, 125), ## 51
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 1, 1, 3, 0, 3),  1, 27, 125), ## 52
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 2, 0, 1, 4),  1, 32, 375), ## 53
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 0, 4, 1, 4),  1,  3,  35), ## 55
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 0, 4, 1, 8),  1,324, 875), ## 56
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 3, 1, 3, 2, 5),  1,  6, 125), ## 57
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 0, 2, 0, 1, 4),  1,  2,  25), ## 58
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 1, 1, 3, 2, 5),  1,  2,  25), ## 59
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 1, 1, 3, 0, 3), -1,  1,  50), ## 60
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 2, 0, 1, 4), -1,  1,  25), ## 61
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 0, 4, 1, 4), -1,  1,  14), ## 63
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 0, 4, 1, 8),  1, 96, 175), ## 64
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 3, 1, 3, 2, 5), -1,  4,  25), ## 65
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 0, 2, 0, 1, 4), -1,  6,  25), ## 66
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 1, 1, 3, 2, 5),  1,  8,  75), ## 67
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 1, 1, 3, 0, 3),  1,  3,  50), ## 68
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 2, 0, 1, 4),  1,  3,  25), ## 69
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 0, 4, 1, 4), -1,  3,  14), ## 71
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 0, 4, 1, 8), -1,  8, 175), ## 72
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 3, 1, 3, 2, 5), -1, 16,  75), ## 73
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 0, 2, 0, 1, 4),  1,  4,  25), ## 74
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 1, 1, 3, 2, 5), -1,  4,  25), ## 75
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 2, 0, 1, 4),  1,  8,  25), ## 77
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 0, 4, 3, 0),  1,  7,  25), ## 78
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 3, 1, 3, 2, 5), -1,  2,  25), ## 81
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 0, 2, 0, 1, 4), -1,  2,  25), ## 82
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 1, 1, 3, 2, 5), -1,  2,  25), ## 83
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 1, 1, 3, 0, 3), -1,  8,  25), ## 84
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 2, 0, 1, 4),  1,  1,  25), ## 85
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 0, 4, 1, 4), -1,  2,   7), ## 87
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 0, 4, 1, 8),  1,  6, 175), ## 88
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 3, 1, 3, 2, 5),  1,  4,  25), ## 89
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 0, 2, 0, 1, 4), -1,  3, 125), ## 90
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 1, 1, 3, 2, 5),  1, 64, 375), ## 91
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 2, 0, 1, 4), -1,  6, 125), ## 93
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 0, 4, 3, 0),  1, 84, 125), ## 94
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 3, 1, 3, 2, 5),  1, 32, 375), ## 97
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 0, 2, 0, 1, 4), -1, 12, 125), ## 98
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 1, 1, 3, 2, 5),  1, 16, 375), ## 99
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 1, 1, 3, 0, 3), -1, 48, 125), ## 100
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 2, 0, 1, 4),  1,  6, 125), ## 101
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 0, 4, 1, 4),  1, 12,  35), ## 103
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 0, 4, 1, 8),  1,  1, 875), ## 104
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 3, 1, 3, 2, 5), -1, 32, 375), ## 105
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 1, 1, 3, 2, 5),  1,  1,  25), ## 106
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 1, 1, 3, 0, 3),  1,  2, 175), ## 107
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 1, 1, 3, 0, 9), -1,  3,   7), ## 108
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 2, 0, 4, 1, 4), -1,  2, 175), ## 109
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 2, 0, 4, 1, 8),  1,  3,   7), ## 110
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 3, 1, 3, 2, 5), -1,  2,  25), ## 111
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 1, 1, 3, 2, 5),  1,  2,  21), ## 112
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 1, 1, 3, 0, 3),  1,  3,  49), ## 113
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 1, 1, 3, 0, 9),  1, 18,  49), ## 114
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 2, 0, 4, 1, 4), -1, 12,  49), ## 115
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 2, 0, 4, 1, 8),  1,  2,  49), ## 116
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 3, 1, 3, 2, 5), -1,  4,  21), ## 117
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 1, 1, 3, 2, 5),  1,  2,   3), ## 118
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 3, 1, 3, 2, 5),  1,  1,   3), ## 123
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 1, 1, 3, 2, 5),  1,  1, 600), ## 124
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 1, 1, 3, 0, 3), -1, 48, 175), ## 125
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 1, 1, 3, 0, 9),  1,  9,  56), ## 126
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 2, 0, 4, 1, 4),  1, 48, 175), ## 127
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 2, 0, 4, 1, 8),  1,  2,   7), ## 128
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 3, 1, 3, 2, 5), -1,  1, 300), ## 129
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 1, 1, 3, 2, 5), -1,  1,  10), ## 130
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 1, 1, 3, 0, 3),  1, 16,  35), ## 131
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 1, 1, 3, 0, 9),  1,  3,  70), ## 132
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 2, 0, 4, 1, 4),  1,  1,  35), ## 133
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 2, 0, 4, 1, 8),  1,  6,  35), ## 134
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 3, 1, 3, 2, 5),  1,  1,   5), ## 135
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 1, 1, 3, 2, 5),  1, 27, 280), ## 136
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 1, 1, 3, 0, 3),  1, 48, 245), ## 137
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 1, 1, 3, 0, 9), -1,  1,1960), ## 138
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 2, 0, 4, 1, 4),  1,108, 245), ## 139
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 2, 0, 4, 1, 8), -1, 18, 245), ## 140
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 3, 1, 3, 2, 5), -1, 27, 140), ## 141
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 0, 2, 0, 1, 8),  1, 24, 125), ## 142
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 1, 1, 3, 2, 5),  1,  3, 125), ## 143
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 1, 1, 3, 0, 9), -1, 11,  25), ## 144
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 2, 0, 1, 8), -1, 12, 125), ## 145
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 0, 4, 1, 4),  1,  2, 175), ## 146
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 0, 4, 1, 8),  1, 33, 175), ## 147
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 3, 1, 3, 2, 5), -1,  6, 125), ## 148
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 0, 2, 0, 1, 8),  1,  4,  15), ## 149
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 1, 1, 3, 2, 5),  1,  2,  15), ## 150
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 2, 2, 0, 1, 8),  1,  8,  15), ## 152
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 3, 1, 3, 2, 5),  1,  1,  15), ## 155
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 0, 2, 0, 1, 8), -1, 64, 375), ## 156
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 1, 1, 3, 2, 5),  1,289,3000), ## 157
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 1, 1, 3, 0, 9),  1, 11, 200), ## 158
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 2, 0, 1, 8),  1, 32, 375), ## 159
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 0, 4, 1, 4),  1,  4, 175), ## 160
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 0, 4, 1, 8),  1, 66, 175), ## 161
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 3, 1, 3, 2, 5), -1,289,1500), ## 162
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 0, 2, 0, 1, 8),  1, 16, 125), ## 163
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 1, 1, 3, 2, 5),  1, 81,1000), ## 164
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 1, 1, 3, 0, 9),  1, 33, 200), ## 165
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 2, 0, 1, 8), -1,  8, 125), ## 166
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 0, 4, 1, 4),  1, 48, 175), ## 167
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 0, 4, 1, 8), -1, 22, 175), ## 168
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 3, 1, 3, 2, 5), -1, 81, 500), ## 169
        		  LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 0, 2, 0, 1, 8),  1,  1,  15), ## 170
        		  LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 1, 1, 3, 2, 5), -1,  8,  15), ## 171
        		  LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 2, 2, 0, 1, 8),  1,  2,  15), ## 173
        		  LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 3, 1, 3, 2, 5), -1,  4,  15), ## 176
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 0, 2, 0, 1, 8), -1, 44, 375), ## 177
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 1, 1, 3, 2, 5), -1, 11, 750), ## 178
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 1, 1, 3, 0, 9), -1,  9,  50), ## 179
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 2, 0, 1, 8),  1, 22, 375), ## 180
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 0, 4, 1, 4),  1, 99, 175), ## 181
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 0, 4, 1, 8), -1,  6, 175), ## 182
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 3, 1, 3, 2, 5),  1, 11, 375), ## 183
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 0, 2, 0, 1, 8),  1, 22, 375), ## 184
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 1, 1, 3, 2, 5), -1, 44, 375), ## 185
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 1, 1, 3, 0, 9),  1,  4,  25), ## 186
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 2, 0, 1, 8), -1, 11, 375), ## 187
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 0, 4, 1, 4),  1, 22, 175), ## 188
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 0, 4, 1, 8),  1, 48, 175), ## 189
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 3, 1, 3, 2, 5),  1, 88, 375), ## 190
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2,10, 1, 1, 3, 0, 9),  1, 16,  25), ## 191
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2,10, 2, 0, 4, 1, 8), -1,  9,  25), ## 192
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2,10, 1, 1, 3, 0, 9),  1,  9,  25), ## 193
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2,10, 2, 0, 4, 1, 8),  1, 16,  25), ## 194
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2,12, 1, 1, 3, 0, 9),  1,  3,   5), ## 195
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2,12, 2, 0, 4, 1, 8),  1,  2,   5), ## 196
        		  LS_jj_me( LS_jj_qn(0, 1,12, 0,12, 1, 1, 3, 0, 9),  1,  2,   5), ## 197
        		  LS_jj_me( LS_jj_qn(0, 1,12, 0,12, 2, 0, 4, 1, 8), -1,  3,   5) ]## 198
    ## d^5  shell;  249 ME in total					    		
    const   LS_jj_d_5 = [ LS_jj_me( LS_jj_qn(0, 0, 0, 1, 1, 1, 1, 3, 1, 4),  1,   2,	5), ## 1
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 1, 1, 2, 0, 4, 0, 3),  1,   1,	5), ## 3
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 1, 1, 3, 1, 3, 1, 4), -1,   2,	5), ## 4
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 1, 1, 3, 1, 4), -1,   7,   30), ## 5
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 2, 0, 4, 2, 5), -1,   8,   15), ## 6
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 3, 1, 3, 1, 4), -1,   7,   30), ## 8
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 1, 1, 3, 1, 4), -1,   4,   15), ## 9
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 2, 0, 4, 2, 5),  1,   7,   15), ## 10
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 3, 1, 3, 1, 4), -1,   4,   15), ## 12
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 1, 1, 1, 3, 1, 4),  1,   1,   10), ## 13
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 1, 2, 0, 4, 0, 3), -1,   4,	5), ## 15
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 1, 3, 1, 3, 1, 4), -1,   1,   10), ## 16
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 3, 0),  1,  16,  375), ## 17
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 1, 4), -1,  14,   75), ## 18
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 2, 2, 0, 0, 3), -1,  21,  125), ## 19
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 2, 0, 4, 2, 5), -1,  28,   75), ## 20
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 3, 0), -1,  16,  375), ## 22
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 1, 4), -1,  14,   75), ## 23
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 3, 0), -1,   7,   75), ## 24
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 1, 4), -1,   1,   30), ## 25
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 2, 0, 0, 3), -1,  12,   25), ## 26
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 0, 4, 2, 5),  1,   4,   15), ## 27
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 3, 0),  1,   7,   75), ## 29
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 1, 4), -1,   1,   30), ## 30
        		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 1, 1, 3, 3, 0),  1,   1,	2), ## 31
        		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 3, 1, 3, 3, 0),  1,   1,	2), ## 36
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 3, 0),  1,   7,   50), ## 38
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 1, 4),  1,   1,	5), ## 39
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 2, 2, 0, 0, 3), -1,   8,   25), ## 40
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 3, 0), -1,   7,   50), ## 43
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 1, 4),  1,   1,	5), ## 44
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 3, 1, 1, 3, 1, 4),  1,   2,	5), ## 46
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 3, 2, 0, 4, 0, 3), -1,   1,	5), ## 49
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 3, 3, 1, 3, 1, 4), -1,   2,	5), ## 51
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 1, 1, 3, 1, 4), -1,   1,   10), ## 53
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 2, 0, 4, 0, 3), -1,   4,	5), ## 56
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 3, 1, 3, 1, 4),  1,   1,   10), ## 58
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 3, 0),  1,  28,  125), ## 59
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 1, 4), -1,   2,   25), ## 60
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 2, 2, 0, 0, 3),  1,   4,  125), ## 61
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 2, 0, 4, 2, 5),  1,   9,   25), ## 62
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 3, 0), -1,  28,  125), ## 64
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 1, 4), -1,   2,   25), ## 65
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 0, 2, 0, 2, 5), -1,  24,  625), ## 66
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 1, 1, 3, 1, 4),  1,   2,  125), ## 67
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 1, 1, 3, 1, 8),  1, 144,  625), ## 68
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 2, 2, 0, 2, 5),  1,  24,  625), ## 69
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 2, 0, 4, 0, 3), -1,   1,  125), ## 71
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 2, 0, 4, 0, 9), -1,  48,  125), ## 72
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 3, 1, 3, 1, 4), -1,   2,  125), ## 73
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 3, 1, 3, 1, 8), -1, 144,  625), ## 74
        		  LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 4, 2, 0, 2, 5), -1,  24,  625), ## 75
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 0, 2, 0, 2, 5),  1,  24,  125), ## 76
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 1, 1, 3, 1, 4), -1,   1,   50), ## 77
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 1, 1, 3, 1, 8), -1,  36,  125), ## 78
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 3, 1, 3, 1, 4), -1,   1,   50), ## 83
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 3, 1, 3, 1, 8), -1,  36,  125), ## 84
        		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 4, 2, 0, 2, 5), -1,  24,  125), ## 85
        		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 0, 2, 0, 2, 5),  1,   1,	6), ## 86
        		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 2, 2, 0, 2, 5),  1,   2,	3), ## 89
        		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 4, 2, 0, 2, 5),  1,   1,	6), ## 95
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 0, 2, 0, 2, 5), -1,   7,  150), ## 96
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 1, 1, 3, 1, 4),  1,   8,   35), ## 97
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 1, 1, 3, 1, 8), -1,  16,  175), ## 98
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 0, 4, 2, 5),  1,   4,   15), ## 100
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 3, 1, 3, 1, 4),  1,   8,   35), ## 103
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 3, 1, 3, 1, 8), -1,  16,  175), ## 104
        		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 4, 2, 0, 2, 5),  1,   7,  150), ## 105
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 0, 2, 0, 2, 5), -1,   8,   75), ## 106
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 1, 1, 3, 1, 4),  1,  81,  490), ## 107
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 1, 1, 3, 1, 8), -1,   4, 1225), ## 108
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 2, 2, 0, 2, 5),  1,   8,   75), ## 109
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 2, 0, 4, 0, 3), -1,  36,  245), ## 111
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 2, 0, 4, 0, 9),  1,  48,  245), ## 112
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 3, 1, 3, 1, 4), -1,  81,  490), ## 113
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 3, 1, 3, 1, 8),  1,   4, 1225), ## 114
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 4, 2, 0, 2, 5), -1,   8,   75), ## 115
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 0, 2, 0, 2, 5),  1,   7,   75), ## 116
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 1, 1, 3, 1, 4),  1,   4,   35), ## 117
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 1, 1, 3, 1, 8), -1,   8,  175), ## 118
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 2, 2, 0, 2, 5), -1,   7,   75), ## 119
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 2, 0, 4, 0, 3), -1,   8,   35), ## 121
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 2, 0, 4, 0, 9), -1,   6,   35), ## 122
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 3, 1, 3, 1, 4), -1,   4,   35), ## 123
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 3, 1, 3, 1, 8),  1,   8,  175), ## 124
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 4, 2, 0, 2, 5),  1,   7,   75), ## 125
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 0, 2, 0, 2, 5),  1,  28,  375), ## 126
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 1, 1, 3, 1, 4), -1,   4,  175), ## 127
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 1, 1, 3, 1, 8),  1, 121, 1750), ## 128
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 0, 4, 2, 5),  1,   2,	3), ## 130
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 3, 1, 3, 1, 4), -1,   4,  175), ## 133
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 3, 1, 3, 1, 8),  1, 121, 1750), ## 134
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 4, 2, 0, 2, 5), -1,  28,  375), ## 135
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 0, 2, 0, 2, 5),  1,  14,   75), ## 136
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 1, 1, 3, 1, 4),  1,   8,   35), ## 137
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 1, 1, 3, 1, 8),  1,   9,  175), ## 138
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 0, 4, 2, 5), -1,   1,   15), ## 140
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 3, 1, 3, 1, 4),  1,   8,   35), ## 143
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 3, 1, 3, 1, 8),  1,   9,  175), ## 144
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 4, 2, 0, 2, 5), -1,  14,   75), ## 145
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 0, 2, 0, 2, 5),  1,  14,  375), ## 146
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 1, 1, 3, 1, 4), -1,   8,  175), ## 147
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 1, 1, 3, 1, 8),  1, 121,  875), ## 148
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 2, 2, 0, 2, 5), -1,  14,  375), ## 149
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 2, 0, 4, 0, 3), -1,  64,  175), ## 151
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 2, 0, 4, 0, 9),  1,  27,  175), ## 152
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 3, 1, 3, 1, 4),  1,   8,  175), ## 153
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 3, 1, 3, 1, 8), -1, 121,  875), ## 154
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 4, 2, 0, 2, 5),  1,  14,  375), ## 155
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 0, 2, 0, 2, 5), -1,  36,  625), ## 156
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 1, 1, 3, 1, 4), -1, 972, 6125), ## 157
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 1, 1, 3, 1, 8), -1,5043,61250), ## 158
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 2, 2, 0, 2, 5),  1,  36,  625), ## 159
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 2, 0, 4, 0, 3), -1,1536, 6125), ## 161
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 2, 0, 4, 0, 9), -1, 578, 6125), ## 162
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 3, 1, 3, 1, 4),  1, 972, 6125), ## 163
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 3, 1, 3, 1, 8),  1,5043,61250), ## 164
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 4, 2, 0, 2, 5), -1,  36,  625), ## 165
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 1, 1, 3, 1, 4), -1,   2,  245), ## 166
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 1, 1, 3, 1, 8), -1,  54,  245), ## 167
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 2, 0, 4, 0, 3),  1,   1,  245), ## 169
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 2, 0, 4, 0, 9),  1, 132,  245), ## 170
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 3, 1, 3, 1, 4),  1,   2,  245), ## 171
        		  LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 3, 1, 3, 1, 8),  1,  54,  245), ## 172
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 1, 1, 3, 1, 4), -1,   3,   70), ## 173
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 1, 1, 3, 1, 8),  1,   5,   14), ## 174
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 2, 0, 4, 2, 5),  1,   1,	5), ## 175
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 3, 1, 3, 1, 4), -1,   3,   70), ## 178
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 3, 1, 3, 1, 8),  1,   5,   14), ## 179
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 1, 1, 3, 1, 4),  1, 121,  560), ## 180
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 1, 1, 3, 1, 8),  1,  15,  112), ## 181
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 2, 0, 4, 2, 5), -1,   3,   10), ## 182
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 3, 1, 3, 1, 4),  1, 121,  560), ## 185
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 3, 1, 3, 1, 8),  1,  15,  112), ## 186
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 1, 1, 3, 1, 4), -1,  25,  112), ## 187
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 1, 1, 3, 1, 8),  1, 243, 2800), ## 188
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 2, 0, 4, 0, 3), -1,   2,	7), ## 190
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 2, 0, 4, 0, 9),  1,  33,  350), ## 191
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 3, 1, 3, 1, 4),  1,  25,  112), ## 192
        		  LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 3, 1, 3, 1, 8), -1, 243, 2800), ## 193
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 1, 1, 3, 1, 4), -1,  27,  112), ## 194
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 1, 1, 3, 1, 8),  1,   1,  112), ## 195
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 2, 0, 4, 2, 5), -1,   1,	2), ## 196
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 3, 1, 3, 1, 4), -1,  27,  112), ## 199
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 3, 1, 3, 1, 8),  1,   1,  112), ## 200
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 1, 1, 3, 1, 4), -1,   9,   98), ## 201
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 1, 1, 3, 1, 8), -1,1369, 7350), ## 202
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 2, 0, 4, 0, 3), -1,   4,   49), ## 204
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 2, 0, 4, 0, 9), -1,1331, 3675), ## 205
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 3, 1, 3, 1, 4),  1,   9,   98), ## 206
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 3, 1, 3, 1, 8),  1,1369, 7350), ## 207
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 1, 1, 3, 1, 4), -1,  99,  560), ## 208
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 1, 1, 3, 1, 8),  1,  11, 1680), ## 209
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 2, 0, 4, 0, 3),  1,  22,   35), ## 211
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 2, 0, 4, 0, 9), -1,   1,  210), ## 212
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 3, 1, 3, 1, 4),  1,  99,  560), ## 213
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 3, 1, 3, 1, 8), -1,  11, 1680), ## 214
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 1, 1, 3, 1, 8),  1,  11,   50), ## 215
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 2, 2, 0, 0, 9), -1,  12,   25), ## 216
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 2, 0, 4, 2, 5),  1,   2,   25), ## 217
        		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 3, 1, 3, 1, 8),  1,  11,   50), ## 219
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 1, 1, 3, 1, 8), -1,  11,  125), ## 220
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 2, 2, 0, 0, 9), -1,  54,  125), ## 221
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 2, 0, 4, 2, 5), -1,  49,  125), ## 222
        		  LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 3, 1, 3, 1, 8), -1,  11,  125), ## 224
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 9, 1, 1, 3, 1, 8), -1,   1,	6), ## 225
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 9, 2, 0, 4, 0, 9), -1,   2,	3), ## 228
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3, 9, 3, 1, 3, 1, 8),  1,   1,	6), ## 229
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 9, 1, 1, 3, 1, 8),  1,   1,	3), ## 230
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 9, 2, 0, 4, 0, 9), -1,   1,	3), ## 233
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 1, 9, 3, 1, 3, 1, 8), -1,   1,	3), ## 234
        		  LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 1, 1, 3, 1, 8), -1,  24,  125), ## 235
        		  LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 2, 2, 0, 0, 9), -1,  11,  125), ## 236
        		  LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 2, 0, 4, 2, 5),  1,  66,  125), ## 237
        		  LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 3, 1, 3, 1, 8), -1,  24,  125), ## 239
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3,11, 1, 1, 3, 1, 8), -1,   6,   25), ## 240
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3,11, 2, 0, 4, 0, 9), -1,  13,   25), ## 241
        		  LS_jj_me( LS_jj_qn(0, 0, 8, 3,11, 3, 1, 3, 1, 8),  1,   6,   25), ## 242
        		  LS_jj_me( LS_jj_qn(0, 2,10, 1,11, 1, 1, 3, 1, 8), -1,   1,	2), ## 243
        		  LS_jj_me( LS_jj_qn(0, 2,10, 1,11, 3, 1, 3, 1, 8), -1,   1,	2), ## 245
        		  LS_jj_me( LS_jj_qn(0, 0,12, 1,11, 1, 1, 3, 1, 8),  1,  13,   50), ## 246
        		  LS_jj_me( LS_jj_qn(0, 0,12, 1,11, 2, 0, 4, 0, 9), -1,  12,   25), ## 247
        		  LS_jj_me( LS_jj_qn(0, 0,12, 1,11, 3, 1, 3, 1, 8), -1,  13,   50), ## 248
        		  LS_jj_me( LS_jj_qn(0, 0,12, 1,13, 2, 0, 4, 0, 9), -1,   1,	1) ]## 249
    ## d^6  shell;  198 ME in total					    		  
    const   LS_jj_d_6 = [ LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 0, 2, 0, 3, 0),  1,  1,  10), ## 1
        		  LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 2, 2, 0, 3, 0),  1,  3,   5), ## 2
        		  LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 4, 2, 0, 3, 0),  1,  3,  10), ## 5
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 0, 2, 0, 3, 0), -1, 21, 250), ## 6
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 2, 2, 0, 3, 0),  1,  7, 125), ## 7
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 2, 0, 4, 1, 4),  1,  8,  25), ## 8
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 3, 1, 3, 0, 3), -1, 64, 125), ## 9
        		  LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 4, 2, 0, 3, 0), -1,  7, 250), ## 10
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 0, 2, 0, 3, 0), -1,  2,   5), ## 11
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 2, 2, 0, 3, 0), -1,  1,  15), ## 12
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 4, 2, 0, 3, 0),  1,  8,  15), ## 15
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 0, 2, 0, 3, 0), -1, 28, 125), ## 16
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 2, 0, 3, 0),  1, 56, 375), ## 17
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 0, 4, 1, 4),  1,  3,  25), ## 18
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 3, 1, 3, 0, 3),  1, 54, 125), ## 19
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 4, 2, 0, 3, 0), -1, 28, 375), ## 20
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 0, 2, 0, 3, 0), -1, 24, 125), ## 21
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 2, 2, 0, 3, 0),  1, 16, 125), ## 22
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 2, 0, 4, 1, 4), -1, 14,  25), ## 23
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 3, 1, 3, 0, 3), -1,  7, 125), ## 24
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 4, 2, 0, 3, 0), -1,  8, 125), ## 25
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 1, 1, 3, 2, 5),  1,  1,   3), ## 26
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 3, 1, 3, 2, 5),  1,  2,   3), ## 28
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 1, 1, 3, 2, 5),  1, 14,  75), ## 30
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 2, 0, 4, 1, 4),  1,  9,  25), ## 31
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 3, 1, 3, 2, 5), -1,  7,  75), ## 32
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 3, 1, 3, 0, 3),  1,  9,  25), ## 33
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 1, 1, 3, 2, 5),  1, 12,  25), ## 34
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 2, 0, 4, 1, 4), -1,  7,  50), ## 35
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 3, 1, 3, 2, 5), -1,  6,  25), ## 36
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 3, 1, 3, 0, 3), -1,  7,  50), ## 37
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 2, 2, 0, 4, 1, 4), -1,  1,   2), ## 39
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 2, 3, 1, 3, 0, 3),  1,  1,   2), ## 41
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 1, 1, 3, 2, 5),  1, 21, 125), ## 42
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 2, 0, 1, 4),  1,112, 375), ## 43
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 0, 4, 3, 0), -1,  6, 125), ## 44
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 3, 1, 3, 2, 5),  1, 42, 125), ## 47
        		  LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 4, 2, 0, 1, 4),  1, 56, 375), ## 49
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 1, 1, 3, 2, 5), -1,  6, 125), ## 50
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 2, 0, 1, 4), -1, 32, 375), ## 51
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 0, 4, 1, 4),  1,  3,  35), ## 53
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 0, 4, 1, 8),  1,324, 875), ## 54
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 3, 1, 3, 2, 5),  1,  3, 125), ## 55
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 3, 1, 3, 0, 3),  1, 27, 125), ## 56
        		  LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 4, 2, 0, 1, 4),  1, 64, 375), ## 57
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 1, 1, 3, 2, 5),  1,  4,  25), ## 58
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 2, 0, 1, 4),  1,  1,  25), ## 59
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 0, 4, 1, 4), -1,  1,  14), ## 61
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 0, 4, 1, 8),  1, 96, 175), ## 62
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 3, 1, 3, 2, 5), -1,  2,  25), ## 63
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 3, 1, 3, 0, 3), -1,  1,  50), ## 64
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 4, 2, 0, 1, 4), -1,  2,  25), ## 65
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 1, 1, 3, 2, 5),  1, 16,  75), ## 66
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 2, 0, 1, 4), -1,  3,  25), ## 67
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 0, 4, 1, 4), -1,  3,  14), ## 69
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 0, 4, 1, 8), -1,  8, 175), ## 70
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 3, 1, 3, 2, 5), -1,  8,  75), ## 71
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 3, 1, 3, 0, 3),  1,  3,  50), ## 72
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 4, 2, 0, 1, 4),  1,  6,  25), ## 73
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 1, 1, 3, 2, 5), -1,  2,  25), ## 74
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 2, 0, 1, 4),  1,  8,  25), ## 75
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 0, 4, 3, 0),  1,  7,  25), ## 76
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 3, 1, 3, 2, 5), -1,  4,  25), ## 79
        		  LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 4, 2, 0, 1, 4),  1,  4,  25), ## 81
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 1, 1, 3, 2, 5), -1,  4,  25), ## 82
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 2, 0, 1, 4), -1,  1,  25), ## 83
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 0, 4, 1, 4), -1,  2,   7), ## 85
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 0, 4, 1, 8),  1,  6, 175), ## 86
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 3, 1, 3, 2, 5),  1,  2,  25), ## 87
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 3, 1, 3, 0, 3), -1,  8,  25), ## 88
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 4, 2, 0, 1, 4),  1,  2,  25), ## 89
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 1, 1, 3, 2, 5),  1, 32, 375), ## 90
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 2, 0, 1, 4), -1,  6, 125), ## 91
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 0, 4, 3, 0),  1, 84, 125), ## 92
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 3, 1, 3, 2, 5),  1, 64, 375), ## 95
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 4, 2, 0, 1, 4), -1,  3, 125), ## 97
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 1, 1, 3, 2, 5),  1, 32, 375), ## 98
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 2, 0, 1, 4), -1,  6, 125), ## 99
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 0, 4, 1, 4),  1, 12,  35), ## 101
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 0, 4, 1, 8),  1,  1, 875), ## 102
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 3, 1, 3, 2, 5), -1, 16, 375), ## 103
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 3, 1, 3, 0, 3), -1, 48, 125), ## 104
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 4, 2, 0, 1, 4),  1, 12, 125), ## 105
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 1, 1, 3, 2, 5),  1,  2,  25), ## 106
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 2, 0, 4, 1, 4), -1,  2, 175), ## 107
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 2, 0, 4, 1, 8),  1,  3,   7), ## 108
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 3, 1, 3, 2, 5), -1,  1,  25), ## 109
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 3, 1, 3, 0, 3),  1,  2, 175), ## 110
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 3, 1, 3, 0, 9), -1,  3,   7), ## 111
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 1, 1, 3, 2, 5),  1,  4,  21), ## 112
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 2, 0, 4, 1, 4), -1, 12,  49), ## 113
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 2, 0, 4, 1, 8),  1,  2,  49), ## 114
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 3, 1, 3, 2, 5), -1,  2,  21), ## 115
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 3, 1, 3, 0, 3),  1,  3,  49), ## 116
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 3, 1, 3, 0, 9),  1, 18,  49), ## 117
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 1, 1, 3, 2, 5),  1,  1,   3), ## 118
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 3, 1, 3, 2, 5),  1,  2,   3), ## 121
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 1, 1, 3, 2, 5),  1,  1, 300), ## 124
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 2, 0, 4, 1, 4),  1, 48, 175), ## 125
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 2, 0, 4, 1, 8),  1,  2,   7), ## 126
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 3, 1, 3, 2, 5), -1,  1, 600), ## 127
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 3, 1, 3, 0, 3), -1, 48, 175), ## 128
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 3, 1, 3, 0, 9),  1,  9,  56), ## 129
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 1, 1, 3, 2, 5), -1,  1,   5), ## 130
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 2, 0, 4, 1, 4),  1,  1,  35), ## 131
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 2, 0, 4, 1, 8),  1,  6,  35), ## 132
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 3, 1, 3, 2, 5),  1,  1,  10), ## 133
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 3, 1, 3, 0, 3),  1, 16,  35), ## 134
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 3, 1, 3, 0, 9),  1,  3,  70), ## 135
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 1, 1, 3, 2, 5),  1, 27, 140), ## 136
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 2, 0, 4, 1, 4),  1,108, 245), ## 137
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 2, 0, 4, 1, 8), -1, 18, 245), ## 138
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 3, 1, 3, 2, 5), -1, 27, 280), ## 139
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 3, 1, 3, 0, 3),  1, 48, 245), ## 140
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 3, 1, 3, 0, 9), -1,  1,1960), ## 141
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 1, 1, 3, 2, 5),  1,  6, 125), ## 142
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 2, 0, 1, 8),  1, 12, 125), ## 143
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 0, 4, 1, 4),  1,  2, 175), ## 144
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 0, 4, 1, 8),  1, 33, 175), ## 145
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 3, 1, 3, 2, 5), -1,  3, 125), ## 146
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 3, 1, 3, 0, 9), -1, 11,  25), ## 147
        		  LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 4, 2, 0, 1, 8), -1, 24, 125), ## 148
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 1, 1, 3, 2, 5),  1,  1,  15), ## 149
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 2, 2, 0, 1, 8),  1,  8,  15), ## 150
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 3, 1, 3, 2, 5),  1,  2,  15), ## 153
        		  LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 4, 2, 0, 1, 8),  1,  4,  15), ## 155
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 1, 1, 3, 2, 5),  1,289,1500), ## 156
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 2, 0, 1, 8), -1, 32, 375), ## 157
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 0, 4, 1, 4),  1,  4, 175), ## 158
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 0, 4, 1, 8),  1, 66, 175), ## 159
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 3, 1, 3, 2, 5), -1,289,3000), ## 160
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 3, 1, 3, 0, 9),  1, 11, 200), ## 161
        		  LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 4, 2, 0, 1, 8),  1, 64, 375), ## 162
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 1, 1, 3, 2, 5),  1, 81, 500), ## 163
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 2, 0, 1, 8),  1,  8, 125), ## 164
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 0, 4, 1, 4),  1, 48, 175), ## 165
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 0, 4, 1, 8), -1, 22, 175), ## 166
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 3, 1, 3, 2, 5), -1, 81,1000), ## 167
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 3, 1, 3, 0, 9),  1, 33, 200), ## 168
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 4, 2, 0, 1, 8), -1, 16, 125), ## 169
        		  LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 1, 1, 3, 2, 5), -1,  4,  15), ## 170
        		  LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 2, 2, 0, 1, 8),  1,  2,  15), ## 171
        		  LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 3, 1, 3, 2, 5), -1,  8,  15), ## 174
        		  LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 4, 2, 0, 1, 8),  1,  1,  15), ## 176
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 1, 1, 3, 2, 5), -1, 11, 375), ## 177
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 2, 0, 1, 8), -1, 22, 375), ## 178
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 0, 4, 1, 4),  1, 99, 175), ## 179
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 0, 4, 1, 8), -1,  6, 175), ## 180
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 3, 1, 3, 2, 5),  1, 11, 750), ## 181
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 3, 1, 3, 0, 9), -1,  9,  50), ## 182
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 4, 2, 0, 1, 8),  1, 44, 375), ## 183
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 1, 1, 3, 2, 5), -1, 88, 375), ## 184
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 2, 0, 1, 8),  1, 11, 375), ## 185
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 0, 4, 1, 4),  1, 22, 175), ## 186
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 0, 4, 1, 8),  1, 48, 175), ## 187
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 3, 1, 3, 2, 5),  1, 44, 375), ## 188
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 3, 1, 3, 0, 9),  1,  4,  25), ## 189
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 4, 2, 0, 1, 8), -1, 22, 375), ## 190
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2,10, 2, 0, 4, 1, 8), -1,  9,  25), ## 191
        		  LS_jj_me( LS_jj_qn(0, 1, 8, 2,10, 3, 1, 3, 0, 9),  1, 16,  25), ## 192
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2,10, 2, 0, 4, 1, 8),  1, 16,  25), ## 193
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2,10, 3, 1, 3, 0, 9),  1,  9,  25), ## 194
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2,12, 2, 0, 4, 1, 8),  1,  2,   5), ## 195
        		  LS_jj_me( LS_jj_qn(0, 1,10, 2,12, 3, 1, 3, 0, 9),  1,  3,   5), ## 196
        		  LS_jj_me( LS_jj_qn(0, 1,12, 0,12, 2, 0, 4, 1, 8), -1,  3,   5), ## 197
        		  LS_jj_me( LS_jj_qn(0, 1,12, 0,12, 3, 1, 3, 0, 9),  1,  2,   5) ]## 198
    ## d^7  shell;  73 ME in total					    		
    const   LS_jj_d_7 = [ LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 2, 0, 4, 2, 5), -1,  8, 15), ## 1
        		          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 3, 1, 3, 1, 4), -1,  7, 15), ## 2
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 2, 0, 4, 2, 5),  1,  7, 15), ## 3
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 3, 1, 3, 1, 4), -1,  8, 15), ## 4
        		          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 3, 0),  1,  8,125), ## 5
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 2, 0, 4, 2, 5), -1, 28, 75), ## 6
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 3, 0), -1,  8,375), ## 7
        		          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 1, 4), -1, 28, 75), ## 8
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 4, 2, 0, 0, 3), -1, 21,125), ## 9
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 3, 0), -1,  7, 50), ## 10
        		          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 0, 4, 2, 5),  1,  4, 15), ## 11
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 3, 0),  1,  7,150), ## 12
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 1, 4), -1,  1, 15), ## 13
        		          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 4, 2, 0, 0, 3), -1, 12, 25), ## 14
                		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 1, 1, 3, 3, 0),  1,  1,  4), ## 15
                		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 3, 1, 3, 3, 0),  1,  3,  4), ## 17
        		          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 3, 0),  1, 21,100), ## 20
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 3, 0), -1,  7,100), ## 22
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 1, 4),  1,  2,  5), ## 23
        		          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 4, 2, 0, 0, 3), -1,  8, 25), ## 24
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 3, 0),  1, 42,125), ## 25
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 2, 0, 4, 2, 5),  1,  9, 25), ## 26
        		          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 3, 0), -1, 14,125), ## 27
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 1, 4), -1,  4, 25), ## 28
                		  LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 4, 2, 0, 0, 3),  1,  4,125), ## 29
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 2, 2, 0, 2, 5),  1, 24,125), ## 30
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 3, 1, 3, 1, 4), -1,  1, 25), ## 32
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 3, 1, 3, 1, 8), -1, 72,125), ## 33
                		  LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 4, 2, 0, 2, 5), -1, 24,125), ## 34
                		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 2, 2, 0, 2, 5),  1,  1,  2), ## 35
                		  LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 4, 2, 0, 2, 5),  1,  1,  2), ## 39
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 2, 0, 2, 5), -1,  7,150), ## 40
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 0, 4, 2, 5),  1,  4, 15), ## 41
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 3, 1, 3, 1, 4),  1, 16, 35), ## 42
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 3, 1, 3, 1, 8), -1, 32,175), ## 43
                		  LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 4, 2, 0, 2, 5),  1,  7,150), ## 44
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 2, 0, 2, 5),  1, 28,375), ## 45
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 0, 4, 2, 5),  1,  2,  3), ## 46
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 3, 1, 3, 1, 4), -1,  8,175), ## 47
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 3, 1, 3, 1, 8),  1,121,875), ## 48
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 4, 2, 0, 2, 5), -1, 28,375), ## 49
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 2, 0, 2, 5),  1, 14, 75), ## 50
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 0, 4, 2, 5), -1,  1, 15), ## 51
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 3, 1, 3, 1, 4),  1, 16, 35), ## 52
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 3, 1, 3, 1, 8),  1, 18,175), ## 53
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 4, 2, 0, 2, 5), -1, 14, 75), ## 54
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 2, 0, 4, 2, 5),  1,  1,  5), ## 55
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 3, 1, 3, 1, 4), -1,  3, 35), ## 56
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 3, 1, 3, 1, 8),  1,  5,  7), ## 57
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 2, 0, 4, 2, 5), -1,  3, 10), ## 58
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 3, 1, 3, 1, 4),  1,121,280), ## 59
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 3, 1, 3, 1, 8),  1, 15, 56), ## 60
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 2, 0, 4, 2, 5), -1,  1,  2), ## 61
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 3, 1, 3, 1, 4), -1, 27, 56), ## 62
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 3, 1, 3, 1, 8),  1,  1, 56), ## 63
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 2, 0, 4, 2, 5),  1,  2, 25), ## 64
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 3, 1, 3, 1, 8),  1, 11, 25), ## 65
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 4, 2, 0, 0, 9), -1, 12, 25), ## 66
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 2, 0, 4, 2, 5), -1, 49,125), ## 67
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 3, 1, 3, 1, 8), -1, 22,125), ## 68
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 4, 2, 0, 0, 9), -1, 54,125), ## 69
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 2, 0, 4, 2, 5),  1, 66,125), ## 70
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 3, 1, 3, 1, 8), -1, 48,125), ## 71
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 4, 2, 0, 0, 9), -1, 11,125), ## 72
                          LS_jj_me( LS_jj_qn(0, 2,10, 1,11, 3, 1, 3, 1, 8), -1,  1,  1) ]## 73
    ## d^8  shell;  19 ME in total					 	       
    const   LS_jj_d_8 = [ LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 2, 2, 0, 3, 0),  1, 2,  5), ## 1
                          LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 4, 2, 0, 3, 0),  1, 3,  5), ## 2
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 2, 2, 0, 3, 0), -1, 3,  5), ## 3
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 4, 2, 0, 3, 0),  1, 2,  5), ## 4
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 3, 1, 3, 2, 5),  1, 1,  1), ## 5
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 0, 4, 3, 0), -1, 6,125), ## 6
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 3, 1, 3, 2, 5),  1,63,125), ## 7
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 4, 2, 0, 1, 4),  1,56,125), ## 8
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 0, 4, 3, 0),  1, 7, 25), ## 9
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 3, 1, 3, 2, 5), -1, 6, 25), ## 10
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 4, 2, 0, 1, 4),  1,12, 25), ## 11
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 0, 4, 3, 0),  1,84,125), ## 12
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 3, 1, 3, 2, 5),  1,32,125), ## 13
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 4, 2, 0, 1, 4), -1, 9,125), ## 14
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 3, 1, 3, 2, 5),  1, 1,  1), ## 15
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 3, 1, 3, 2, 5),  1, 1,  5), ## 16
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 4, 2, 0, 1, 8),  1, 4,  5), ## 17
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 3, 1, 3, 2, 5), -1, 4,  5), ## 18
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 4, 2, 0, 1, 8),  1, 1,  5) ]## 19
    ## d^9  shell;  2 ME in total 					    	      
    const   LS_jj_d_9 = [ LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 3, 1, 3, 3, 0),  1, 1, 1),  ## 1
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 4, 2, 0, 2, 5),  1, 1, 1) ] ## 2
    ## d^10  shell;  1 ME in total  LS_jj_qn(				    	      
    const   LS_jj_d_10= [ LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 4, 2, 0, 3, 0),  1, 1, 1) ] ## 1
    									    
end # module			
