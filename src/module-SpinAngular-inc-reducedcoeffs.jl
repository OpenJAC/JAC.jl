
    
    """
    `SpinAngular.completlyReducedCfp(j::AngularJ64, Nr::Int64, Q::AngularJ64, J::AngularJ64, Nrp::Int64, Qp::AngularJ64, Jp::AngularJ64)`  
        ... returns the (completely) reduced coefficient of fractional parentage (j alpha Q J||| a^(qj)||| j alpha' Q' J')
            for j = 1/2, 3/2, 5/2 and 7/2, and where alpha == Nr refers to some additional quantum number for j=9/2 subshell
            states. In all other cases Nr = Nrp = 0 is expected by this routine. A coeff::Float64 is returned
    """
    function completlyReducedCfp(j::AngularJ64, Nr::Int64, Q::AngularJ64, J::AngularJ64, Nrp::Int64, Qp::AngularJ64, Jp::AngularJ64)
        #
        # Define a tuple for comparison in the subsequent lists; apply symmetries to reduced the list of coefficients
        println("!!! Here, apply symmetries to reduced the list of coefficients !!!")
        wc = (Nr, twice(Q), twice(J), Nrp, twice(Qp), twice(Jp))
        #
        # Distinguish due to different j-values
        if       j == AngularJ64(1//2)
            if      wc == (0, 0, 0,  0, 0, 0)        wa = ( 0,     0,     0)             # !! Please, check
            end
        elseif   j == AngularJ64(3//2)
            if      wc == (0, 2, 0,  0, 1, 3)        wa = (-1,    12,     1)
            elseif  wc == (0, 0, 4,  0, 1, 3)        wa = ( 1,    20,     1)
            end
        elseif   j == AngularJ64(5//2)
            if      wc == (0, 2, 5,  0, 3, 0)        wa = ( 0,     0,     0)
            elseif  wc == (0, 2, 5,  0, 1, 4)        wa = ( 0,     0,     0)
            elseif  wc == (0, 2, 5,  0, 1, 8)        wa = ( 0,     0,     0)
            elseif  wc == (0, 0, 3,  0, 1, 4)        wa = ( 0,     0,     0)
            elseif  wc == (0, 0, 3,  0, 1, 8)        wa = ( 0,     0,     0)
            elseif  wc == (0, 0, 9,  0, 1, 4)        wa = ( 0,     0,     0)
            elseif  wc == (0, 0, 9,  0, 1, 8)        wa = ( 0,     0,     0)
            end
        ##  elseif   j == AngularJ64(5//2)  For j = 5/2, 7/2 and 9/2, we should generate the corresponding code from running
        ##  elseif   j == AngularJ64(7//2)  the routine of ANCO and printing the results in a proper format, I guess
        ##  elseif   j == AngularJ64(9//2)
        else     error("stop a")
        end
        
        wb = wa[1] * sqrt(wa[2]/wa[3])
        return( wb )
    end

    
    """
    `SpinAngular.completelyReducedWkk(kq::Int64, kj::Int64, j::AngularJ64, 
                                      Nr::Int64, Q::AngularJ64, J::AngularJ64, Nrp::Int64, Qp::AngularJ64, Jp::AngularJ64)`  
        ... returns the (completely) reduced coefficient of fractional parentage (j alpha Q J||| a^(qj)||| j alpha' Q' J')
            for j = 1/2, 3/2, 5/2 and 7/2, and where alpha == Nr refers to some additional quantum number for j=9/2 subshell
            states. In all other cases Nr = Nrp = 0 is expected by this routine. A coeff::Float64 is returned
    """
    function completelyReducedWkk(kq::Int64, kj::Int64, j::AngularJ64, 
                                  Nr::Int64, Q::AngularJ64, J::AngularJ64, Nrp::Int64, Qp::AngularJ64, Jp::AngularJ64)
        #
        # Define a tuple for comparison in the subsequent lists; apply symmetries to reduced the list of coefficients
        println("!!! Here, apply symmetries to reduced the list of coefficients !!!")
        wc = (Nr, twice(Q), twice(J), Nrp, twice(Qp), twice(Jp))
        #
        # Distinguish due to different j-values
        if       j == AngularJ64(1//2)   &&     (kq,kj) == (0,0)
            if      wc == (0, 0, 0,  0, 0, 0)           wa = ( 0,     0,     0)             # !! Please, check
            end
        elseif   j == AngularJ64(3//2)   &&     (kq,kj) == (0,0)
            if      wc == (0, 2, 0,  0, 1, 3)           wa = (-1,    12,     1)
            elseif  wc == (0, 0, 4,  0, 1, 3)           wa = ( 1,    20,     1)
            end
        elseif   j == AngularJ64(3//2)   &&     (kq,kj) == (1,0)
            if      wc == (0, 2, 0,  0, 1, 3)           wa = (-1,    12,     1)
            elseif  wc == (0, 0, 4,  0, 1, 3)           wa = ( 1,    20,     1)
            end
        ##  elseif   j == AngularJ64(5//2)  For j = 3/2, 5/2, 7/2 and 9/2, again, we should generate the corresponding code from running
        ##  elseif   j == AngularJ64(7//2)  the routine of ANCO and printing the results in a proper format, I guess
        ##  elseif   j == AngularJ64(9//2)
        else     error("stop a")
        end
        
        wb = wa[1] * sqrt(wa[2]/wa[3])
        return( wb )
    end

    
    """
    `SpinAngular.subshellTerms(j::AngularJ64)`  
        ... returns a list of subshell terms for the given j; a termList::Array{SpinAngular.SubshellTerm,1} is returned.
    """
    function  subshellTerms(j::AngularJ64)
        # First define some short-cuts to simplify the subsequent input
        h0  = AngularJ64(0//2);   h1  = AngularJ64(1//2);   h2  = AngularJ64(2//2);   h3  = AngularJ64(3//2);   h4  = AngularJ64(4//2)    
        h5  = AngularJ64(5//2);   h6  = AngularJ64(6//2);   h7  = AngularJ64(7//2);   h8  = AngularJ64(8//2);   h9  = AngularJ64(9//2)    
        h10 = AngularJ64(10//2);  h11 = AngularJ64(11//2);  h12 = AngularJ64(12//2);  h13 = AngularJ64(13//2);  h14 = AngularJ64(14//2)
        h15 = AngularJ64(15//2);  h16 = AngularJ64(16//2);  h17 = AngularJ64(17//2);  h18 = AngularJ64(18//2);  h19 = AngularJ64(19//2)
        h20 = AngularJ64(20//2);  h21 = AngularJ64(21//2);  h22 = AngularJ64(22//2);  h23 = AngularJ64(23//2);  h24 = AngularJ64(24//2)
        h25 = AngularJ64(25//2);  h26 = AngularJ64(26//2);  h27 = AngularJ64(27//2)
        #
        # Distinguish due to different j-values
        if       j == AngularJ64(1//2)
            wa = [ SubshellTerm(h1, h0, 1, h1, 0), SubshellTerm(h1, h1, 0, h0, 0) ]
            #
        elseif   j == AngularJ64(3//2)
            wa = [ SubshellTerm(h3, h1, 1, h3, 0), SubshellTerm(h3, h2, 0, h0, 0), SubshellTerm(h3, h0, 2, h4, 0) ]
            #
        elseif   j == AngularJ64(5//2)
            wa = [ SubshellTerm(h5, h2, 1, h5, 0), SubshellTerm(h5, h0, 3, h3, 0), SubshellTerm(h5, h0, 3, h9, 0),
                   SubshellTerm(h5, h3, 0, h0, 0), SubshellTerm(h5, h1, 2, h4, 0), SubshellTerm(h5, h1, 2, h8, 0) ]
            #
        elseif   j == AngularJ64(7//2)
            wa = [ SubshellTerm(h7, h1, 3, h3, 0), SubshellTerm(h7, h1, 3, h5, 0), SubshellTerm(h7, h3, 1, h7, 0),
                   SubshellTerm(h7, h1, 3, h9, 0), SubshellTerm(h7, h1, 3,h11, 0), SubshellTerm(h7, h1, 3,h15, 0),
                   SubshellTerm(h7, h4, 0, h0, 0), SubshellTerm(h7, h2, 2, h4, 0), SubshellTerm(h7, h0, 4, h8, 0),
                   SubshellTerm(h7, h0, 4, h4, 0), SubshellTerm(h7, h2, 2, h8, 0), SubshellTerm(h7, h0, 4,h10, 0),
                   SubshellTerm(h7, h2, 2,h12, 0), SubshellTerm(h7, h0, 4,h16, 0) ]
            #
        elseif   j == AngularJ64(9//2)
            wa = [ SubshellTerm(h9, h0, 5, h1, 0), SubshellTerm(h9, h2, 3, h3, 0), SubshellTerm(h9, h2, 3, h5, 0), 
                   SubshellTerm(h9, h0, 5, h5, 0), SubshellTerm(h9, h2, 3, h7, 0), SubshellTerm(h9, h0, 5, h7, 0), 
                   SubshellTerm(h9, h4, 1, h9, 0), SubshellTerm(h9, h2, 3, h9, 0), SubshellTerm(h9, h0, 5, h9, 0), 
                   SubshellTerm(h9, h2, 3,h11, 0), SubshellTerm(h9, h0, 5,h11, 0), SubshellTerm(h9, h2, 3,h13, 0),
                   SubshellTerm(h9, h0, 5,h13, 0), SubshellTerm(h9, h2, 3,h15, 0), SubshellTerm(h9, h0, 5,h15, 0), 
                   SubshellTerm(h9, h2, 3,h17, 0), SubshellTerm(h9, h0, 5,h17, 0), SubshellTerm(h9, h0, 5,h19, 0), 
                   SubshellTerm(h9, h2, 3,h21, 0), SubshellTerm(h9, h0, 5,h25, 0), SubshellTerm(h9, h5, 0, h0, 0), 
                   SubshellTerm(h9, h1, 4, h0, 0), SubshellTerm(h9, h3, 2, h4, 0), SubshellTerm(h9, h1, 4, h4, 0),  
                   SubshellTerm(h9, h1, 4, h6, 0), SubshellTerm(h9, h3, 2, h8, 0), SubshellTerm(h9, h1, 4, h8, 1), 
                   SubshellTerm(h9, h1, 4, h8, 2), SubshellTerm(h9, h1, 4,h10, 0), SubshellTerm(h9, h3, 2,h12, 0), 
                   SubshellTerm(h9, h1, 4,h12, 1), SubshellTerm(h9, h1, 4,h12, 2), SubshellTerm(h9, h1, 4,h14, 0), 
                   SubshellTerm(h9, h3, 2,h16, 0), SubshellTerm(h9, h1, 4,h16, 0), SubshellTerm(h9, h1, 4,h18, 0),  
                   SubshellTerm(h9, h1, 4,h20, 0), SubshellTerm(h9, h1, 4,h24, 0) ]
        else     error("stop a")
        end
        
        return(wa)
    end
