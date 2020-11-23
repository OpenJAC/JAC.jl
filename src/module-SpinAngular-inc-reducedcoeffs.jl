
    
    """
    `SpinAngular.completlyReducedCfpByIndices(aTerm::QspaceTerm, bTerm::QspaceTerm)`  
        ... returns the (completely) reduced coefficient of fractional parentage (j alpha Q J||| a^(qj)||| j alpha' Q' J')
            for j = 1/2, 3/2, 5/2 and 7/2, and where alpha == Nr refers to some additional quantum number for j=9/2 subshell
            states. In all other cases Nr = Nrp = 0 is expected by this routine. A coeff::Float64 is returned
    """
    function completlyReducedCfpByIndices(aTerm::QspaceTerm, bTerm::QspaceTerm)
        #
        # Distinguish due to different j-values
        if       aTerm.j != bTerm.j               error("stop a")
        elseif   aTerm.j == AngularJ64(1//2)      wa = reducedCfp1half[aTerm.index, bTerm.index]
        elseif   aTerm.j == AngularJ64(3//2)      wa = reducedCfp3half[aTerm.index, bTerm.index]
        elseif   aTerm.j == AngularJ64(5//2)      wa = reducedCfp5half[aTerm.index, bTerm.index]
        elseif   aTerm.j == AngularJ64(7//2)      wa = reducedCfp7half[aTerm.index, bTerm.index]
        elseif   aTerm.j == AngularJ64(9//2)      wa = reducedCfp9half[aTerm.index, bTerm.index]
        else     error("stop b")
        end
        @show wa
        
        wb = wa[1] * sqrt(wa[2]/wa[3])
        return( wb )
    end

    
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
    `SpinAngular.qspaceTerms(j::AngularJ64)`  
        ... returns a list of Q-space terms for the given j; a termList::Array{SpinAngular.QspaceTerm,1} is returned.
    """
    function  qspaceTerms(j::AngularJ64)
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
            wa = [ QspaceTerm(h1, h0, h1, 0,  1), QspaceTerm(h1, h1, h0, 0,  2) ]
            #
        elseif   j == AngularJ64(3//2)
            wa = [ QspaceTerm(h3, h1, h3, 0,  1), QspaceTerm(h3, h2, h0, 0,  2), QspaceTerm(h3, h0, h4, 0,  3) ]
            #
        elseif   j == AngularJ64(5//2)
            wa = [ QspaceTerm(h5, h2, h5, 0,  1), QspaceTerm(h5, h0, h3, 0,  2), QspaceTerm(h5, h0, h9, 0,  3),
                   QspaceTerm(h5, h3, h0, 0,  4), QspaceTerm(h5, h1, h4, 0,  5), QspaceTerm(h5, h1, h8, 0,  6) ]
            #
        elseif   j == AngularJ64(7//2)
            wa = [ QspaceTerm(h7, h1, h3, 0,  1), QspaceTerm(h7, h1, h5, 0,  2), QspaceTerm(h7, h3, h7, 0,  3),
                   QspaceTerm(h7, h1, h9, 0,  4), QspaceTerm(h7, h1,h11, 0,  5), QspaceTerm(h7, h1,h15, 0,  6),
                   QspaceTerm(h7, h4, h0, 0,  7), QspaceTerm(h7, h2, h4, 0,  8), QspaceTerm(h7, h0, h8, 0,  9),
                   QspaceTerm(h7, h0, h4, 0, 10), QspaceTerm(h7, h2, h8, 0, 11), QspaceTerm(h7, h0,h10, 0, 12),
                   QspaceTerm(h7, h2,h12, 0, 13), QspaceTerm(h7, h0,h16, 0, 14) ]
            #
        elseif   j == AngularJ64(9//2)
            wa = [ QspaceTerm(h9, h0, h1, 0,  1), QspaceTerm(h9, h2, h3, 0,  2), QspaceTerm(h9, h2, h5, 0,  3), 
                   QspaceTerm(h9, h0, h5, 0,  4), QspaceTerm(h9, h2, h7, 0,  5), QspaceTerm(h9, h0, h7, 0,  6), 
                   QspaceTerm(h9, h4, h9, 0,  7), QspaceTerm(h9, h2, h9, 0,  8), QspaceTerm(h9, h0, h9, 0,  9), 
                   QspaceTerm(h9, h2,h11, 0, 10), QspaceTerm(h9, h0,h11, 0, 11), QspaceTerm(h9, h2,h13, 0, 12),
                   QspaceTerm(h9, h0,h13, 0, 13), QspaceTerm(h9, h2,h15, 0, 14), QspaceTerm(h9, h0,h15, 0, 15), 
                   QspaceTerm(h9, h2,h17, 0, 16), QspaceTerm(h9, h0,h17, 0, 17), QspaceTerm(h9, h0,h19, 0, 18), 
                   QspaceTerm(h9, h2,h21, 0, 19), QspaceTerm(h9, h0,h25, 0, 20), QspaceTerm(h9, h5, h0, 0, 21), 
                   QspaceTerm(h9, h1, h0, 0, 22), QspaceTerm(h9, h3, h4, 0, 23), QspaceTerm(h9, h1, h4, 0, 24),  
                   QspaceTerm(h9, h1, h6, 0, 25), QspaceTerm(h9, h3, h8, 0, 26), QspaceTerm(h9, h1, h8, 1, 27), 
                   QspaceTerm(h9, h1, h8, 2, 28), QspaceTerm(h9, h1,h10, 0, 29), QspaceTerm(h9, h3,h12, 0, 30), 
                   QspaceTerm(h9, h1,h12, 1, 31), QspaceTerm(h9, h1,h12, 2, 32), QspaceTerm(h9, h1,h14, 0, 33), 
                   QspaceTerm(h9, h3,h16, 0, 34), QspaceTerm(h9, h1,h16, 0, 35), QspaceTerm(h9, h1,h18, 0, 36),  
                   QspaceTerm(h9, h1,h20, 0, 37), QspaceTerm(h9, h1,h24, 0, 38) ]
        else     error("stop a")
        end
        
        return( wa )
    end

    
    """
    `SpinAngular.subshellTerms(j::AngularJ64, N::Int64)`  
        ... returns a list of subshell terms for the given j^N; a termList::Array{SpinAngular.SubshellTerm,1} is returned.
    """
    function  subshellTerms(j::AngularJ64, N::Int64)
        # First define some short-cuts to simplify the subsequent input
        h0  = AngularJ64(0//2);   h1  = AngularJ64(1//2);   h2  = AngularJ64(2//2);   h3  = AngularJ64(3//2);   h4  = AngularJ64(4//2)    
        h5  = AngularJ64(5//2);   h6  = AngularJ64(6//2);   h7  = AngularJ64(7//2);   h8  = AngularJ64(8//2);   h9  = AngularJ64(9//2)    
        h10 = AngularJ64(10//2);  h11 = AngularJ64(11//2);  h12 = AngularJ64(12//2);  h13 = AngularJ64(13//2);  h14 = AngularJ64(14//2)
        h15 = AngularJ64(15//2);  h16 = AngularJ64(16//2);  h17 = AngularJ64(17//2);  h18 = AngularJ64(18//2);  h19 = AngularJ64(19//2)
        h20 = AngularJ64(20//2);  h21 = AngularJ64(21//2);  h22 = AngularJ64(22//2);  h23 = AngularJ64(23//2);  h24 = AngularJ64(24//2)
        h25 = AngularJ64(25//2);  h26 = AngularJ64(26//2);  h27 = AngularJ64(27//2)
        #
        # Distinguish due to different j-values
        if       j == AngularJ64(1//2)        &&  N in [0, 2]
            wa = [ SubshellTerm(h1, h1, N, 0, h0, 0) ]
        elseif   j == AngularJ64(1//2)        &&  N in [1]
            wa = [ SubshellTerm(h1, h0, N, 1, h1, 0) ]
            #
            #
        elseif   j == AngularJ64(3//2)        &&  N in [0, 4]
            wa = [ SubshellTerm(h3, h2, N, 0, h0, 0) ]
        elseif   j == AngularJ64(3//2)        &&  N in [1, 3]
            wa = [ SubshellTerm(h3, h1, N, 1, h3, 0) ]
        elseif   j == AngularJ64(3//2)        &&  N in [2]
            wa = [ SubshellTerm(h3, h2, N, 0, h0, 0), SubshellTerm(h3, h0, N, 2, h4, 0) ]
            #
            #
        elseif   j == AngularJ64(5//2)        &&  N in [0, 6]
            wa = [ SubshellTerm(h5, h3, N, 0, h0, 0) ]
        elseif   j == AngularJ64(5//2)        &&  N in [1, 5]
            wa = [ SubshellTerm(h5, h2, N, 1, h5, 0) ]
        elseif   j == AngularJ64(5//2)        &&  N in [2, 4]
            wa = [ SubshellTerm(h5, h3, N, 0, h0, 0), SubshellTerm(h5, h1, N, 2, h4, 0), SubshellTerm(h5, h1, N, 2, h8, 0) ]
        elseif   j == AngularJ64(5//2)        &&  N in [3]
            wa = [ SubshellTerm(h5, h2, N, 1, h5, 0), SubshellTerm(h5, h0, N, 3, h3, 0), SubshellTerm(h5, h0, N, 3, h9, 0) ]
            #
            #
        elseif   j == AngularJ64(7//2)        &&  N in [0, 8]
            wa = [ SubshellTerm(h7, h4, N, 0, h0, 0) ]
        elseif   j == AngularJ64(7//2)        &&  N in [1, 7]
            wa = [ SubshellTerm(h7, h3, N, 1, h7, 0) ]
        elseif   j == AngularJ64(7//2)        &&  N in [2, 6]
            wa = [ SubshellTerm(h7, h3, N, 1, h7, 0), SubshellTerm(h7, h2, N, 2, h4, 0), SubshellTerm(h7, h2, N, 2, h8, 0), 
                   SubshellTerm(h7, h2, N, 2,h12, 0) ]
        elseif   j == AngularJ64(7//2)        &&  N in [3, 5]
            wa = [ SubshellTerm(h7, h3, N, 1, h7, 0), SubshellTerm(h7, h1, N, 3, h3, 0), SubshellTerm(h7, h1, N, 3, h5, 0), 
                   SubshellTerm(h7, h3, N, 1, h9, 0), SubshellTerm(h7, h1, N, 3,h11, 0), SubshellTerm(h7, h1, N, 3,h15, 0) ]
        elseif   j == AngularJ64(7//2)        &&  N in [4]
            wa = [ SubshellTerm(h7, h4, N, 0, h0, 0), SubshellTerm(h7, h2, N, 2, h4, 0), SubshellTerm(h7, h0, N, 4, h8, 0),
                   SubshellTerm(h7, h2, N, 2,h12, 0), SubshellTerm(h7, h0, N, 4, h4, 0), SubshellTerm(h7, h0, N, 4, h8, 0),
                   SubshellTerm(h7, h0, N, 4,h10, 0), SubshellTerm(h7, h0, N, 4,h16, 0) ]
            #
            #
        elseif   j == AngularJ64(9//2)        &&  N in [0, 10]
            wa = [ SubshellTerm(h9, h5, N, 0, h0, 0) ]
        elseif   j == AngularJ64(9//2)        &&  N in [1, 9]
            wa = [ SubshellTerm(h9, h4, N, 1, h9, 0) ]
        elseif   j == AngularJ64(9//2)        &&  N in [2, 8]
            wa = [ SubshellTerm(h9, h5, N, 0, h0, 0), SubshellTerm(h9, h3, N, 2, h4, 0), SubshellTerm(h9, h3, N, 2, h8, 0), 
                   SubshellTerm(h9, h3, N, 2,h12, 0), SubshellTerm(h9, h3, N, 2,h16, 0) ]
        elseif   j == AngularJ64(9//2)        &&  N in [3, 7]
            wa = [ SubshellTerm(h9, h4, 1, N, h9, 0), SubshellTerm(h9, h2, N, 3, h3, 0), SubshellTerm(h9, h2, N, 3, h5, 0), 
                   SubshellTerm(h9, h2, N, 3, h7, 0), SubshellTerm(h9, h2, N, 3, h9, 0), SubshellTerm(h9, h2, 3, N,h11, 0), 
                   SubshellTerm(h9, h2, N, 3,h13, 0), SubshellTerm(h9, h2, N, 3,h15, 0), SubshellTerm(h9, h2, 3, N,h17, 0), 
                   SubshellTerm(h9, h2, 3, N,h21, 0) ]
        elseif   j == AngularJ64(9//2)        &&  N in [4, 6]
            wa = [ SubshellTerm(h9, h5, N, 0, h0, 0), SubshellTerm(h9, h3, N, 2, h4, 0), SubshellTerm(h9, h3, N, 2, h8, 0), 
                   SubshellTerm(h9, h3, N, 2,h12, 0), SubshellTerm(h9, h3, N, 2,h16, 0), SubshellTerm(h9, h1, 4, N, h0, 0), 
                   SubshellTerm(h9, h1, N, 4, h4, 0), SubshellTerm(h9, h1, 4, N, h6, 0), SubshellTerm(h9, h1, N, 4, h8, 1), 
                   SubshellTerm(h9, h1, 4, N, h8, 2), SubshellTerm(h9, h1, N, 4,h10, 0), SubshellTerm(h9, h1, 4, N,h12, 1), 
                   SubshellTerm(h9, h1, N, 4,h12, 2), SubshellTerm(h9, h1, N, 4,h14, 0), SubshellTerm(h9, h1, N, 4,h16, 0), 
                   SubshellTerm(h9, h1, N, 4,h18, 0), SubshellTerm(h9, h1, 4, N,h20, 0), SubshellTerm(h9, h1, N, 4,h24, 0) ]
        elseif   j == AngularJ64(9//2)        &&  N in [5]
            wa = [ SubshellTerm(h9, h4, 1, N, h9, 0), SubshellTerm(h9, h2, N, 3, h3, 0), SubshellTerm(h9, h2, N, 3, h5, 0), 
                   SubshellTerm(h9, h2, N, 3, h7, 0), SubshellTerm(h9, h2, N, 3, h9, 0), SubshellTerm(h9, h2, 3, N,h11, 0), 
                   SubshellTerm(h9, h2, N, 3,h13, 0), SubshellTerm(h9, h2, N, 3,h15, 0), SubshellTerm(h9, h2, 3, N,h17, 0), 
                   SubshellTerm(h9, h2, 3, N,h21, 0), SubshellTerm(h9, h0, 5, N, h1, 0), SubshellTerm(h9, h0, 5, N, h5, 0), 
                   SubshellTerm(h9, h0, N, 5, h7, 0), SubshellTerm(h9, h0, N, 5, h9, 0), SubshellTerm(h9, h0, N, 5,h11, 0), 
                   SubshellTerm(h9, h0, 5, N,h13, 0), SubshellTerm(h9, h0, N, 5,h15, 0), SubshellTerm(h9, h0, N, 5,h17, 0), 
                   SubshellTerm(h9, h0, N, 5,h19, 0), SubshellTerm(h9, h0, N, 5,h25, 0) ]
        else     error("stop a")
        end
        
        return(wa)
    end
    
    
    
    #######################################################################################################################
    #######################################################################################################################

    reducedCfp1half = Matrix{Tuple{Int64,Int64,Int64}}(undef,2,2);    for i=1:2  for j=1:2  reducedCfp1half[i,j] = (0,0,0)   end  end  
    reducedCfp3half = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);    for i=1:3  for j=1:3  reducedCfp3half[i,j] = (0,0,0)   end  end
    reducedCfp5half = Matrix{Tuple{Int64,Int64,Int64}}(undef,6,6);    for i=1:6  for j=1:6  reducedCfp5half[i,j] = (0,0,0)   end  end
    reducedCfp7half = Matrix{Tuple{Int64,Int64,Int64}}(undef,14,14);  for i=1:14 for j=1:14 reducedCfp7half[i,j] = (0,0,0)   end  end
    
    
    reducedCfp1half[1,1] = (-1,      4,     1)
    reducedCfp1half[1,2] = (-1,      4,     1)
    
    
    reducedCfp3half[1,1] = (-1,      4,     1);     reducedCfp3half[1,2] = (-1,      4,     1)
 
 
 
 
    W_00_one_half = Matrix{Tuple{Int64,Int64,Int64}}(undef,2,2);       for i=1:2  for j=1:2  W_00_one_half[i,j] = (0,0,0)   end  end
    
    
