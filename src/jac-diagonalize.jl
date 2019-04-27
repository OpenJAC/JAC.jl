
"""
`Basics.diagonalize()`  ... diagonalizes matrices of various kinds and by different methods 

  + `("matrix: Julia, eigfact", matrix::Array{Float64,2})`  ... to apply the standard eigfact() method from Julia for a quadratic, full matrix;
                                                                nothing about the symmetry of the matrix is assumed here; an eigen::JAC.Eigen
                                                                is returned.
"""
function Basics.diagonalize(sa::String, matrix::Array{Float64,2})
    if       sa == "matrix: Julia, eigfact" 
        # Use the standard eigfact() method from Julia for a quadratic, full matrix   
        wa = eigen( matrix )
        ##x println("diagonalize-ab: values  = ", wa[:values] )
        ##x println("diagonalize-ac: vectors = ", wa[:vectors] )
        vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
        for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d])    end
        wc = Basics.Eigen( wa.values, vectors )
        return( wc )
    elseif   sa == "matrix: Julia method B"     error("Not yet implemented")
    else     error("Unsupported keystring = $sa")
    end
end


"""
  + `("generalized eigenvalues: Julia, eigfact", matrixA::Array{Float64,2}, matrixB::Array{Float64,2})`  
                                ... to apply the standard eigfact() method from Julia for a generalized eigenvalue problem with two quadratic   
                                    (full) matrices; nothing about the symmetry of the matrix is assumed here; an eigen::JAC.Eigen is returned.
"""
function Basics.diagonalize(sa::String, matrixA::Array{Float64,2}, matrixB::Array{Float64,2})
    if       sa == "generalized eigenvalues: Julia, eigfact" 
        # Use the standard eigfact() method from Julia for two quadratic, full matrices   
        wa = eigen( matrixA, matrixB )
        ##x println("diagonalize-ab: values  = ", wa[:values] )
        ##x println("diagonalize-ac: vectors = ", wa[:vectors] )
        ##x vectors = Vector{Float64}[];    wb = wa[:vectors];    d = size(wb)[1]
        vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
        # for  i = 0:d-1    push!(vectors, real(wb[i*d+1:i*d+d]) )    end
        for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d] )    end
        ##x wc = JAC.Eigen( wa[:values], vectors )
        # wc = JAC.Eigen( real(wa.values), vectors )
        wc = Basics.Eigen( wa.values, vectors )
        return( wc )
    else     error("Unsupported keystring = $sa")
    end
end

