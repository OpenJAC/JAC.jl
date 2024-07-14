

#==

Task:  Evaluate symbolically the recoupling coefficient for the (recoupling of) three
-----  angular momenta by means of a simple Racah expression.
       This input has been adapted from ../examples/example-Ga.jl
==#

using SymEngine

j1   = Basic(:j1);    j2 = Basic(:j2);    j3 = Basic(:j3);    J12 = Basic(:J12);    J23 = Basic(:J23);    J = Basic(:J)
m1   = Basic(:m1);    m2 = Basic(:m2);    m3 = Basic(:m3);    M12 = Basic(:M12);    M23 = Basic(:M23);    M = Basic(:M)

w3ja = W3j(J12, j3, J, M12, m3, -M)                                 # Specify a Wigner 3-j symbol algebraically.
w3jb = W3j(j1, j2, J12, m1, m2, -M12)       
w3jc = W3j(j2, j3, J23, m2, m3, -M23)     
w3jd = W3j(j1, J23, J, m1, M23, -M) 

# Specify the overall Racah expression in terms if its summations, phases, weight and the product of Wigner 3-j symbols.
rex  = RacahExpression(RacahExpression(), 
                       summations = [m1, m2, m3, M12, M23],
                       phase      = -J12 + 2*j3 - 2*M - 2*j1 - M12 - M23 + J23,
                       weight     = (2*J+1) * sqrt( (2*J12+1)*(2*J23+1) ), 
                       w3js       = [w3ja, w3jb, w3jc, w3jd] ) 

rex  = RacahAlgebra.evaluate(rex);   @show rex
rex  = RacahAlgebra.evaluate(rex)


#==
Main output:  
------------

>> Apply sum rule for three W3j -- Sum(m4,m5,m6) ...
>> Apply sum rule for two W3j -- Sum(np,nq) (-1)^(-np-nq) ...
>>> two W3j -- Sum(np,nq) (-1)^(-np-nq)::::  phase = -J + 2*J23 + 2*M + j1 - j2 - j3   without = Basic[M]   
    zeroTerms = Basic[2*J + 2*J23 + 2*j1, -M + M23 + m1, 2*j1 + 2*m1, 2*J23 + 2*M23, 2*J - 2*M, M - M23 - m1, 2*j1 - 2*m1, 2*J23 - 2*M23, 2*J + 2*M, 2*J12 + 2*j1 + 2*j2, 2*J23 + 2*j2 + 2*j3, 2*J + 2*J12 + 2*j3]
>>> two W3j -- Sum(np,nq) (-1)^(-np-nq):  newPhase = J + 2*J23 - 3*j1 - j2 - j3 - 4*m1  done !!!!
>>> two W3j -- Sum(np,nq) (-1)^(-np-nq):  newPhase = -J - j1 - j2 - j3  now shorter !!!!
rex = (-1)^(-J - j1 - j2 - j3)  (sqrt((1 + 2*J12)*(1 + 2*J23)))  W6j{j1, J23, J; j3, J12, j2}  

(-1)^(-J - j1 - j2 - j3)  (sqrt((1 + 2*J12)*(1 + 2*J23)))  W6j{j1, J23, J; j3, J12, j2}  

==#
