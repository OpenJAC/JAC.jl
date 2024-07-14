
#==

Task:  Calculate the K-LL Auger rates (spectrum) of atomic neon.
-----  This input has been adapted from ../examples/example-De.jl
==#

setDefaults("unit: energy", "eV")   
setDefaults("unit: rate", "1/s")   
setDefaults("method: continuum, Galerkin")  
setDefaults("method: normalization, Alok")
grid          = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
augerSettings = AutoIonization.Settings(AutoIonization.Settings(), calcAnisotropy = true, printBefore = true, 
                                        lineSelection = LineSelection(true, indexPairs=[(1,0)])  )
    
comp = Atomic.Computation(Atomic.Computation(), name="K-LL Auger rates of neon.", grid=grid, nuclearModel=Nuclear.Model(10.), 
                          initialConfigs  = [Configuration("1s 2s^2 2p^6")],
                          finalConfigs    = [Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), 
                                             Configuration("1s^2 2p^6")], 
                          processSettings = augerSettings )

perform(comp)


#==
Main output:  
------------
 
>> Auger computations for line No = 10   ...
>> (Re-) Define a storage array for various B-spline matrices:
>> Continuum B-spline-Galerkin orbital for energy=2.7415e+01,  kappa=-1 [mpt=1102, r[mtp]=8.0967e+00, smallest eigenvalue=-2.2071e-13].
>> Radial potential with effective charge Zbar=2.0000e+00 (Delta-Zbar=0.0000e+00) at r=1.0061e+01 a.u.
>> Normalization with Coulomb functions:   r = 8.096663160516204,   iPhase = -1.2594105154673,   cPhase = 0.1484543599400503 
Compute (CoulombInteraction()) Auger matrix of dimension 3 x 1 in the continuum- and initial-state bases for the transition [1- ...] and for partial wave s_1/2 ... done. 
 
  Auger rates and intrinsic angular parameters: 

  ----------------------------------------------------------------------------------------------------------
       i-level-f           i--J^P--f           Energy      Electron energy    Auger rate        alpha_2      
                                                [eV]             [eV]            [1/s]                        
  ----------------------------------------------------------------------------------------------------------
         1 --    1       1/2 + --> 2 +      -2.631778e+03    8.106703e+02    5.672459e+10       0.0000e+00    
         1 --    2       1/2 + --> 1 +      -2.631778e+03    8.105888e+02    1.254009e+08      -0.0000e+00    
         1 --    3       1/2 + --> 0 +      -2.631778e+03    8.105505e+02    6.070283e+09       0.0000e+00    
         1 --    4       1/2 + --> 2 +      -2.631778e+03    8.072747e+02    2.470948e+14       0.0000e+00    
         1 --    5       1/2 + --> 0 +      -2.631778e+03    8.041767e+02    2.952512e+13       0.0000e+00    
         1 --    6       1/2 + --> 2 -      -2.631778e+03    7.841101e+02    1.509872e+13       0.0000e+00    
         1 --    7       1/2 + --> 1 -      -2.631778e+03    7.840327e+02    9.210682e+12      -0.0000e+00    
         1 --    8       1/2 + --> 0 -      -2.631778e+03    7.839936e+02    3.095883e+12       0.0000e+00    
         1 --    9       1/2 + --> 1 -      -2.631778e+03    7.720310e+02    7.379028e+13      -0.0000e+00    
         1 --   10       1/2 + --> 0 +      -2.631778e+03    7.460056e+02    2.419943e+13       0.0000e+00    
  ----------------------------------------------------------------------------------------------------------
 
  Auger lifetimes, total rates and widths:
 
  --------------------------------------------------------------------------------------------------------
     Level       J^P        Lifetime        Total rate                          Widths                    
                              [sec]            [1/s]            Hartrees         Kaysers           eV   
  --------------------------------------------------------------------------------------------------------
         1      1/2 +     2.487081e-15     4.020778e+14      9.725797e-03    2.134566e+03    2.646524e-01    
  --------------------------------------------------------------------------------------------------------

==#
