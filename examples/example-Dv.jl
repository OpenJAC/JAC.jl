
println("Dv) Test of the PairAnnihilation1Photon module with ASF from an internally generated initial- and final-state multiplet.")

@warn("\n\n !!! This example does not work properly at present !!! \n\n")
#
wa = Atomic.Computation("xx",  Nuclear.Model(26.), JAC.Radial.Grid("grid: exponential"), false, [JAC.NoProperty],
     Configuration[],
     Configuration[Configuration("1s^2")], 
     Configuration[], 
     Configuration[Configuration("1s") ],
     Einstein.Settings(),  Hfs.Settings(),  IsotopeShift.Settings(),   PlasmaShift.Settings(), LandeZeeman.Settings(),
     JAC.PairA1P,  PairAnnihilation1Photon.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [0.], true, false, Tuple{Int64,Int64}[])  )

wb = perform(wa)


#     [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")], 

