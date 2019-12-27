
println("Dv) Test of the PairAnnihilation1Photon module with ASF from an internally generated initial- and final-state multiplet.")

@warn("\n\n !!! This example does not work properly at present !!! \n\n")
#
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(36.) )

wb = perform(wa)

