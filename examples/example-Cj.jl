
#
println("Cj) Test of the GreenFunction module for a given list of bound-state configurations.")
#
setDefaults("print summary: open", "zzz-GreenFunction.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[Green], 
                        configs=[Configuration("[Ne] 3s^2 3p^6")],
                        greenSettings=GreenFunction.Settings(GreenFunction.SingleCSFwithoutCI(), Basics.DeExciteSingleElectron(), 
                                                             5, [0, 1, 2], [LevelSymmetry(1//2,Basics.plus), LevelSymmetry(1//2,Basics.plus)], 
                                                             true, false, Int64[]) )

wb = perform(wa)
setDefaults("print summary: close", "")


