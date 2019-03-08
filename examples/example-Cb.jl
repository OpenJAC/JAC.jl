#
println("Cb) Test of the Hfs module for HFS A,B parameters and hyperfine representation with ASF from an internally generated multiplet.")
#
JAC.define("print summary: open", "zzz-Hfs.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0); properties=[JAC.HFS],
                        configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                        hfsSettings=Hfs.Settings(true, true, true, true, true, true, false, Int64[] ) )

wb = perform(wa)
JAC.define("print summary: close", "")



##                        hfsSettings=Hfs.Settings(true, true, false, true, true, true, false, Int64[] ) )
