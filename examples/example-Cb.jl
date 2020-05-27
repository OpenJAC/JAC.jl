#
println("Cb) Test of the Hfs module for HFS A,B parameters and hyperfine representation with ASF from an internally generated multiplet.")
#
setDefaults("unit: energy", "Hz")
setDefaults("print summary: open", "zzz-Hfs.sum")
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), 
                        nuclearModel=Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0), properties=[HFS()],
                        # grid=JAC.Radial.Grid(true),
                        configs=[Configuration("[Ne] 3s"), Configuration("[Ne] 3p"), Configuration("[Ne] 3d")],
                        hfsSettings=Hfs.Settings(true, true, false, false, false, false, false, Int64[] ) )

wb = perform(wa)
setDefaults("print summary: close", "")


#==
wa = Atomic.Computation("xx",  Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0); properties=[HFS()],
                        # grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 0., NoPoints = 600),
                        configs=[Configuration("[Ne] 3s"), Configuration("[Ne] 3p"), Configuration("[Ne] 3d")],
                        hfsSettings=Hfs.Settings(true, true, false, false, false, false, false, Int64[] ) )

wa = Atomic.Computation("xx",  Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0); properties=[HFS()],
                        # grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 0., NoPoints = 600),
                        configs=[Configuration("[Ne] 3s^2"), Configuration("[Ne] 3s 3p"), Configuration("[Ne] 3p^2")],
                        hfsSettings=Hfs.Settings(true, true, false, false, false, false, false, Int64[] ) )

wa = Atomic.Computation("xx",  Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0); properties=[HFS()],
                        # grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 0., NoPoints = 600),
                        configs=[Configuration("[Ne] 3s^2 3d"), Configuration("[Ne] 3d^3")],
                        hfsSettings=Hfs.Settings(true, true, false, false, false, false, false, Int64[] ) )
==#
