
#
println("Cc) Test of the IsotopeShift module with ASF from an internally generated multiplet.")
#
setDefaults("print summary: open", "zzz-IsotopeShift.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0); properties=[JAC.Isotope],
                        # grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 0., NoPoints = 600),
                        configs=[Configuration("[Ne] 3s^2"), Configuration("[Ne] 3s 3p"), Configuration("[Ne] 3p^2")],
                        isotopeSettings=IsotopeShift.Settings(true, true, true, true, false, Int64[], "method-1") )

wb = perform(wa)
setDefaults("print summary: close", "")


#==
wa = Atomic.Computation("xx",  Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0); properties=[JAC.Isotope],
                        # grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 0., NoPoints = 600),
                        configs=[Configuration("[Ne] 3s"), Configuration("[Ne] 3p"), Configuration("[Ne] 3d")],
                        isotopeSettings=IsotopeShift.Settings(true, true, true, true, false, Int64[], "method-1") )

==#
