#
println("Ea) Three-step cascade after 1s-3p photo-excitation of Si^+: Configuration model.")
#
JAC.define("print summary: open", "zzz-Cascade.sum")
JAC.define("method: continuum, asymptotic Coulomb")    ## JAC.define("method: continuum, Galerkin")
JAC.define("method: normalization, pure sine")         ## JAC.define("method: normalization, pure Coulomb")    JAC.define("method: normalization, pure sine")

wa = Cascade.Computation("Si+ (1s^-1) cascade", Nuclear.Model(18.), 
                         JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600), AsfSettings(),
                         JAC.Cascade.AverageSCA(), [JAC.Radiative],          #  JAC.Cascade.AverageSCA()
                         [Configuration("1s^1 2s^2 2p^6 3s^2 3p^2")], [(1, 1.)], 1, 0, Shell[], Shell[], JAC.Cascade.Step[])
wb = perform(wa; output=false)
JAC.define("print summary: close", "")

