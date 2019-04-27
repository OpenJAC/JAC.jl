#
println("Ea) Three-step cascade after 1s-3p photo-excitation of Si^+: Configuration model.")
#
setDefaults("print summary: open", "zzz-Cascade.sum")
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

wa = Cascade.Computation("Si+ (1s^-1) cascade", Nuclear.Model(18.), 
                         JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600), AsfSettings(),
                         JAC.Cascade.AverageSCA(), [JAC.Radiative],          #  JAC.Cascade.AverageSCA()
                         [Configuration("1s^1 2s^2 2p^6 3s^2 3p^2")], [(1, 1.)],  0, 0, Shell[], Shell[], JAC.Cascade.Step[])
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

## using JLD
## save("zzz-Cascade.jld", wb)
