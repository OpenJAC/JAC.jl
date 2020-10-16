
println("Di) Test of the RayleighCompton module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
#
setDefaults("print summary: open", "zzz-RayleighCompton.sum")

if  true
    # Rayleigh scattering on the ground-configuration levels of B-like neon
    asfSettings      = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    rayleighSettings = RayleighCompton.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [10., 100.], 2., true, true, true, 
                                                ExpStokes(), [SolidAngle(1.0, 0.0)], LineSelection())
    grid             = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")  
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs=[Configuration("1s^2 2s^2 2p")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p")], 
                            process = Compton(), processSettings = rayleighSettings )
    @show wa
    wb = perform(wa)
end

setDefaults("print summary: close", "")
