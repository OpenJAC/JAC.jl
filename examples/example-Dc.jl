#
println("Dc) Test of the PhotoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoIonization.sum")
setDefaults("method: continuum, asymptotic Coulomb")  ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure Coulomb")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

# Test April 2020
if  true
defaultSettings = PhotoIonization.Settings()
setDefaults("unit: energy", "eV")
e1 = convertUnits("energy: to atomic",  600.)
e2 = convertUnits("energy: to atomic",  1000.)

photoSettings = PhotoIonization.Settings(defaultSettings, gauges=[UseCoulomb, UseBabushkin], photonEnergies=[e1, e2], 
                                         calcAnisotropy=true, calcPartialCs=true, printBefore=true)
                                         
grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.3e-2, NoPoints = 600)
wa   = Atomic.Computation(Atomic.Computation(), name="Photoionization of Ar: cross sections + beta's", 
              grid=grid, nuclearModel=Nuclear.Model(18.),
              initialConfigs=[Configuration("[Ne] 3s^2 3p^6")],
              finalConfigs  =[Configuration("[Ne] 3s^2 3p^5")], 
              process = JAC.Photo,  processSettings=photoSettings)
@show wa
perform(wa)                                         

error("stop a")
end

wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(false), nuclearModel=Nuclear.Model(20.), 
                        initialConfigs=[Configuration("[Ne] 3s^2 3p^6")],
                        finalConfigs  =[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6") ], 
                        process = JAC.Photo, 
                        processSettings=PhotoIonization.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [50., 60., 70., 80., 90., 100.], false, true, true, true, 
                                        false, Tuple{Int64,Int64}[], JAC.ExpStokes(1., 0., 0.)) )

@show wa
## wb = perform(wa)
setDefaults("print summary: close", "")

