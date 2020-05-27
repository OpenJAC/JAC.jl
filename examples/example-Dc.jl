#
println("Dc) Test of the PhotoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
setDefaults("method: normalization, pure Coulomb")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

if  false
    # For tutorial
    defaultSettings = PhotoIonization.Settings()
    setDefaults("unit: energy", "eV")
    e1 = convertUnits("energy: to atomic",  600.)
    e2 = convertUnits("energy: to atomic",  1000.)

    photoSettings = PhotoIonization.Settings(defaultSettings, gauges=[UseCoulomb, UseBabushkin], photonEnergies=[e1], 
                                             calcAnisotropy=true, calcPartialCs=true, printBefore=true)
                                         
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.3e-2, NoPoints = 1000)
    wa   = Atomic.Computation(Atomic.Computation(), name="Photoionization of Ar: cross sections + beta's", 
                              grid=grid, nuclearModel=Nuclear.Model(18.),
                              initialConfigs=[Configuration("[Ne] 3s^2 3p^6")],
                              finalConfigs  =[Configuration("[Ne] 3s^2 3p^5")], 
                              process = Photo(),  processSettings=photoSettings)
    @show wa
    perform(wa)                                         

elseif false
    # Test case from FAC/User Guide
    defaultSettings = PhotoIonization.Settings()
    setDefaults("unit: energy", "eV")
    setDefaults("unit: cross section", "Mbarn")
    e1 = 1.01 * 75.20911382458269  * 27.21
    e2 = 1.20 * 75.20911382458269  * 27.21
    e3 = 2.00 * 75.20911382458269  * 27.21
    e4 = 3.00 * 75.20911382458269  * 27.21
    e5 = 4.00 * 75.20911382458269  * 27.21
    e6 = 6.00 * 75.20911382458269  * 27.21

    photoSettings = PhotoIonization.Settings(defaultSettings, gauges=[UseCoulomb, UseBabushkin], photonEnergies=[e1, e2, e3, e4, e5, e6], ### 
                                             calcAnisotropy=false, calcPartialCs=false, printBefore=true)
                                         
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 0.3e-2, NoPoints = 1200)
    wa   = Atomic.Computation(Atomic.Computation(), name="Photoionization of Li-like iron", 
                              grid=grid, nuclearModel=Nuclear.Model(26.),
                              initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                              finalConfigs  =[Configuration("1s^2")], 
                              process = Photo(),  processSettings=photoSettings)
    @show wa
    perform(wa)                                         

elseif true
    # Test case 
    defaultSettings = PhotoIonization.Settings()
    setDefaults("unit: energy", "eV")
    setDefaults("unit: cross section", "Mbarn")
    e1 =  40.0
    e2 =  60.0
    e3 =  80.0
    e4 = 100.0
    e5 = 120.0

    photoSettings = PhotoIonization.Settings(defaultSettings, gauges=[UseCoulomb, UseBabushkin], photonEnergies=[e1, e2, e3, e4, e5], ### , e2, e3, e4, e5, e6
                                             calcAnisotropy=false, calcPartialCs=false, printBefore=true)
                                         
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.0e-2, NoPoints = 900)
    wa   = Atomic.Computation(Atomic.Computation(), name="Photoionization of neutral neon", 
                              grid=grid, nuclearModel=Nuclear.Model(10.),
                              initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                              finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6")], 
                              process = Photo(),  processSettings=photoSettings)
    @show wa
    perform(wa)                                         
end

setDefaults("print summary: close", "")

