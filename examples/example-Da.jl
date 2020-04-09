#
println("Da) Test of the PhotoEmission module with ASF from an internally generated initial- and final-state multiplet.")
#
# setDefaults("unit: rate", "a.u.")
setDefaults("unit: rate", "1/s")
setDefaults("print summary: open", "zzz-radiative.sum")

if  true

    grid = Radial.Grid(true)
    setDefaults("standard grid", grid)
    defaultsSettings = PhotoEmission.Settings()
    photoSettings = PhotoEmission.Settings(defaultsSettings, multipoles=[E1, M1, E2, M2], gauges=[UseCoulomb, UseBabushkin], printBefore=true)
    
    comp = Atomic.Computation(Atomic.Computation(), name="Energies and Einstein coefficients for the spectrum Fe X",  
              grid=grid, nuclearModel=Nuclear.Model(26.);
              initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],
              finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")], 
              process = Radiative, processSettings = photoSettings ); 
    @show comp
    perform(comp)          
              
elseif false
    #== wa = Atomic.Computation("xx",  Nuclear.Model(36.);
                            initialConfigs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs  =[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")], 
                            process = JAC.Radiative, 
                            processSettings=PhotoEmission.Settings([E1, M1, E2, M2], [JAC.UseCoulomb, JAC.UseBabushkin], true, true, 
                                true, Tuple{Int64,Int64}[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0), (15,0), (16,0)], 0., 0., 10000. ) )  ==#
elseif true                           
    wa = Atomic.Computation(Atomic.Computation(), name="Energies and Einstein coefficients for the spectrum Fe X", 
                            grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                            initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],
                            finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")], 
                            process = Radiative, 
                            processSettings = PhotoEmission.Settings([E1, M1, E2, M2], [JAC.UseCoulomb, JAC.UseBabushkin], 
                                              false, true, false, Tuple{Int64,Int64}[], 0., 0., 10000. ) );
    wb = @time( perform(wa) )
end
setDefaults("print summary: close", "")


## , Configuration("[Ne] 3s^2 3p^4 3d")  Configuration("[Ne] 3s 3p^6")
