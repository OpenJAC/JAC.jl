
#
println("Ed) Photoabsorption of Ne^+: AverageSCA model.")
using JLD
#
if  false
    # Calculation of the Ne^+  2s^2 2p^5; 2s --> 2p, 3p, 4p, 5p excitation cross sections
    setDefaults("print summary: open", "Ne-plus-photoexcitation.sum")
    setDefaults("method: continuum, asymptotic Coulomb")                                                      ## "method: continuum, Galerkin"
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, NoPoints=1110)
    name = "Photoabsorption calculations for Ne^+ for energies = [10, 20] a.u."
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                               scheme=Cascade.PhotonExcitationScheme([PhotoExc()], [E1], 0.1, 1000.0, 1, [Shell("2s"), Shell("2p")], 
                                                                     [Shell("2s"), Shell("2p"), Shell("3p"), Shell("4p"), Shell("5p")]),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif false
    # Calculation of the Ne^+  2s^2 2p^5; photo-ionization cross sections
    setDefaults("print summary: open", "Ne-plus-photoionization.sum")

    name = "Photoionization of Si- "
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, NoPoints=1110)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                               scheme=Cascade.PhotonIonizationScheme([Photo()], 1, [20.0, 30.0]),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    
else
    # Simulation of the Ne^+  2s^2 2p^5; photo-ionization cross sections
    setDefaults("print summary: open", "Ne-plus-photoabsorption-simulation.sum")

    data = [JLD.load("zzz-cascade-ionizing-computations-2020-06-23T09.jld"), JLD.load("zzz-cascade-excitation-computations-2020-06-23T09.jld")]
    ##x JLD2.@load "zzz-cascade-ionizing-computations-2020-06-23T09.jld"
    ##x resIon  = results
    ##x JLD2.@load "zzz-cascade-excitation-computations-2020-06-23T09.jld"
    ##x resExc = results
    ##x data = [resExc]

    name = "Simulation of photoabsorption for Ne^+"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.PhotoAbsorption(), 
                              settings=Cascade.SimulationSettings(0., 0., 0., 0., 30., [(1, 1.0)]), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")

end
