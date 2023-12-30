#
println("Gb) Tests of (empirical) impact-ionization cross sections.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "S-plus.sum")
grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  true
    #
    # Calculate partial EII cross section for K-subshell of Ne IX in different models 
    # BEBmodel(), BEDmodel(), RelativisticBEBmodel(), RelativisticBEDmodel(), MUIBEDmodel()
    approx      = ImpactIonization.RelativisticBEBmodel()
    iEnergies   = [1200., 1400., 1600., 1800., 2000., 3000., 4000., 5000., 6000., 8000., 
                   10000., 20000., 40000., 60000., 80000., 100000.] ## unit: eV. The incident energies should be > epsilon_subshell.
    shells      = Basics.generateShellList(1,1, [0])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("1s^2")]
    name        = "EII cross section for K-subshell of Ne IX."
    nucModel    = Nuclear.Model(10.0)
    eiiSettings = ImpactIonization.Settings(approx, iEnergies, true, true, selection)
    ##X basis       = Basics.performSCF(configs, nucModel, grid, AsfSettings()) 
    ##x display     = ImpactIonization.displayCrossSections(basis, grid, nucModel, shells, eiiSettings, configs)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    @show comp
    perform(comp)
    #
elseif false
    #
    # Calculate partial EII cross section for M-subshell of Kr XIX in different models 
    # BEBmodel(), BEDmodel(), RelativisticBEBmodel(), RelativisticBEDmodel(), MUIBEDmodel()
    approx      = ImpactIonization.RelativisticBEBmodel()
    iEnergies   = [800, 1000, 1200., 1400., 1600., 1800., 2000., 3000., 4000., 5000., 6000., 8000., 
                   10000., 20000., 40000., 60000., 80000., 100000.] ## unit: eV. The incident energies should be > epsilon_subshell.
    shells      = Basics.generateShellList(3,3, [0,1,2])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("[Ar]")]
    name        = "EII cross section for K-subshell of Ne IX."
    nucModel    = Nuclear.Model(36.0)
    eiiSettings = ImpactIonization.Settings(approx, iEnergies, false, true, true, selection)
    ##X basis       = Basics.performSCF(configs, nucModel, grid, AsfSettings()) 
    ##x display     = ImpactIonization.displayCrossSections(basis, grid, nucModel, shells, eiiSettings, configs)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    @show comp
    perform(comp)
    #
elseif false
    #
    # Calculate partial EII cross section for M-subshells of U in different models
    approx      = ImpactIonization.MUIBEDmodel()
    iEnergies   = [6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000., 30000., 40000., 60000., 70000., 
                   80000., 90000., 100000., 120000., 140000., 160000., 180000.,  200000., 500000., 1000000., 2000000., 
                   3000000., 5000000., 10000000.] ## unit: eV. The incident energies should be > epsilon_subshell.
    shells      = Basics.generateShellList(3,3, [0,1,2])
    configs     = [Configuration("[Rn] 5f^3 6d 7s^2")]
    name        = "EII cross section for M-subshells of U."
    nucModel    = Nuclear.Model(92.0)
    eiiSettings = ImpactIonization.Settings(approx, iEnergies, false, true, false, selection)
    ##x basis       = Basics.performSCF(configs, nucModel, grid, AsfSettings()) 
    ##x display     = ImpactIonization.displayCrossSections(basis, grid, nucModel, iEnergies, shells, eiiSettings, configs)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    @show comp
    perform(comp)
    #
elseif false
    #
    # Calculate partial EII cross section for for L-subshells of Bi in different models
    approx      = ImpactIonization.BEBmodel() 
    iEnergies   = [10000.0, 12000., 14000., 16000., 18000., 20000., 30000., 40000., 50000., 60000., 70000., 80000., 
                   100000., 200000., 400000., 600000., 800000., 1000000., 2000000., 4000000., 6000000., 8000000., 
                   10000000., 100000000, 1000000000] ## unit: eV. The incident energies should be > epsilon_subshell.                
    shells      = Basics.generateShellList(2,2, [0,1])  
    selection   = ShellSelection(true, shells=shells)
    configs     = [Configuration("[Xe] 4f^14 5d^10 6s^2 6p^3")] 
    name        = "EII cross section for L-subshells of Bi."
    nucModel    = Nuclear.Model(83.0)
    ##x basis       = Basics.performSCF(configs, nucModel, grid, AsfSettings()) 
    eiiSettings = ImpactIonization.Settings( approx, iEnergies, false, true, false, selection )
    ##x display     = ImpactIonization.displayCrossSections(basis, grid, nucModel, iEnergies, shells, eiiSettings, configs)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    @show comp
    perform(comp)
    #
end
 
