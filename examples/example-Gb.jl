#
println("Gb) Tests of (empirical) impact-ionization cross sections.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "S-plus.sum")
grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  false
    #
    # Calculate partial EII cross section for K-subshell of Ne IX in different approximations 
    # BEBapproximation(), BEDapproximation(), RelativisticBEBapproximation(), RelativisticBEDapproximation(), MUIBEDapproximation()
    approx      = ImpactIonization.BEBapproximation()
    iEnergies   = [300., 400., 500., 600., 700., 800., 900., 1000., 1200., 1400., 1600., 1800., 2000., 3000., 4000., 5000., 6000., 8000., 
                   10000., 20000., 40000., 60000., 80000., 100000.] ## unit: eV. The incident energies should be > epsilon_subshell.
    eEnergies   = [10., 50., 100.]                                  ## unit: eV.
    shells      = Basics.generateShellList(2,2, [0,1])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("1s^2 2s^2 2p^5")]
    name        = "EII cross section for L-subshell of Si VI."
    nucModel    = Nuclear.Model(14.0)
    eiiSettings = ImpactIonization.Settings(approx, iEnergies, eEnergies, false, true, true, selection)
    ##X basis       = Basics.performSCF(configs, nucModel, grid, AsfSettings()) 
    ##x display     = ImpactIonization.displayCrossSections(basis, grid, nucModel, shells, eiiSettings, configs)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    @show comp
    perform(comp)
    #
elseif false
    #
    # Calculate partial EII cross section for for L-subshells of Bi in different approximations
    approx      = ImpactIonization.RelativisticBEBapproximation() 
    iEnergies   = [14000., 16000., 18000., 20000., 30000., 40000., 50000., 60000., 70000., 80000., 
                   100000., 200000., 400000., 600000., 800000., 1000000., 2000000., 4000000., 6000000., 8000000., 
                   10000000., 100000000, 1000000000] ## unit: eV. The incident energies should be > epsilon_subshell.
    eEnergies   = [10.,  50., 100.]                   
    shells      = Basics.generateShellList(2,2, [0,1])  
    selection   = ShellSelection(true, shells=shells)
    configs     = [Configuration("[Xe] 4f^14 5d^10 6s^2 6p^3")] 
    name        = "EII cross section for L-subshells of Bi."
    nucModel    = Nuclear.Model(83.0)
    ##x basis       = Basics.performSCF(configs, nucModel, grid, AsfSettings()) 
    eiiSettings = ImpactIonization.Settings( approx, iEnergies, eEnergies, false, true, true, selection )
    ##x display     = ImpactIonization.displayCrossSections(basis, grid, nucModel, iEnergies, shells, eiiSettings, configs)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    @show comp
    perform(comp)
    #
elseif true
    #
    # Calculate partial EII cross section for M-subshells of U in different approximations
    approx      = ImpactIonization.FittedBEDapproximation()
    iEnergies   = [6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000., 30000., 40000., 60000., 70000., 
                   80000., 90000., 100000., 120000., 140000., 160000., 180000.,  200000., 500000., 1000000., 
                   10000000., 100000000., 1000000000.] ## unit: eV. The incident energies should be > epsilon_subshell.
    eEnergies   = [10., 50., 100.]  
    shells      = Basics.generateShellList(3,3, [0,1,2])
    selection   = ShellSelection(true, shells=shells)
    configs     = [Configuration("[Rn] 5f^3 6d 7s^2")]
    name        = "EII cross section for M-subshells of U."
    nucModel    = Nuclear.Model(92.0)
    eiiSettings = ImpactIonization.Settings(approx, iEnergies, eEnergies, false, true, true, selection)
    ##x basis       = Basics.performSCF(configs, nucModel, grid, AsfSettings()) 
    ##x display     = ImpactIonization.displayCrossSections(basis, grid, nucModel, iEnergies, shells, eiiSettings, configs)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    @show comp
    perform(comp)
    #
end
 
