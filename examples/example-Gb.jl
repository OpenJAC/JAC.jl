#
println("Gb) Tests of (empirical) impact-ionization cross sections.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "S-plus.sum")
grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  true
    # Calculate partial impact ionization cross sections in different approximations
    approx      = ImpactIonization.BEBapproximation()  ## BEBapproximation(), RelativisticBEBapproximation()
    iEnergies   = [100., 1000.]                        ## in user-specified units
    fEnergies   = [10.,  50., 100.]                    ## in user-specified units
    shells      = Basics.generateShellList(1,2, [0,1])
    selection   = ShellSelection(true, shells=shells)
    configs     = [Configuration("[Ar] 4s^2")]
    name        = "First EII cross section computations."
    #
    eiiSettings = ImpactIonization.Settings( approx, iEnergies, fEnergies, false, true, true, selection )
    comp        = Empirical.Computation(name, Nuclear.Model(20.0), grid, configs, eiiSettings)
    #
    @show comp
    perform(comp)
    #
elseif false
    #
end
 
