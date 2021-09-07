#
println("Hb) Tests of the StrongField module to calculate energy and momentum distributions of photoelectron in SFA.")
#

dataName        = "Neon_StrongFieldIonization" #If files are exported below (plot, data), all filenames will include dataName


#########################################################
#-------------------PHYSICAL SETTINGS--------------------
#########################################################


#_____________1. Atomic target______________________________________________

#An initial and final level for the strong-field ionization process will be computed
#In the strong-field computation, the highest lying one-electron  orbital will be taken as "active electron"

#Specify the initial atomic and final ionic states
initialName     = "Li ground-state"
nuclearModel    = Nuclear.Model(3.39336) #Adjust the nuclear charge so that the ionization potential of the selected active electron fits (Li: 3.39336, Ne: 10.463599, Ar: 18.50435055, Kr: 36.529677, Xe: 54.5504443)
refConfigInitial= [Configuration("[He] 2s")]

finalName       = "He ground-state"
refConfigFinal  = [Configuration("[He]")]

hydrogenic      = false             #True: Use hydrogenic wave functions for the initial state (quantum numbers n,l,m and the ionization potential are taken from initialLevel)
hydrogenic1s    = false             #True (AND hydrogenic = true): Same as hydrogenic, but quantum numbers n,l,m of initial level are ignored and a 1s state is used instead (ionization potential is taken from initialLevel)

#___________________________________________________________________________


#_____________2. Laser field and electron continuum_________________________

#---Define laser beam parameters---------------
wavelength      = 800.                                                          #nm
intensity       = 5e14                                                          #W/cm^2
CEP             = 0.                                                            #carrier-envelope phase 0 <= CEP < 2*pi
envelope        = Pulse.SinSquaredEnvelope(8)                                   #InfiniteEnvelope(), SinSquaredEnvelope(np), GaussianEnvelope(np) with the number of optical cycles np
polarization    = Basics.RightElliptical(0.36)                                  #RightCircular(), LeftCircular(), 
                                                                                #RightElliptical(0 <= ellipticity <= 1), LeftElliptical(0 <= ellipticity <= 1)
                                                                                #NOTE: Elliptical polarization is only implemented for SinSquaredEnvelope so far
                                                                                
#---Choose electron continuum in the laser field-------------------
volkov          = StrongField.FreeVolkov()                                #FreeVolkov(), CoulombVolkov(Z) with the charge of the ion Z, or DistortedVolkov()

#-----AUTOMATIC: derived quantities and settings-------------------
omega           = convertUnits("energy: from wavelength [nm] to atomic", wavelength)
intensity       = convertUnits("intensity: from W/cm^2 to atomic", intensity)
A0              = Pulse.computeFieldAmplitude(intensity, omega)
beam            = Pulse.PlaneWaveBeam(A0, omega, CEP)                            #Arguments: A0, omega, CEP #Only PlaneWaveBeam possible

#___________________________________________________________________________


#_____________3. Observable that will be computed__________________________

observable      = StrongField.SfaEnergyDistribution(pi/2, pi/2, 10, 10*omega)          #Arguments: theta, phi, number of energies, maximum energy
#observable      = StrongField.SfaMomentumDistribution(pi/2, 30, 100, 10*omega)         #Arguments: theta, number of azimuthal angles, number of energies, maximum energy
#observable      = StrongField.SfaAzimuthalAngularDistribution(pi/2, 100, 2.21*omega)   #Arguments: theta, number of azimuthal angles, energy
#observable      = StrongField.SfaPolarAngularDistribution(pi/2, 100, 2.5*omega)        #Arguments: phi, number of azimuthal angles, energy

#___________________________________________________________________________


#########################################################
#-----------------COMPUTATIONAL SETTINGS-----------------
#########################################################

asfSettings     = AsfSettings(AsfSettings(), generateScf=true)
rGrid           = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0) #Radial grid for atomic computation and matrix elements

multipoles      = [E1]              #Only "dipole" ( [E1] ) SFA implemented so far
gauges          = [UseCoulomb]      #Only Coulomb (velocity) gauge implemented so far
mAverage        = true              #True: Average over initialState projections mj or ml and sum over final state spin projection msp
coupledBasis    = true              #True: compute in |j l mj> angular momentum basis; False: compute in |l ml> angular momentum basis

printAmplitudes = true              #True: Print all requested amplitudes during computation

exportData      = true              #True: Probabilities will be saved to dataName-observable-j.csv as (E_p,P) or (p_x,p_y,P)
exportPlot      = true              #True: Energy distributions or Angular distributions will be plotted and saved to dataName-observable.pdf
plotEnergyScale = "omega"           #Energy scaling used in plot ("atomic" or "omega")
plotProbabilityScale = "linear"     #Energy scaling used in plot ("linear" or "log" - "log" here means log10)
exportRadialWF  = true              #True: Radial wave functions of initial states used in StrongField computations will be saved to dataName-initial_radial_wave_function-j.csv
plotRadialWF    = true              #True: Radial wave functions of initial states used in StrongField computations will be plotted and saved to dataName-initial_radial_wave_function.pdf


##########################################################################################################################################################


##############################################################
#-----------------COMPUTATION AND DATA EXPORT-----------------
##############################################################

#_____________1. Atomic computations: Compute the initial and final levels for the SFA computation_____
AtomicComp              = Atomic.Computation(Atomic.Computation(), name=initialName, grid=rGrid, nuclearModel=nuclearModel, configs=refConfigInitial,  asfSettings=asfSettings )
AtomicData              = perform(AtomicComp, output=true)
initialLevel            = AtomicData["multiplet:"].levels[1]
#
AtomicComp              = Atomic.Computation(Atomic.Computation(), name=finalName, grid=rGrid, nuclearModel=nuclearModel, configs=refConfigFinal,  asfSettings=asfSettings )
AtomicData              = perform(AtomicComp, output=true)
finalLevel              = AtomicData["multiplet:"].levels[1]


#_____________2. StrongField computation(s)___________
sfaSettings         = StrongField.Settings(multipoles, gauges, printAmplitudes, coupledBasis, hydrogenic, hydrogenic1s, mAverage)
SFIComp1            = StrongField.Computation(observable, nuclearModel, rGrid, initialLevel, finalLevel, beam, envelope, polarization, volkov, sfaSettings)
SFIData1            = StrongField.perform(SFIComp1, output=true)

compArray           = [SFIComp1] #Multiple computations can be combined here into an array and will be plotted and exported below
dataArray           = [SFIData1]


#_____________3. Export data and generate plots___________

if exportData
    StrongField.exportData(compArray, dataArray, dataName)
end

if exportPlot
    StrongField.plot(compArray, dataArray, plotEnergyScale, plotProbabilityScale, dataName)
end

if exportRadialWF
    StrongField.exportRadialWavefunctions( compArray, dataName, plotRadialWF )
end
