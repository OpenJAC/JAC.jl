
# Compile and print the examples (files) for developing the JAC tools  
println("Examples, tests & development of JAC, ordered by a few main branches: \n")
println("A)  Examples, tests & development of the electronic structure part.")
println("B)  Examples, tests & development of atomic amplitudes")
println("C)  Examples, tests & development of atomic properties")
println("D)  Examples, tests & development of basic atomic processes")
println("E)  Examples, tests & development of composed atomic processes")
println("F)  Examples, tests & development of atomic cascades")
println("G)  Examples, tests & development of symbolic evaluations of Racah expressions")
println("H)  Examples, tests & development of semiempirical computations")
println("I)  Examples, tests & development of atomic response computations")
println("J)  Examples, tests & development of plasma computations \n")
print(  "Enter a Letter (A, B, ...) to select the associated branch of examples; ... enter Q to quit:  ")
char = read(stdin, Char) 
println(" ")

if      char == 'A'
    println("A)  Examples, tests & development of the electronic structure part.")
    println("-------------------------------------------------------------------")
    println("Aa) Apply & test several radial e-e potentials.")
    println("Ab) Apply & test the  a SCF field: B-spline primitives and one-particle spectra in a local potential.")
    println("Ac) Apply & test the CI part for an internally generated neon multiplet without Breit interaction.")
    println("Ad) Apply & test for the frequency-independent Breit interaction for an internally generated neon multiplet.")
    println("Ae) Apply & test for the QED model corrections to the level structure of atoms and ions.")
    println("Af) Generate, normalize & test for continuum orbitals in a local potential.")
    println("Ag) Apply & test the jj-LS transformation of levels from a given multiplet.")
    println("Ah) Apply & test a mean-field basis, mean-field multiplets as well as a configuration-interaction (CI) expansions.")
    println("Ai) Apply & test for restricted-active-space (RAS) expansions.")
    println("Aj) Apply & test for a Green (-function) expansion.")
    println("Ak) Apply & test the computation of spin-angular coefficients.")
    println("Al) Apply & test for parallel computating methods/techniques with Julia.")
    #
elseif  char == 'B'
    println("B)  Examples, tests & development of atomic amplitudes")
    println("------------------------------------------------------")
    println("Ba) Apply & test the dipole, em multipole and momentum-transfer amplitudes.")
    println("Bb) Apply & test the parity non-conservation, Schiff moment and anapole moment amplitudes.")
    #
elseif  char == 'C'
    println("C)  Examples, tests & development of atomic properties")
    println("------------------------------------------------------")
    println("Ca) Apply & test the Einstein module with ASF from an internally generated multiplet.")
    println("Cb) Apply & test the Hfs module for HFS A,B parameters and hyperfine representation with ASF from an internally generated multiplet.")
    println("Cc) Apply & test the IsotopeShift module with ASF from an internally generated multiplet.")
    println("Cd) Apply & test the LandeZeeman module with ASF from an internally generated multiplet.")
    println("Ce) Apply & test the FormFactor module with ASF from an internally generated multiplet.")
    println("Cf) Apply & test the reduced 1- and 2-particle density matrices & natural orbitals.")
    println("Cg) Apply & test the PlasmaShift module with ASF from an internally generated multiplet.")
    println("Ch) Apply & test the RadiativeOpacity module with ASF from an internally generated multiplet.")
    println("Ci) Apply & test the AlphaVariation module with ASF from an internally generated multiplet.")
    println("Cj) Apply & test the MultipolePolarizibility module with ASF from an internally generated multiplet.")
    println("Ck) Apply & test the DecayYield module with ASF from an internally generated multiplet.")
    #
elseif  char == 'D'
    println("D)  Examples, tests & development of basic atomic processes")
    println("-----------------------------------------------------------")
    println("Da) Apply & test the PhotoEmission module with ASF from an internally generated initial- and final-state multiplet.")
    println("Db) Apply & test the PhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dc) Apply & test the PhotoIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dd) Apply & test the PhotoRecombination module with ASF from an internally generated initial- and final-state multiplet.")
    println("De) Apply & test the AutoIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Df) Apply & test the Dielectronic module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    println("Dg) Apply & test the RayleighCompton module with ASF from an internally generated initial- and final-state multiplets.")
    println("Dh) Apply & test the MultiPhotonDeExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    #
    println("Di) Apply & test the DoubleAutoIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dj) Apply & test the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dk) Apply & test the RadiativeAuger module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dl) Apply & test the ImpactExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dm) Apply & test the InternalRecombination module with ASF from an internally generated initial and final-state multiplet.")
    println("Dn) Apply & test the TwoElectronOnePhoton module with ASF from an internally generated initial and final-state multiplet.")
    println("Dp) Apply & test the ParticleScattering module with ASF from an internally generated initial and final-state multiplet.")
    println("Dq) Apply & test the InternalConversion module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dr) Apply & test the CoulombIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Ds) Apply & test the CoulombExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    #
elseif  char == 'E'
    println("E)  Examples, tests & development of composed atomic processes")
    println("--------------------------------------------------------------")
    println("Ea) Apply & test the BeamPhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    println("Eb) Apply & test the  PhotoExcitationFluores module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    println("Ec) Apply & test the PhotoExcitationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    println("Ed) Apply & test the Photoionization and PlasmaShift module with ASF from an internally generated initial- and final-state multiplet.")
    println("Ef) Apply & test the Auger and PlasmaShift modules with ASF from an internally generated initial- and final-state multiplet.")
    #
    println("Eg) Test of the ImpactExcitationAutoion module with ASF from an internally generated initial-, intermediate and final-state multiplet.")
    println("Ei) Test of the MultiPhotonIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Ej) Test of the MultiPhotonDoubleIon module with ASF from an internally generated initial- and final-state multiplet.")
    println("Ek) Test of the PairAnnihilation1Photon.jl module with ASF from an internally generated initial- and final-state multiplet.")
    println("El) Test of the PhotoIonizationFluores module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    println("Em) Test of the PhotoIonizationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    #
elseif  char == 'F'
    println("F)  Examples, tests & development of atomic cascades")
    println("----------------------------------------------------")
    println("Fa) Compute & test a three-step cascade model following the 1s-3p photo-excitation of Si^+.")
    println("Fb) Simulate & test a four-step cascade model following the 1s-3p photo-excitation of Si^+.")
    println("Fc) Three-step cascade computations and simulations for the decay of the neon 1s^-1 3p hole states: AverageSCA model.")
    println("Fd) Photoabsorption of Ne^+: AverageSCA model.")
    println("Fe) Compute & test a dielectronic recombination cascade and DR plasma rate coefficients.")
    println("Ff) Compute & test a radiative recombination cascade and RR plasma rate coefficients.")
    println("Fg) Compute & test an expansion opacity (cascade, not yet).")
    println("Fh) Compute & test a photoabsorption cascade computations.")
    println("Fi) Compute & test an electron-impact excitation (cascade) for lithium-like Ne.")
    #
elseif  char == 'G'
    println("G)  Examples, tests & development of symbolic evaluations of Racah expressions")
    println("------------------------------------------------------------------------------")
    println("Ga) Symbolic evaluation by means of special values & recursion relations of the Wigner n-j symbols.")
    println("Gb) Symbolic simplification of Racah expressions by means of sum rules.")
    println("Gc) Symbolic simplification of recoupling coefficients by means of sum rules.")
    println("Gd) Symbolic simplification of Racah expressions including spherical harmonics and Wigner rotation matrices.")
    println("Ge) Symbolic simplification of spherical tensor operators, matrix elements and spherical amplitudes.")
    #
elseif  char == 'H'
    println("H)  Examples, tests & development of semiempirical computations")
    println("---------------------------------------------------------------")
    println("Ha) Apply & test the semi-empirical computation of ADK rates.")
    println("Hb) Compute & test (empirical) impact-ionization cross sections.")
    #
elseif  char == 'I'
    println("I)  Examples, tests & development of atomic response computations")
    println("-----------------------------------------------------------------")
    println("Ia) Apply & test the HighHarmonic module to compute high-harmonic spectra in single-electron approximation.")
    println("Ib) Apply & test the StrongField module to calculate energy and momentum distributions of photoelectron in SFA.")
    println("Ic) Apply & test the StrongField module to demonstrate the Coulomb asymmetry in ATI azimuthal angular distributions.")
    println("Id) Apply & test the StrongField2 module to compute in ATI photoelectron momentum distributions.")
    #
elseif  char == 'J'
    println("J)  Examples, tests & development of plasma computations")
    println("--------------------------------------------------------")
    println("Ja) Tests of average-atom computations.")
    #
elseif  char == 'Q'
    return( nothing )
else    
    println("Unknown Letter; redo ...")
end

