

The following features have been recently implemented and up-loaded into JAC's master code.
===========================================================================================



## 2021

* **Strong-field ionization (together with Birger Böning):** A first module (and examples) have been worked out to 
    apply the strong-field approximation (SFA) in a partial-wave reprensentation of the initial and Volkov states; 
    this makes the SFA suitable for dealing with realistic electronic structures of atoms and with different gauges
    of the radiation field. *(September'21)*
* **Extended and improved User Guide & Compendium:**  *(August'21)*
* **Write-up of atomic Green functions:** Cf. S. Fritzsche  & A. Surzhykov,
    Molecules 2021, 26, 2660; https://doi.org/10.3390/molecules26092660   *(May'21)*
* **Extended dielectronic recombination (DR) calculations:** Apart from the DR resonance strength, the
    DR (plasma) rate coefficients can be calculated as well.*(March'21)*
* **Write-up of atomic cascade computations:** Cf. S. Fritzsche, P. Palmeri & S. Schippers,
    Symmetry 2021, 13(3), 520; https://doi.org/10.3390/sym13030520   *(March'21)*
* **Julia implementation of the spin-angular integration (together with Gediminas Gaigalas)** The former 
    dependence of JAC upon the (Fortran) ANCO module for the spin-angular coefficients has now been removed
    and replaced by the (Julia) SpinAngular module. This simplifies significantly the installation of the 
    JAC toolbox and will support further developments on cascades, reduced density matrices, etc.
    *(March'21)*
* **Isotope-shift calculations:** Set-up of a first module for isotope-shift calculations *(February'21)*



## 2020

* **Extended and improved User Guide & Compendium:**  *(October'20)*
* **Enlarged number of second-order processes:** Several second-order processes are now partly supported,
    including photoexcitation & fluorescence,  photoexcitation & autoionization, Rayleigh scattering 
    of light, multi-photon excitation & decay, double-Auger decay, photo-double ionization, 
    radiative-Auger decay. *(October'20)*
* **Enlarged number of atomic processes:** The presently supported atomic processes now include 
    photoemission, photoexcitation, photoionization, photorecombination, Auger & autoionization,
    dielectronic recombination. *(September'20)*
* **Atomic cascades:** The presently supported cascade computations include decay cascades as well as
    their excitation by photoionization, photoexcitation or electron-capture. These excitation also
    strongly affect the complexity of a cascade. *(September'20)*
* **Autoionization & dielectronic recombination (DR):** See `? AutoIonization.Settings` or `Dielectronic.Settings`
    for the calculation of Auger rates and DR resonance strength. *(August'20)*
* **Photoionization cross sections:** See `? PhotoIonization.Settings` for the calculation of photoionization
    cross sections and angular parameters. *(July'20)*
* **Generation & use of continuum orbitals:** See `? setDefaults` how one can choose different methods
    for the generation ("method: continuum, Galerkin", "method: continuum, asymptotic Coulomb", ...)
    and normalization ("method: normalization, pure Coulomb", ...) of continuum orbitals. *(June'20)*
* **New level & line selectors:** See `? LevelSelection`, `? LineSelection` and `PathwaySelection`; this
    enables one to easily select individual levels, lines or pathways in the computation of properties
    and processes. *(June'20)*
* **New and extended tutorials:** A good number of new tutorials have been implemented and extended,
    inluding the symbolic evaluation of expressions from Racah's algebra. *(May'20)*
* **Improved continuum orbitals:** Use of Galerkin method for generating continuum orbitals with 
    well-defined energy. *(April'20)*
* **Extended and improved User Guide & Compendium:**  *(January-March'20)*



## 2019

* **New and extended tutorials:** A good number of new tutorials have been implemented and extended,
    inluding SCF computations and different representations of atomic states. *(30.11.19)*
* **Green (function) expansions:** Implememtation of simple approximate Green expansions, though still without
    any continuum interactions. *(26.11.19)*
* **Restricted active-space (RAS) expansions:** Automatic generation of RAS expansions for selected levels
    or symmetries. *(22.11.19)*
* **Atomic representations:** Different representations of wave and Green functions can now be distinguished;
    cf. User Guide, section 4.1. *(15.10.19)*
* **Racah algebra:** Symbolic evaluation of expressions from Racah's algebra, based on special value and
    sum rules; cf. User Guide, chapter 15. *(10.8.19)*
* **jjJ - LSJ transformation:** LSJ expansion & notations of atomic levels with one (nonrelativistic) 
    open shell. *(12.7.19)*
* **Atomic form factors:** Standard and modified form factors for atoms in spherical-symmetric levels
    *(15.6.19)*
* **Debye-Hückel plasma potential:** Implementation of the plasma shifts of energy levels into the
    CI computation of a multiplet; cf. User Guide, section 6.1.e.   *(26.3.19)*
* **Breit interaction:**  Implementation of the frequency-independent Breit interaction amplitudes 
    that can be accessed via `AsfSettings(..., breitCI=true, ...)` in all configuration-interaction
    computations.  *(4.1.19)*
