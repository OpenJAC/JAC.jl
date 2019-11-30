

The following features have been recently implemented and up-loaded into JAC's master code.
===========================================================================================



## 2020



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
