
# Jena Atomic Calculator (JAC) for the computation of atomic representations, processes and cascades


**JAC**, the **Jena Atomic Calculator** JAC is a (relativistic) electronic structure code for the 
computation of (atomic many-electron) interaction amplitudes, properties as well as a large number of 
excitation and decay processes, such as photoexcitation and emission, photoionization, autoionization
as well as a few selected electron-impact processes. JAC is open-source Julia package suitable for 
-- almost all -- open-shell atoms and ions across the whole periodic table.


A primary goal in developing JAC was to establish  a **general and easy-to-use toolbox for the atomic physics 
community**, including an interface that is equally accessible for spectroscopic scientists, theoreticians and 
code developers. Apart from its simple use, however, we wish to provide a modern code design and a reasonable 
detailed documentation of the code. In particular, we aim that many (typical) computations can be described 
in a descriptive  language similar to how the given task is summarized in a spoken or written language. 
Many further details about this toolbox can be found on github: https://github.com/OpenJAC/JAC.jl  
from where the code, a first Readme file, News as well as a detailed *User Guide, Compendium & Theoretical
Background to JAC* can be freely downloaded by the user. Beside of the central features of JAC, 
the Readme also outlines a few of the current limitations with regard to the shell structure or the size 
and speed of computations. The code in this repository is distributed under the MIT licence.


Once you have installed Julia on your computer, the JAC package can be installed like any other package in your 
environment by either using the package manager (pressing <Alt> ]) and by typing

pkg> add https://github.com/OpenJAC/JAC.jl                          or by typing:
> import Pkg; Pkg.add("https://github.com/OpenJAC/JAC.jl")

at the Julia prompt, and later called by 
> using JAC

as well as updated like any other package in Julia. --- Indeed, that's all. The JAC package itself makes use 
of standard Julia packages, such as SpecialFunctions, GSL and QuadGK, but these dependencies are maintained 
automatically by Julia's package manager.


For a quickstart of the JAC toolbox, since simple *demos* (demonstration files) have been prepared:

+ demo-A-FeX-lifetimes.jl                  ... to compute transition probabilities & lifetimes for Fe^9+ ions.
+ demo-B-Ne-KLL-Auger.jl                   ... to compute the K-LL Auger rates of atomic neon.
+ demo-C-Lshell-photoionization-neon.jl    ... to compute the 2s and 2p photoionization cross sections for neon
+ demo-D-Li-like-DR.jl
+ demo-E-Ar-form-factors.jl 
+ demo-F-recoupling-symbolic-evaluation.jl

which can be called either "line-wise" or by simply typing at the Julia prompt:
> using JAC 
> include("../demos/demo-A-FeX-lifetimes.jl")

These files also display the (final) part of the outcome of the computations -- as a comment
as well as for quick review & comparison.  --- Good luck and enjoy using JAC.

