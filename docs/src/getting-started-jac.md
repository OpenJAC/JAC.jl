# Getting started with JAC (in REPL)

!!! info
    JAC user guide pdf .... link

```
using JAC
```

### Welcome to **JAC**, the **Jena Atomic Calculator**

... that provides various tools for performing atomic (structure) calculations of different kinds and complexities. Apart from the computation of atomic (many-electron) amplitudes, properties and processes, **JAC supports interactive, restricted-active space (RAS) and cascade computations**. It also help perform a few simple *hydrogenic* and *semi-empirical* estimates as well as simplify symbolic expressions from Racah's algebra. --- Let's first use  `? JAC` in order to obtain more information about this toolbox:

```
? JAC
```

```@docs ; canonical=false
JAC
```
H'm, this tells us a lot of details which we still need to better understand. To quickly list the atomic properties, that have been (partly) considered in JAC, we can use `? Details.properties`   or some other of the listed calls:

```
? Details.properties
```

```@docs ; canonical=false
Details.properties
```

In the design of JAC, we first of all **aim for a precise language** that (i) is simple enough for both, seldom and a more frequent use of this package, (ii) highlights the underlying physics and (iii) avoids most technical slang that is often unnecessary but quite common to many other codes. An intuitive picture about the level or hyperfine structure of an atom, its properties as well as possible excitation and/or decay processes should (always) come first in order to generate the desired data: By making use of suitable data types (`struct`), **we indeed wish to introduce a language close to the underlying formalism.** --- While JAC is overall based on a rather large number $(> 300)$ of such types, a few simple examples are:

  + (atomic) `Shell`:                 $\quad$1s, 2s, 2p, ...
  + `Subshell`:                       $\quad$1s_1/2, 2s_1/2, 2p_1/2, 2p_3/2, ...
  + (electron) `Configuration`:       $\quad$1s^2 2s^2 2p^6 3s $\quad$  or $\quad$  [Ne] 3s, ...
  + `Level`:                          $\quad$1s^2 2s^2  ^1S_0, ...
  
and many other terms (types) that we shall explain later.  


Let us simply start, for instance, with specifying and assigning the $1s$ and $2p$ shells:

```@repl
using JAC   # hide
w1s = Shell("1s")
w2p = Shell("2p")
```

Similarly, we can readily specify and assign any (relativistic) subshell:

```@repl
using JAC   # hide
Subshell("2p_1/2"),   Subshell("2p_3/2")  
```

In JAC, we make use of these `Shell`'s and `Subshell`'s whenever they will naturally occur in describing the level structure or the excitation, decay or occupation of an atom, and this both at input and output. If you have *forgotten* how to specify such a subshell (constructor), simply *ask*:

```
? Subshell
```

```@docs ; canonical=false
Subshell
```

Of course, we can interactively also specify any electron configuration:

```
? Configuration
```

```
using JAC # hide
wc1 = Configuration("1s^2 2s^2 2p^5")
wc2 = Configuration("[Ar] 4s^2 3d^5")
```
!!! info 
    For specific processes users can find the list of types and functions in the API Reference

This input just shows three (very) simple examples and how the details of some computation can be readily specified in line with our basic understanding of the atomic shell model. One can use  `? Details.datatypes`  in order to see a more complete list of most data structures that are speficic to the JAC module ... and which will give you a very **first impression about the size of the JAC program**.

```
? Details.datatypes
```

```@docs ; canonical=false
Details.datatypes
```

This list gives further details why Julia (and JAC) is a very suitable and powerful framework for running -- many-electron -- atomic computations. 

Of course, there are many other features that make Julia & JAC as powerful as it is: For example, the user may pre-define and overwrite the **units** in which he wishes to communicate with JAC. These units determine how (most of) the input data are interpreted as well as output data are displayed in tabulations or to screen. The current defaults settings for the units can be seen by typing:

```@repl
using JAC   # hide
Basics.display("settings")
```

which show that energies are taken/printed in eV, rates in 1/s, etc. Apart from modifying these defaults directly in the source code, the can be *overwritten* by the user at any time of the program executation. This is done by means of the function

```
? Defaults.setDefaults
```

```@docs
Defaults.setDefaults
```
which enables one to re-define various **global values** of JAC. If we wish to enter/display energies in **Kaysers** or cross sections in atomic units, we can simply type:

```@repl
using JAC   # hide
Defaults.setDefaults("unit: energy", "Kayser")
Defaults.setDefaults("unit: cross section", "a.u.")
```

Here, again `nothing` is returned but the corresponding global constants are now changed.

```@repl    
using JAC   # hide
Basics.display("settings")
```
Apart from the default units, one can similarly *overwrite* the method that is use for the generation and normalization of continuum orbitals and several others. Although called *global*, the corresponding values can be accesses just in two ways. (i) The **global constants**, such as the electron mass, the speed of light, the fine-structure constant $\alpha$, etc., are accessed via the function:

```
? Defaults.getDefaults
```

```@docs
Defaults.getDefaults
```

```@repl
using JAC   # hide
Defaults.getDefaults("alpha")
Defaults.getDefaults("electron rest energy")
Defaults.getDefaults("unit: energy")
```

(ii) These **global values** are frequently applied in order to -- internally or externally -- convert physical numbers into units of the same dimension. This is done by the function:

```
? Defaults.convertUnits
```

```@docs
Defaults.convertUnits
```

This function is called at many places within JAC to generate tables where all physical data are printed out in the pre-specified units:

```@repl
using JAC   # hide
Defaults.convertUnits("energy: from atomic", 1.0)
```

With the given user-selection, this is equivalent to:

```@repl
using JAC # hide
Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
```

In JAC, the call of this function is often combined with some proper formatting of the results, such as:

```@repl
import Pkg  ; Pkg.add("Printf") ; # hide
using JAC   # hide
using Printf
@sprintf("%.4e", Defaults.convertUnits("energy: from atomic", 1.0))
```

