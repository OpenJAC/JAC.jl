## Demo files
!!! tips "Tips"
    - The demo files are designed with Pluto notebook. A detaild guideline of Pluto Installation is available at (https://plutojl.org/#install).
    - Users can also find the files in the **demo** folder JAC package.

- [SCF and CI](https://github.com/OpenJAC/JAC.jl/blob/master/demos/B1-Pluto-compute-SCF%2BCI-carbon-III.jl)
- [Transition Probabilities and Lifetime](https://github.com/OpenJAC/JAC.jl/blob/master/demos/demo-A-FeX-lifetimes.jl)
- [Dielectronic Recombination](https://github.com/OpenJAC/JAC.jl/blob/master/demos/demo-D-H-like-DR.jl)
- [Auger rates](https://github.com/OpenJAC/JAC.jl/blob/master/demos/demo-B-Ne-KLL-Auger.jl)

## Compute the low-lying levels of C``^{2+}`` ``1s^{2}(2s^{2}+2s2p+2p^{2})`` : 

Perform SCF and configuration interaction calculations for these low-lying levels.

Let us first invoke JAC, either from the own source-code basis or simply by **using JAC**



The low-lying levels (level structure) of beryllium-like ions, and especially of C``^{2+}``, has been calculated in many case studies
in the literature. While the level structure of these ions is still quite simple, it exhibits a considerable admixture of the 
``2s^22p^2`` configuration already for the ``1s^{2}2s^{2}`` ``^{1}S_{0}`` ground level.

We here show how the low-lying levels of C``^{2+}`` can be readily calculated in JAC by either following the default settings or
by specifying further details for both, the SCF and configuration-interaction (CI) computations. As usual, we first need to 
specify a radial grid as well as the nuclear model for the subsequent computations:

```@repl
using JAC

grid     = Radial.Grid(true)
nucModel = Nuclear.Model(6., "Fermi")

multiplet = SelfConsistent.performSCF([Configuration("1s^2 2s^2")], nucModel, grid, AsfSettings())

multiplet2 = SelfConsistent.performSCF([Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2p^2")], nucModel, grid, AsfSettings())
```

For a quick computation of the ground level of C``^{2+}`` ions, we can simply use the **standard settings** as given by **AsfSettings()**:

