### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ f18a9488-13b8-11f0-342b-f99afaf3ac19
begin
	using Pkg
	Pkg.develop(path="/home/fritzsch/fri/JAC.jl")
end

# ╔═╡ fd4ef417-311f-43aa-a601-bf8a9d66b890
using JAC

# ╔═╡ 79374db5-1fbd-4499-b8d9-c7b840bdee52
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(150px, 10%);
    	padding-right: max(150px, 10%);
	}
</style>
"""

# ╔═╡ 8d491602-2e68-4ce9-b949-859f91cfee9d
md"""
# Compute the low-lying levels of C``^{2+}`` ``1s^{2}(2s^{2}+2s2p+2p^{2})`` : 
Perform SCF and configuration interaction calculations for these low-lying levels.
"""

# ╔═╡ 75061d53-4fb7-45a6-b8bb-c199d52479da
md"""
Let us first invoke JAC, either from the own source-code basis or simply by **using JAC**
"""

# ╔═╡ 956d7bd6-d59f-4daa-88b6-9b8879918d79
md"""
The low-lying levels (level structure) of beryllium-like ions, and especially of C``^{2+}``, has been calculated in many case studies
in the literature. While the level structure of these ions is still quite simple, it exhibits a considerable admixture of the 
``2s^22p^2`` configuration already for the ``1s^{2}2s^{2}`` ``^{1}S_{0}`` ground level.

We here show how the low-lying levels of C``^{2+}`` can be readily calculated in JAC by either following the default settings or
by specifying further details for both, the SCF and configuration-interaction (CI) computations. As usual, we first need to 
specify a radial grid as well as the nuclear model for the subsequent computations:
"""

# ╔═╡ 5b8c341c-4923-4146-bdb1-5b6a7cc20c17
begin
	grid     = Radial.Grid(true)
	nucModel = Nuclear.Model(6., "Fermi")
end

# ╔═╡ 1945551c-c5e1-4a75-897d-ef07205bb341

	md"""
	For a quick computation of the ground level of C``^{2+}`` ions, we can simply use the **standard settings** as given by **AsfSettings()**:
	"""

# ╔═╡ 561451e9-62b6-4a02-9d33-1393c575e890
multiplet = SelfConsistent.performSCF([Configuration("1s^2 2s^2")], nucModel, grid, AsfSettings())

# ╔═╡ 9787ed44-f14f-45a7-b2df-a40c9c01679a
multiplet2 = SelfConsistent.performSCF([Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2p^2")], nucModel, grid, AsfSettings())

# ╔═╡ 081bb1f4-fe80-4428-a234-e27970ae33be
md"""
From the comparison of the two ground-state energies, we see that the admixture of the ``2p^{2}`` configuration has lowered the (total) ground
state energy by about 1.8 eV, a rather remarkable admixture, as the ``^{3}P_{0}`` is just 6.6 eV above of the ground level.

Further control about these electronic computations can be obtained by modifying the (so-called) **settings**. 
In general, all computations of the electronic structure, properties and processes as well as all more advanced computations
can be controlled quite in details by various settings that are associated to the different computational requests. 
The SCF and configuration interaction calculations are controlled by **AsfSettings** that specify all details for the generation of the ASF.
We can first have a look at the internal representation of these settings by typing **AsfSettings** in the LiveDocs (or ? AsfSettings in the Julia REPL). \
We can also invoke the default values of these settings by calling the constructor ManyElectron.AsfSettings() or just AsfSettings():
"""

# ╔═╡ 0e57775c-b714-4795-9630-2e82da20871b
defaultAsfSettings = AsfSettings()

# ╔═╡ cd37b365-a6e7-4328-a513-7f9c29d4644e
md"""
From this list, we easily see that the self-consistent field is by default based on a (mean) Dirac-Fock-Slater potential with parameter ``x_\alpha = 1.0``, a choice that we could overwrite by a mean-Core-Hartree or any of the pre-defined potentials. At present, however, no full treatment of the exchange interaction has yet been implemented in this first release of the program. The standard settings also show that the SCF is usually based on on a maximum number `maxIterationsScf = 24`,  the accuracy `accuracyScf = 1.0e-6` in order to terminate the SCF computations and that the individual orbitals are improved due to the standard subshell order, `shellSequenceScf = Subshell[]`.

For the configuration-interaction (CI) parameters, the treatment of the Breit and QED interaction is of particular interest. At present, the defaults does not included neither Breit interactions nor QED. Such QED estimates can either be neglected (NoneQed()) or estimated by using an effective Hamiltonian approach due to Shabaev and coworkers (QedPetersburg()) or effective potential approach (QedSydney(); Flambaum et al.) However, further tests need to be done to better understand the reliability of these QED estimates to the level structure and state represetation of the ASF.

As seen from the settings above, moreover, there are special features in order to select individual levels for the CI computations, either in terms of their (relative) level No within the given multiplet or in terms of their level symmetry, i.e. their total angular momentum and parity, respectively. The `levelSelectionCI::LevelSelection` hereby tells whether (and which) selections were made; apparently, no selection of level numbers of symmetries is made by default though this can be overwritten. The selection of individual symmetries, in particular, may considerably reduced the computational effort as the Hamiltonian matrix need then to be calculated and diagonalized only for the selected symmetries.

In principle, these standard settings can be easily re-defined within the code by simply modifying the constructor AsfSetings() with no additional arguments. Alternatively, we can easily overwrite those parameters in some given (instance of) AsfSetting which we just wish to modify. This is achieved by
"""

# ╔═╡ a9b71283-8fea-4631-9e23-8b164bae31fd
asfSettings = AsfSettings(defaultAsfSettings; generateScf=true, jjLS=LSjjSettings(true), levelSelectionCI=LevelSelection(true, symmetries= [LevelSymmetry(0,"+"), LevelSymmetry(1,"-")]) )

# ╔═╡ fa6e2add-1be1-47be-87b1-12cdea6b602d
md"""
Finally, we can compute also a slightly enlarged level structure of C``^{2+}`` by including two futher configuration, for admixture and the enlargement of the level basis:
"""

# ╔═╡ e331a05c-2eed-4506-a055-9e1746449711
begin
	configs    = [Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2p^2"), Configuration("1s^2 3s^2"),  Configuration("1s^2 3p^2")]
	multiplet3 = SelfConsistent.performSCF(configs,  nucModel, grid, asfSettings)
end

# ╔═╡ 1cee605d-f846-41bf-970c-1816ed2e2d9d
md"""
We finish this (simple) tutorial by enlarging the configuration basis for the low-lying levels but by restricting the CI computations to the level symmetries **``J^{P}=0^{+}``** and **``1^{+}``**. This is achieved by specifying the settings to: 
"""

# ╔═╡ d70082f6-fa88-4b54-8306-46efe25138fc
md"""
# Still under construction !!!
"""

# ╔═╡ Cell order:
# ╟─79374db5-1fbd-4499-b8d9-c7b840bdee52
# ╟─8d491602-2e68-4ce9-b949-859f91cfee9d
# ╟─75061d53-4fb7-45a6-b8bb-c199d52479da
# ╟─f18a9488-13b8-11f0-342b-f99afaf3ac19
# ╠═fd4ef417-311f-43aa-a601-bf8a9d66b890
# ╟─956d7bd6-d59f-4daa-88b6-9b8879918d79
# ╠═5b8c341c-4923-4146-bdb1-5b6a7cc20c17
# ╟─1945551c-c5e1-4a75-897d-ef07205bb341
# ╠═561451e9-62b6-4a02-9d33-1393c575e890
# ╠═9787ed44-f14f-45a7-b2df-a40c9c01679a
# ╟─081bb1f4-fe80-4428-a234-e27970ae33be
# ╠═0e57775c-b714-4795-9630-2e82da20871b
# ╟─cd37b365-a6e7-4328-a513-7f9c29d4644e
# ╠═a9b71283-8fea-4631-9e23-8b164bae31fd
# ╟─fa6e2add-1be1-47be-87b1-12cdea6b602d
# ╠═e331a05c-2eed-4506-a055-9e1746449711
# ╟─1cee605d-f846-41bf-970c-1816ed2e2d9d
# ╟─d70082f6-fa88-4b54-8306-46efe25138fc
