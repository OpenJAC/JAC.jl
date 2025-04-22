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

# ╔═╡ 2d8989c9-f09d-4794-b352-1ee66e7d8764
begin
	Defaults.convertUnits("energy: from atomic", 1.0)
	Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
	using Printf
	Printf.@sprintf("%.4e", Defaults.convertUnits("energy: from atomic", 1.0))
end

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
# Welcome to JAC, the Jena Atomic Calculator
"""

# ╔═╡ 75061d53-4fb7-45a6-b8bb-c199d52479da
md"""
... that provides various tools for performing atomic (structure) calculations of different kinds and complexities. Apart from the computation of atomic (many-electron) amplitudes, properties and processes, JAC supports interactive, restricted-active space (RAS) and cascade computations. It also help perform a few simple hydrogenic and semi-empirical estimates as well as simplify symbolic expressions from Racah's algebra, see the documentation for further details.
"""

# ╔═╡ 956d7bd6-d59f-4daa-88b6-9b8879918d79
md"""
In the design of JAC, we first of all **aim for a precise language** that (i) is simple enough for both, seldom and a more frequent use of this package, (ii) highlights the underlying physics and (iii) avoids most technical slang that is often unnecessary but quite common to many other codes. An intuitive picture about the level or hyperfine structure of an atom, its properties as well as possible excitation and/or decay processes should (always) come first in order to generate the desired data: By making use of suitable data types (struct), **we indeed wish to introduce a language close to the underlying formalism.** --- While JAC is overall based on a rather large number of such types, a few simple examples are: \
	+ (atomic) Shell: ``1s, 2s, 2p, ...`` \
	+ Subshell: ``1s_{1/2}, 2s_{1/2}, 2p_{1/2}, 2p_{3/2}, ...`` \
	+ (electron) Configuration: ``1s^2 2s^2 2p^6 3s`` or ``[Ne] 3s, ...`` \
	+ Level: ``1s^2 2s^2\;\: ^1S_0``, ...  \

and many other terms (types) that we shall explain later. --- Let us simply start, for instance, with specifying and assigning the and shells:
"""

# ╔═╡ b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
begin
	w1s = Shell("1s")
	w2p = Shell("2p")
	Subshell("2p_1/2"),   Subshell("2p_3/2")  
end

# ╔═╡ 4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
md"""
In JAC, we make use of these `Shell`'s and `Subshell`'s whenever they will naturally occur in describing the level structure or the excitation, decay or occupation of an atom, and this both at input and output. If you have forgotten how to specify such a subshell (constructor), looks at the *LiveDocs* for them. --- Of course, we can interactively also specify any electron configuration:
"""

# ╔═╡ bbd3211e-8ff3-497e-960c-2c8ab04fd7f9
begin
	wc1 = Configuration("1s^2 2s^2 2p^5")
	wc2 = Configuration("[Ar] 4s^2 3d^5")
end

# ╔═╡ 759da271-0adf-4ae7-9642-8620c4c34690
md"""
This input just shows three (very) simple examples and how the details of some computation can be readily specified in line with our basic understanding of the atomic shell model.You can look at the documentation for a complete list of most data structures that are speficic to the JAC module. --- The use of a proper terminology and data structures has been found essential for developing the JAC module. Although we presently support just a (more or less small) number of frequently requested tasks in atomic structure and collision theory, we tried to define data types that are flexible enough to further extend these tools in the future. Following the Julia's standard conventions, all types (struct) are named in CamelCase notation. The power of these data types lays not in their number but in the consistency of how they are declared and utilized internally. It is the aim of the documentation as well as these (Pluto) notebooks to show the user how one can benefit from a proper set of such data types (`struct`'s) without that one needs to know all of them nor the details how they are defined. 
"""

# ╔═╡ 055f695a-1307-4a81-917a-f3efe9b1ab1f
md"""
**Use of physical units:** Of course, there are many other features that make Julia & JAC as powerful as it is: For example, the user may pre-define and overwrite the units in which he wishes to communicate with JAC. These units determine how (most of) the input data are interpreted as well as output data are displayed in tabulations or to screen. The current defaults settings for the units can be seen by typing:
"""

# ╔═╡ 13fe0a92-e138-4f16-9051-a80302b3d46a
Basics.display("settings")

# ╔═╡ 41599617-1b88-4815-90a5-79f967c9ea04
md"""
which show that energies are taken/printed in eV, rates in 1/s, etc. Apart from modifying these defaults directly in the source code, the can be overwritten by the user at any time of the program executation. This is done by means of the function `Defaults.setDefaults` which enables one to re-define various global values of JAC. If we wish to enter/display energies in Kaysers or cross sections in atomic units, we can simply type: 
"""

# ╔═╡ 06742b88-807b-4a7d-af3f-1cbc20660a3f
begin
	Defaults.setDefaults("unit: energy", "Kayser")
	Defaults.setDefaults("unit: cross section", "a.u.")
	Basics.display("settings")
end

# ╔═╡ 378d2c6c-2ebe-4ebe-86ea-0c27df79037d
md"""
Apart from the default units, one can similarly overwrite the method that is use for the generation and normalization of continuum orbitals and several others. Although called global, the corresponding values can be accesses just in two ways. (i) The global constants, such as the electron mass, the speed of light, the fine-structure constant, etc., are accessed via the function `Defaults.getDefaults` 
"""

# ╔═╡ d937ddee-3122-4b2d-8c10-1a76f4541ffa
begin
	Defaults.getDefaults("alpha")
	Defaults.getDefaults("electron rest energy")
	Defaults.getDefaults("unit: energy")
end

# ╔═╡ 92bbbefa-44da-4f2f-b548-7e828557107d
md"""
(ii) These global values are frequently applied in order to -- internally or externally -- convert physical numbers into units of the same dimension. This is done by the function `Defaults.convertUnits`. In JAC, the call of this function is often combined with some proper formatting of the results.
"""

# ╔═╡ Cell order:
# ╟─79374db5-1fbd-4499-b8d9-c7b840bdee52
# ╟─8d491602-2e68-4ce9-b949-859f91cfee9d
# ╠═75061d53-4fb7-45a6-b8bb-c199d52479da
# ╟─f18a9488-13b8-11f0-342b-f99afaf3ac19
# ╠═fd4ef417-311f-43aa-a601-bf8a9d66b890
# ╟─956d7bd6-d59f-4daa-88b6-9b8879918d79
# ╠═b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
# ╟─4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
# ╠═bbd3211e-8ff3-497e-960c-2c8ab04fd7f9
# ╟─759da271-0adf-4ae7-9642-8620c4c34690
# ╟─055f695a-1307-4a81-917a-f3efe9b1ab1f
# ╠═13fe0a92-e138-4f16-9051-a80302b3d46a
# ╟─41599617-1b88-4815-90a5-79f967c9ea04
# ╠═06742b88-807b-4a7d-af3f-1cbc20660a3f
# ╟─378d2c6c-2ebe-4ebe-86ea-0c27df79037d
# ╠═d937ddee-3122-4b2d-8c10-1a76f4541ffa
# ╟─92bbbefa-44da-4f2f-b548-7e828557107d
# ╠═2d8989c9-f09d-4794-b352-1ee66e7d8764
