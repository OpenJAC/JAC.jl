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
# Obtaining data from the periodic table
"""

# ╔═╡ 75061d53-4fb7-45a6-b8bb-c199d52479da
md"""
Let us first invoke JAC, either from the own source-code basis or simply by **using JAC**
"""

# ╔═╡ 956d7bd6-d59f-4daa-88b6-9b8879918d79
md"""
Most atomic computations require input that is specific for some given element or isotope. This information can often be easily read from the Periodic Table of Elements and have to be provided to JAC, wherever necessary. In other situations, it may be helpful to have these data also available within the same toolbox. To obtain for instance the nuclear charge of some element with known symbol, one can call the function `PeriodicTable.getAtomicNumber` [cf. *LiveDocs*] where symbols are defined in Julia by the syntax: 
"""

# ╔═╡ b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
[:He, :Ne, :Ar]

# ╔═╡ 4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
PeriodicTable.getAtomicNumber(:U)

# ╔═╡ 759da271-0adf-4ae7-9642-8620c4c34690
md"""
Similarly, basis information about the elements can be obtained from the function `PeriodicTable.getData`
"""

# ╔═╡ d5e827cf-ca43-4b6a-ad79-a1128e8dd187
begin
	@show PeriodicTable.getData("mass", :U), PeriodicTable.getData("mass", 92)
	@show PeriodicTable.getData("1st IP", :He), PeriodicTable.getData("1st IP", :Ne)
	@show PeriodicTable.getData("polarizibility", :He), PeriodicTable.getData("polarizibility", :Ne)
	@show PeriodicTable.getData("ground configuration", :He), PeriodicTable.getData("ground configuration", :Ne)
end

# ╔═╡ 469176ba-343e-462e-bf0c-892471c7e8e6
md"""
These data might be useful especially if some tabulations or legend of some plots need to be prepare. Apart from those data, which can be read of directly from various representations of the periodic table, we shall provide also some isotope-selected data, although no attempt will be made to compete here with any useful tabulations from the literature. For many cases, and as seen below, this function is simply not yet implemented properly, cf. `PeriodicTable.getIsotopeData` and `Semiempirical.estimate`.  
"""

# ╔═╡ c89e90c3-0d42-4e36-b8bc-75e5adb2d7e7


# ╔═╡ a460393c-4c10-473d-9c27-0e5d032e2f58


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
# ╠═b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
# ╠═4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
# ╟─759da271-0adf-4ae7-9642-8620c4c34690
# ╠═d5e827cf-ca43-4b6a-ad79-a1128e8dd187
# ╟─469176ba-343e-462e-bf0c-892471c7e8e6
# ╠═c89e90c3-0d42-4e36-b8bc-75e5adb2d7e7
# ╠═a460393c-4c10-473d-9c27-0e5d032e2f58
# ╟─d70082f6-fa88-4b54-8306-46efe25138fc
