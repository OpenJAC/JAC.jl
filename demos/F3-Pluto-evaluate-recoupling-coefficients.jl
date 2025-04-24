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
using JAC, SymEngine

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
# Evaluate re-coupling coefficients
"""

# ╔═╡ 75061d53-4fb7-45a6-b8bb-c199d52479da
md"""
Let us first invoke JAC, either from the own source-code basis or simply by **using JAC**
"""

# ╔═╡ 956d7bd6-d59f-4daa-88b6-9b8879918d79
md"""
**Note:** The Julia package `SymEngine` is needed to perform symbolic simplifications of Racah algebra expressions in JAC but, by default, is not automatically loaded.

In various research areas, the quantum mechanical description of many-particle structures and processes often requires an explicit transformation of the angular momenta (of the subsystems) due to different **coupling schemes**. Here, a quite simple example refers to the transformation (re-coupling) of three angular momenta  $j_1, j_2, j_3$. For this example, we already saw that the standard Clebsch-Gordan expansions may give rapidly rise to complex and cumbersome expression, and which are very prone for making errors. In general, many of these transformations can be expressed in terms of recoupling coefficients, a formal generalization of the well-known *Clebsch-Gordan* coefficients. Often, these recoupling coefficients need to be evaluated over and over again. Here, we introduce and explain a notation which makes the application and evaluation of general **recoupling coefficients** much easier. As always, we will need some Basic variables:
"""

# ╔═╡ 0b6a4fba-1d54-4c5b-9345-84aa32b3df70
begin
	j1  = Basic(:j1);    j2  = Basic(:j2);    j3  = Basic(:j3);    j4  = Basic(:j4);    j5  = Basic(:j5)
	j6  = Basic(:j6);    j7  = Basic(:j7);    j8  = Basic(:j8);    j9  = Basic(:j9);    j10 = Basic(:j10)
	j11 = Basic(:j11);   j12 = Basic(:j12);   j13 = Basic(:j13);   j14 = Basic(:j14);   j15 = Basic(:j15)
	j16 = Basic(:j16);   j17 = Basic(:j17);   j18 = Basic(:j18);   j19 = Basic(:j19);   j20 = Basic(:j20)
	j21 = Basic(:j21);   j22 = Basic(:j22);   j23 = Basic(:j23);   j24 = Basic(:j24);   j25 = Basic(:j25)
	J   = Basic(:J)
end

# ╔═╡ c04a5888-42b4-4bfb-9bb4-3c534c48ceeb
md"""
Let us consider again the recoupling coefficients $< (j_1, j_2) J_{12}, j_3: JM| j_1, (j_2,j_3) J_{23}: JM >$ for the re-coupling of three angular momenta. To avoid the explicit use of repeated Clebsch-Gordan expansions, we here introduce the notation of a coupling sequence, and which enables us to enter the coupling of each side of the coefficient separately. For this, we implemented the (data) type `RacahAlgebra.Csq`, search in *LiveDocs*. The struct RacahAlgebra.Csq then helps to express the left- and right-hand side of the recoupling coefficient above as:
"""

# ╔═╡ b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
begin
	leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( j1, j2, j12), j3, J)
	rightCsq = RacahAlgebra.Csq( j1, RacahAlgebra.Csq( j2, j3, j23), J)
end

# ╔═╡ 4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
md"""
and, obviously, could be easily extended towards much more complex coupling sequences. -- As before, we can evaluate this re-coupling coefficient by
"""

# ╔═╡ a7347554-9d1b-456d-a5f4-559f1d560f48
RacahAlgebra.evaluate(leftCsq, rightCsq)

# ╔═╡ 151c2037-f583-450b-a816-d5a0cfa12801
md"""
What need to be done next to make the result more obvious ??
"""

# ╔═╡ cfd7ccf9-36bd-4049-b93f-6b4bb7b47f31
begin
	leftCsq2  = RacahAlgebra.Csq( RacahAlgebra.Csq(j1,j2,j5), RacahAlgebra.Csq(j3,j4,j6), j7 )
	rightCsq2 = RacahAlgebra.Csq( j1, RacahAlgebra.Csq( RacahAlgebra.Csq(j2,j3,j8), j4, j9), j7 )
	rex2      = RacahAlgebra.evaluate(leftCsq2, rightCsq2)
end

# ╔═╡ Cell order:
# ╟─79374db5-1fbd-4499-b8d9-c7b840bdee52
# ╟─8d491602-2e68-4ce9-b949-859f91cfee9d
# ╟─75061d53-4fb7-45a6-b8bb-c199d52479da
# ╟─f18a9488-13b8-11f0-342b-f99afaf3ac19
# ╠═fd4ef417-311f-43aa-a601-bf8a9d66b890
# ╟─956d7bd6-d59f-4daa-88b6-9b8879918d79
# ╠═0b6a4fba-1d54-4c5b-9345-84aa32b3df70
# ╟─c04a5888-42b4-4bfb-9bb4-3c534c48ceeb
# ╠═b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
# ╟─4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
# ╠═a7347554-9d1b-456d-a5f4-559f1d560f48
# ╟─151c2037-f583-450b-a816-d5a0cfa12801
# ╠═cfd7ccf9-36bd-4049-b93f-6b4bb7b47f31
