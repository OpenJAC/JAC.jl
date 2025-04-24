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
# Evaluate Wigner n-j symbols by recursions & special values
"""

# ╔═╡ 75061d53-4fb7-45a6-b8bb-c199d52479da
md"""
Let us first invoke JAC, either from the own source-code basis or simply by **using JAC**
"""

# ╔═╡ 956d7bd6-d59f-4daa-88b6-9b8879918d79
md"""
Note: The Julia package `SymEngine` is **needed** in order to symbolically manipulate and simplify Racah algebra expressions in JAC **but, by default, is not automatically loaded.**

In atomic and quantum many-particle physics, (the coupling of) angular momenta and spherical tensors play a very crucial role for the evaluation and spin-angular integration of many-particle matrix elements. This coupling often leads to algebraic expressions which, more often than not, are written in terms of generalized Clebsch-Gordan coefficients and/or Wigner nj symbols as well as the Wigner rotation matrices and spherical harmonics. Although the evaluation and simplification of such expressions is, at least in principle, a rather straigthforward task, it may become (extremely) cumbersome if complex systems or physical scenarios are considered.

After the pioneering work by Wigner in the late 1930s, in particular Guilio Racah developed a powerful machinery in the 1940s, known as **Racah algebra (techniques)** today, in order to deal with the angular momenta of (rotationally invariant) quantum many-particle systems. Briefly speaking, this machinery includes a number of strategies to simplify (so-called) **Racah expressions**. In JAC, these expressions are defined by a summation over Wigner nj symbols of different kind(s) as well as (various integrals over) the spherical harmonics and Kronecker and triangular deltas, cf. User Guide, Section 15.1. Here, the summation may formally run over any number of angular momenta and Wigner nj symbols. Since, in the theory of angular momentum, most symbols obey a rather high symmetry, the complexity of such Racah expressions increases very rapidly as more Wigner symbols, rotation matrices and/or spherical harmonics appear in the (product) terms over which one needs to sum.

There are different strategies, that can be (successively) applied in order to **simplify such Racah expressions algebracially**, i.e. if the angular momenta are not specified numerically. These strategies include:

* Use of known **special values**. This strategy replaces single Wigner nj symbols or spherical harmomic by a (much) simpler expression that, in particular, does not contain any implicit summation. In practice, each Wigner nj symbol can be analysed and perhaps replaced independently by some special-value rule, if available and if useful in the given context.
* Use of **orthogonality relations** for the Wigner nj symbols or spherical harmonics.
* Use of known algebraic **sum rules** for the Wigner symbols; in fact, the orthogonality relations can be treated as particular sum rules, and this is done also in JAC.

To make use of these strategies, we need a suitable **representation** of the various symbols from the theory of angular momenta as well as for Racah expressions as a whole. Therefore, let us first have a look for the (internal) representation of the Wigner nj symbol `W3j(ja,jb,jc,ma,mb,mc)`; search for `W3j` and `W9j` in *LiveDocs*.

"""

# ╔═╡ b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
md"""
However, the real power of Racah's algebra become apparent only, if products of these Wigner symbols appear in some summation, and if multiple summations/integrations can be performed algebraically. To deal with such (product) expressions or a sum of such expressions, JAC defines a (more complex struct) for a `RacahExpression` as a whole.

As already seen from the arguments, such expressions enables one to comprise lists (arrays) of Wigner symbols of different types, an overall phase, weight and the formal summation over quantum numbers. Here, we shall postpone the further discussion of these `RacahExpression` but will need them to just handle the output of the recursion relations or the symmetry-representations of the Wigner symbols.

In order to make use of JAC's extension for simplifying expressions from Racah' algebra, we shall need a number of symbolic variables, whose type `Basic` is derived from `SymEngine` (although not much in-line documentation is provided by this package); in this and the two subsequent tutorials, we shall define all necessary Basic variables nearby to where we will need them. We define a number of such (Basic) variables to facilitate our later discussion:
"""

# ╔═╡ 72bbbb7f-2c43-4ffb-b0b6-458759430612
begin
	a  = Basic(:a);     b  = Basic(:b);    c  = Basic(:c);   d  = Basic(:d);     ee  = Basic(:ee);    f  = Basic(:f)
	X  = Basic(:X);     Y  = Basic(:Y);    Z  = Basic(:Z)  
	j  = Basic(:j);     m  = Basic(:m);    typeof(m)
end

# ╔═╡ 3544b096-cc01-4836-a932-3c438bc0e2b6
md"""
The Wigner nj symbols are known to fullfill certain recursion relations which enable one to step-up or step-down in the quantum numbers. These recursion relations have been applied to make the underlying physics more obvious or to compute arrays of Wigner symbols more efficiently. For the Wigner 3j symbols, especially, different recursions are known and are distinguished in JAC by the (abstract) type `RacahAlgebra.AbstractRecursionW3j`. We can apply these recusions to any given Wigner 3j symbol:
"""

# ╔═╡ 4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
begin
	ja  = Basic(:ja);    jb = Basic(:jb);    jc = Basic(:jc);    ma = Basic(:ma);    mb = Basic(:mb);    mc = Basic(:mc)
	w3j = W3j(ja, jb, jc, ma, mb, mc)
end

# ╔═╡ 759da271-0adf-4ae7-9642-8620c4c34690
md"""
by just typing:
"""

# ╔═╡ a350cff5-5cb7-4571-8923-c734faabc63f
RacahAlgebra.recursionW3j(w3j, RacahAlgebra.RecursionW3jMagnetic())

# ╔═╡ c3fb5b25-ef6b-4049-a759-c59ac81ac9ba
md"""
In general, all recursion relations give rise to an array of `RacahExpressions`, though without any additional summation or other Wigner symbols; of course, we can apply any other recursion as well:
"""

# ╔═╡ f020366b-a4d7-45c0-b296-d9ec20c2fec7
RacahAlgebra.recursionW3j(w3j, RacahAlgebra.RecursionW3jOneStep())

# ╔═╡ 8465236f-f6e7-4e27-a223-1397a4f8cd17
RacahAlgebra.recursionW3j(w3j, RacahAlgebra.RecursionW3jLouck())

# ╔═╡ 6545c3c3-88db-49fe-8ee7-f7daf8c45ce9
md"""
If not all angular momenta (quantum numbers) of a given Wigner nj symbol are independent of each other (but defined by some fixed relation, for instance, $b = a + 1$), one often attempts to re-write the Wigner symbol in terms of some **special value**. In depends on the particular application of how useful such a (formal) simplification is in practice. In JAC, the *search for such special values* is facilitated by calling `RacahAlgebra.evaluate()`:
"""

# ╔═╡ 4dfa610a-8fe3-477b-8571-7253aeeff3c0
begin
	w3ja = W3j(j+3//2, j, 3//2, m, -m-3//2, 3//2)
	wa  = RacahAlgebra.evaluate(w3ja)
end

# ╔═╡ 2985ca4d-579c-4920-9289-d9796069ff64
md"""
The last step appears to be quite simple, as the special value of `W3j(3/2 + j, j, 3/2; m, -3/2 - m, 3/2)` can be found in many texts. Owing to the symmetry of the Wigner nj symbols, however, these special values are not always simple to recognize. We here note that the Wigner 3-j symbols have 12 classical equivalent forms (symmetries), the Wigner 6-j symbols already 24, and the Wigner 9-j symbols even 72 **equivalent representations (symmetries)**, and without that their value would change if the proper phase is taken into account. These numbers of equivalent representation (symmetries) numbers do not yet account for the **Regge symmetries**, in which the quantum numbers may change by $\pm 1/2$. --- We can easily re-write all Wigner symbols by just choosing randomly any of the symmetric form (including the proper phase) and, then, try to simplify it again:
"""

# ╔═╡ 56c18b94-adc7-4174-b2b0-732cf51215a5
rex = RacahAlgebra.equivalentForm(w3ja);    RacahAlgebra.evaluate(rex, special=true)

# ╔═╡ 6a740b65-e885-443d-83a9-6593905dffb7
md"""
In contrast to the literature, where the special values are typically shown for one standard form of the Wigner symbols (as, for instance, in the initial definition above), JAC also finds the special value for all equivalent forms as it *internally cycles through all symmetries*, and by keeping the corresponding phase and weight factors (cf. later) into account. The result of such an evaluation is always a `RacahExpression`, though it is printed in a -- more or less -- neat format.

Apart from the classical symmetries, (so-called) **Regge symmetries** are known for the Wigner 3-j and 6-j symbols and could be readily implemented into JAC, if necessary. In practice, various special values of the Wigner nj symbols are recognized also by Mathematica and, perhaps, by a few other computer-algebra systems. We include them here for the sake of convinience, while the **main emphasis is clearly placed upon the sum rule evaluation, a very special and convenient feature of JAC.**

A few other examples for special values of the Wigner 3-j, 6-j and 9-j symbols are:
"""

# ╔═╡ 3eaf5441-fd0a-4212-aa58-846f354b4b7b
w3jb = W3j(j, j, 0, m, -m, 0);    RacahAlgebra.evaluate(w3jb)

# ╔═╡ 60ded5b7-c7e0-4efa-88e0-02eb991a87d0
md"""
Alternatively, we may call this evaluation as part of a RacahExpression as in the following example:
"""

# ╔═╡ d5fc16e4-1592-4de4-ae76-c191e35eb802
begin
	w6ja = W6j( a, b, c, 2, c-2, b-2)
	rexa = RacahAlgebra.equivalentForm(w6ja)
	RacahAlgebra.evaluate(rexa, special=true)
end

# ╔═╡ ab2d37dd-bf42-47db-9007-f79814262501
begin
	g   = Basic(:g);  h = Basic(:h)
	w9j = W9j(a, b, c, d, ee, f, g, h, 0)
	RacahAlgebra.evaluate(w9j)
end

# ╔═╡ 8037ce41-59a6-421d-9a9d-5638af13cf06
md"""
The last example show that a Wigner 9j symbol always simplifies to a 6j symbol if one of the quantum numbers is 0. In the next notebook we shall demonstrate how **sum rules** can be applied to more complex `RacahExpression`, the real power of *Racah's algebra*.
"""

# ╔═╡ 85cb1ef0-cc08-493e-a82e-57d849caf216
w6j = W6j( a, b, c, 2, c-2, b+2);   RacahAlgebra.evaluate(w6j)

# ╔═╡ 64394599-eeaa-4332-991f-ecae6ca9426b
# ╠═╡ disabled = true
#=╠═╡
begin
	w6j = W6j( a, b, c, 3//2, c-1//2, b+1//2)
	RacahAlgebra.evaluate(w6j)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─79374db5-1fbd-4499-b8d9-c7b840bdee52
# ╠═8d491602-2e68-4ce9-b949-859f91cfee9d
# ╟─75061d53-4fb7-45a6-b8bb-c199d52479da
# ╟─f18a9488-13b8-11f0-342b-f99afaf3ac19
# ╠═fd4ef417-311f-43aa-a601-bf8a9d66b890
# ╟─956d7bd6-d59f-4daa-88b6-9b8879918d79
# ╟─b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
# ╠═72bbbb7f-2c43-4ffb-b0b6-458759430612
# ╟─3544b096-cc01-4836-a932-3c438bc0e2b6
# ╠═4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
# ╟─759da271-0adf-4ae7-9642-8620c4c34690
# ╠═a350cff5-5cb7-4571-8923-c734faabc63f
# ╟─c3fb5b25-ef6b-4049-a759-c59ac81ac9ba
# ╠═f020366b-a4d7-45c0-b296-d9ec20c2fec7
# ╠═8465236f-f6e7-4e27-a223-1397a4f8cd17
# ╠═6545c3c3-88db-49fe-8ee7-f7daf8c45ce9
# ╠═4dfa610a-8fe3-477b-8571-7253aeeff3c0
# ╟─2985ca4d-579c-4920-9289-d9796069ff64
# ╠═56c18b94-adc7-4174-b2b0-732cf51215a5
# ╟─6a740b65-e885-443d-83a9-6593905dffb7
# ╠═3eaf5441-fd0a-4212-aa58-846f354b4b7b
# ╠═64394599-eeaa-4332-991f-ecae6ca9426b
# ╠═85cb1ef0-cc08-493e-a82e-57d849caf216
# ╟─60ded5b7-c7e0-4efa-88e0-02eb991a87d0
# ╠═d5fc16e4-1592-4de4-ae76-c191e35eb802
# ╠═ab2d37dd-bf42-47db-9007-f79814262501
# ╟─8037ce41-59a6-421d-9a9d-5638af13cf06
