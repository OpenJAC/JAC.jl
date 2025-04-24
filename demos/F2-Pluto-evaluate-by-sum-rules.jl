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
# Simplify Racah expressions by means of sum rules
"""

# ╔═╡ 75061d53-4fb7-45a6-b8bb-c199d52479da
md"""
Let us first invoke JAC, either from the own source-code basis or simply by **using JAC**
"""

# ╔═╡ 956d7bd6-d59f-4daa-88b6-9b8879918d79
md"""
**Note:** The Julia package `SymEngine` is needed to perform symbolic simplifications of Racah algebra expressions in JAC but, by default, is not automatically loaded.

As mentioned before, the (data) type `RacahExpression` is very central to applying the techniques from Racah's algebra; see *LiveDocs*. This data type enables one to comprise -- less or more -- sophisticated expressions into a single Julia variable and to attempt its simplification by a set of internal (sum) rules. As seen from the definition of this (data) struct, the different 'delta' and Wigner symbols of such an expression are kept and maintained separately, so that the known sum rules can be applied more readily. Moreover, the last constructor of a `RacahExpression` shows that it quite simple to overwrite or extent an already existing RacahExpression, starting from a simple 𝟙:
"""

# ╔═╡ b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
RacahExpression()

# ╔═╡ 4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
md"""
Perhaps, the simplest sum rules refer to the orthogonality of the Wigner 3-j and 6-j symbols; for example, the Wigner 6-j symbols fullfill the following orthogonality relations which can be displayed (just for illustration here) by looking up `RacahAlgebra.sumRulesForTwoW6j`. Both of the shown rules show that a summation over the product of two Wigner 6-j symbols can be re-written just in terms of some quantum numbers and triangle conditions. Note that $[a,b,...] = (2a+1)\: (2b+1)\: ...$. More general, all the implemented sum rules are displayed as inline comments in the code, although not as docstrings (apart from this particular function here).

*Typically, only some standard form of each sum rule is shown in the literature*, and many of these sum rules are just displayed in quite specialized texts about angular momenta. Likely, the most comprehensive compilation of these (and many other) rules can be found in the **monograph by Varshalovich et al. (1988)**. --- In general, however, one needs to recognize all the symmetries of a Racah expressions, implying all the phases and possible (weight) factors that arise from these symmetries. In JAC, this is realized by cycling automatically through all symmetric forms of the Wigner n-j (n = 3,6,9) symbols. In a later step, we also plan to take the spherical harmonics and the Wigner rotation matrices into account as well into the internal representation of a RacahExpression.

Again, let us first declare some Basic variables which we can later apply to define our first `RacahExpression`:
"""

# ╔═╡ e9325c21-5136-49bd-b9e4-77924560b4e2
begin
	a = Basic(:a);    b = Basic(:b);    c = Basic(:c);    d  = Basic(:d);    ee  = Basic(:ee);    f  = Basic(:f)
	g = Basic(:g);    h = Basic(:h);    k = Basic(:k);    l  = Basic(:l);    p   = Basic(:p);     q  = Basic(:q);     
	r = Basic(:r);    s = Basic(:s);    t = Basic(:t);    X  = Basic(:X);    Y   = Basic(:Y);     Z  = Basic(:Z);
end

# ╔═╡ 571405c1-0f1f-4cb5-88dd-31e4f9afc74d
begin
	aw6j = W6j(X, Y, Z, a, b ,c);    bw6j = W6j(X, Y, Z, a, b ,c)
	rex  = RacahExpression( [X, Y, Z], Integral[], Basic(0), Basic((2*X+1)*(2*Y+1)*(2*Z+1)), 
                            Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j], W9j[], Ylm[], Djpq[] )
end

# ╔═╡ b0619394-8dc7-43d0-b419-387eaaf9fe6b
md"""
As before, we can simply evaluate this expression which attempts to apply one of the -- more than 45 implemented -- sum rules in order to reduce either the number of Wigner symbols and/or the number of summation indices:
"""

# ╔═╡ d557472b-3056-4b6b-830c-dcfcbba34a27
RacahAlgebra.evaluate(rex)

# ╔═╡ 5b12989c-8069-4f42-96af-963881db0df2
md"""
This example looks perhaps quite *over-simplified* as we could use exactly the *orthogonaly relation* from above to get this result. However, the same simplification also works if we first randomly re-write the given Racah expression and then attempt its simplification again.
"""

# ╔═╡ b47294f2-626b-48c2-80fb-e326ea88d692
begin
	rex2 = RacahAlgebra.equivalentForm(rex);   @show rex2
	RacahAlgebra.evaluate(rex2)
end

# ╔═╡ b6575a7f-43a2-4bc6-ae0c-efa6e182eb20
md"""
You may test this simplification several times for (randomly) different equivalent forms of `rex` and, likely, will receive slightly different results with regard to the number of symbols and summations. This is related to the **phase issue**, which refers to the fact that it is not easy to always recognize how the overall phase can be re-written internally so that a particular sum rule applies. Here, we note that the application of any sum rule always requests that all other parts of the given Racah expression, including its overall phase, must be independent of those parts which are to be removed. Further (formal) improvement on this **phase issue** might be possible but, sometimes, these equivalences need to be recognized and corrected manually.

Of course, we can simplify also less obvious Racah expressions, such as:
"""

# ╔═╡ 4e7dbd4d-28c7-47c2-befb-d756e1c413ed
begin
	cw6j = W6j(a, b, X, c, d, p);    dw6j = W6j(c, d, X, b, a, q)
	rex3 = RacahExpression( [X], Integral[], Basic(X), Basic(2*X+1), Kronecker[], Triangle[], W3j[], W6j[cw6j, dw6j], W9j[], Ylm[], Djpq[] )
	@show rex3
	rex4  = RacahAlgebra.equivalentForm(rex3)
	RacahAlgebra.evaluate(rex4)
end

# ╔═╡ 6b391ff1-7605-4679-9673-5229c07562a4
begin
	ew6j = W6j(a, f, X, ee, b, s);   fw9j = W9j(a, f, X, d, q, ee, p, c, b)
	rex5 = RacahExpression( [X], Integral[], Basic(0), Basic(2*X+1), Kronecker[], Triangle[], W3j[], W6j[ew6j], W9j[fw9j], Ylm[], Djpq[] )
	@show rex
	rex6 = RacahAlgebra.equivalentForm(rex5)
	RacahAlgebra.evaluate(rex6)
end

# ╔═╡ 76823179-81b9-4c02-afbe-f415be01996d
md"""
Apart from these quite simple expressions, much more complex ones rapidly arise if angular momenta are coupled together or re-coupled in order to allow the simplification of many-particle matrix elements.

In the next example, we shall consider the (so-called) re-coupling coefficients $ < (j_1, j_2) J_{12}, j_3: JM| j_1, (j_2,j_3) J_{23}: JM >$ which is known to be independent of $M$. The expression of this re-coupling coefficient can be written down quite easily by applying twice a Clesch-Gordan expansion on both sides of the 'overlap matrix element'. Simple manipulations gives immediately rise to the `RacahExpression`:
"""

# ╔═╡ 998e3082-bb78-4add-a564-8cde78e136e9
begin
	j1   = Basic(:j1);    j2 = Basic(:j2);    j3 = Basic(:j3);    J12 = Basic(:J12);    J23 = Basic(:J23);    J = Basic(:J)
	m1   = Basic(:m1);    m2 = Basic(:m2);    m3 = Basic(:m3);    M12 = Basic(:M12);    M23 = Basic(:M23);    M = Basic(:M)
	w3ja = W3j(J12, j3, J, M12, m3, -M);        w3jb = W3j(j1, j2, J12, m1, m2, -M12)       
	w3jc = W3j(j2, j3, J23, m2, m3, -M23);      w3jd = W3j(j1, J23, J, m1, M23, -M)   
end

# ╔═╡ 85b4f8e2-38c4-43e2-b97e-a9c8e2050479
begin
	rex7 = RacahExpression( [m1, m2, m3, M12, M23], Integral[], -J12 + 2*j3 - 2*M - 2*j1 - M12 - M23 + J23, 
          	  (2*J+1) * sqrt( (2*J12+1)*(2*J23+1) ), Kronecker[], Triangle[], [w3ja, w3jb, w3jc, w3jd], W6j[], W9j[], Ylm[], Djpq[] )

end

# ╔═╡ 7eab1e4b-ab7d-4012-9b5d-9f96d2d989ea
md"""
which includes a five-fold summation (three further summations, for instance, for $m_1', m_2', m_3'$, and can be simplified because of the assumed normalization of the $|j_p m_p > $ states of all subsystems). Here, the simplification of this RacahExpression is already harder to see but can be obtained by calling:
"""

# ╔═╡ 5bfe0fea-a3a9-4ca8-8e8e-33b3f6048e27
RacahAlgebra.evaluate(rex7)

# ╔═╡ f9f9639f-68e0-4ee5-8cb4-9328ac4d76b3
md"""
The given recoupling coefficient is obviously independent of $M$ and just given by a Wigner 6j symbol times some rather trivial (delta) factors. A closer inspection of the Wigner symbol also enables one to express the phase in a slightly more convinient form.
"""

# ╔═╡ 759da271-0adf-4ae7-9642-8620c4c34690


# ╔═╡ Cell order:
# ╟─79374db5-1fbd-4499-b8d9-c7b840bdee52
# ╟─8d491602-2e68-4ce9-b949-859f91cfee9d
# ╟─75061d53-4fb7-45a6-b8bb-c199d52479da
# ╟─f18a9488-13b8-11f0-342b-f99afaf3ac19
# ╠═fd4ef417-311f-43aa-a601-bf8a9d66b890
# ╟─956d7bd6-d59f-4daa-88b6-9b8879918d79
# ╠═b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
# ╟─4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
# ╠═e9325c21-5136-49bd-b9e4-77924560b4e2
# ╠═571405c1-0f1f-4cb5-88dd-31e4f9afc74d
# ╟─b0619394-8dc7-43d0-b419-387eaaf9fe6b
# ╠═d557472b-3056-4b6b-830c-dcfcbba34a27
# ╟─5b12989c-8069-4f42-96af-963881db0df2
# ╠═b47294f2-626b-48c2-80fb-e326ea88d692
# ╟─b6575a7f-43a2-4bc6-ae0c-efa6e182eb20
# ╠═4e7dbd4d-28c7-47c2-befb-d756e1c413ed
# ╠═6b391ff1-7605-4679-9673-5229c07562a4
# ╟─76823179-81b9-4c02-afbe-f415be01996d
# ╠═998e3082-bb78-4add-a564-8cde78e136e9
# ╠═85b4f8e2-38c4-43e2-b97e-a9c8e2050479
# ╟─7eab1e4b-ab7d-4012-9b5d-9f96d2d989ea
# ╠═5bfe0fea-a3a9-4ca8-8e8e-33b3f6048e27
# ╟─f9f9639f-68e0-4ee5-8cb4-9328ac4d76b3
# ╠═759da271-0adf-4ae7-9642-8620c4c34690
