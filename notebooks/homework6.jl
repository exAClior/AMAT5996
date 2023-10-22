### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ b07f1119-5968-4368-9da4-d549555901cc
begin
    using Pkg
    Pkg.activate("..")
    using Symbolics, LinearAlgebra, Plots
end

# ╔═╡ cf105f34-7086-11ee-1031-01e4aa101f91
md"
# Problem 1
In this problem we will see how to implement the
algorithm mentioned in the paper:
[A Tutorial on Matrix Perturbation Theory](https://arxiv.org/abs/2002.05001v2).
Notations will follow the paper.
We use `Symbolics.jl` to help us document solving process.
Let's first define the variables we will use.
"

# ╔═╡ 01d7dcb5-dfef-4647-83e9-3a5bdba57637
begin
@variables E₁ E₂ aᵣ aᵢ bᵣ bᵢ λ;
a = aᵣ + im * aᵢ;
b = bᵣ + im * bᵢ;
nothing
end

# ╔═╡ 831e47f6-5c36-4ed9-b0be-7a4e2fa6605c
md"
Our Hamiltonian could be divided into two parts.
$A_0$ represents the Hamiltonian before perturbation.
$A_1$ represents the perturbation.
"

# ╔═╡ e2d00203-32d2-4cc3-95c3-2ec62ff4610d
begin
A₀ = Matrix([E₁ 0 0; 0 E₁ 0; 0 0 E₂])
A₁ = Matrix([0 0 a; 0 0 b; conj(a) conj(b) 0])
Aₑ = A₀ + A₁
end


# ╔═╡ 212bfe72-1510-4c0a-b3fa-4d13774787cd
md"
Let's first try to solve the eigenvalues exactly.
We will use `det` to get the characteristic equation.
Then we notice that $E₁$ is a root of this equation.
We then divide $λ - E₁$ from the equation to get a lower order polynomial which
could be solved by quadratic formula.
"

# ╔═╡ 5a03812e-5a26-4f37-9f0f-5a264b31d3a3
    chara_eqn = det(Aₑ - λ .* Matrix(I, 3, 3))

# ╔═╡ 17c9e2a7-78b1-4079-a518-c533c7bff10a
    substitute(chara_eqn, Dict(λ => E₁))

# ╔═╡ 962a75c1-2ac5-42df-acf5-596c498bc92c
 reduced_chara_eqn = -1 * simplify(chara_eqn / (λ - E₁))

# ╔═╡ 4ae212d5-00bc-4f8a-acf8-1dbbd9a45f98
begin
    const_term = substitute(reduced_chara_eqn, Dict(λ => 0))
    linear_term = substitute(simplify((reduced_chara_eqn - const_term) / λ), Dict(λ => 0))
    quadratic_term = simplify((reduced_chara_eqn - const_term - linear_term * λ) / λ^2)
    nothing
end

# ╔═╡ b786c77f-deaa-49ec-b3f8-95553d1f01cb
b2m4ac = real(simplify(linear_term^2 - 4 * quadratic_term * const_term))

# ╔═╡ 85f56594-4c27-4b3c-bcf0-e4915a6a6704
root2 = simplify((-linear_term + sqrt(b2m4ac)) / (2 * quadratic_term))

# ╔═╡ 85b5922d-979f-40ec-b6c4-db179e3f1526
root3 = simplify((-linear_term - sqrt(b2m4ac)) / (2 * quadratic_term))

# ╔═╡ a9a93cd0-f21a-4382-93a1-c64c385cb36c
md"
In summary, we have the following three eigenvalues of the perturbed Hamiltonian.
1. $\lambda_1 = E_1$
2. $\lambda_2 = \frac{1}{2} \left( E_1 + E_2 + \sqrt{\left(  E_1 - E_2 \right)^{2} + 4 \left(  a_i^{2} + a_r^{2} + b_i^{2} + b_r^{2} \right)} \right)$
3. $\lambda_3 = \frac{1}{2} \left( E_1 + E_2 - \sqrt{\left(  E_1 - E_2 \right)^{2} + 4 \left(  a_i^{2} + a_r^{2} + b_i^{2} + b_r^{2} \right)} \right)$
Note, we could apply the approximation $\sqrt{1+x} \approx 1 + x/2$ when $x \ll 1$. And get
1. $\lambda_1 = E_1$
2. $\lambda_2 \approx  \left( E_1 - \frac{ \left(  a_i^{2} + a_r^{2} + b_i^{2} + b_r^{2} \right)}{(E_2 - E_1)} \right)$
3. $\lambda_3 \approx  \left( E_2 + \frac{ \left(  a_i^{2} + a_r^{2} + b_i^{2} + b_r^{2} \right)}{(E_2 - E_1)} \right)$

Let's try to solve for the eigenvalues using perturbation approach.
Firstly, let's blindly apply the non-degenerate perturbation method to second order.
Let's do first order perturbation. We know the Hamiltonian before perturbation has eigenvalues E₁, E₁, E₂
and eigenvectors [1, 0, 0], [0, 1, 0], [0, 0, 1] respectively.
"

# ╔═╡ af4d6281-d0a8-4d6f-9421-62644c5a14d3
A₀

# ╔═╡ d83e1e50-da8c-4702-aaa7-8e7cce8bdc90
md"
We can define $V_0$ and $\Lambda_0$ accordingly.
"

# ╔═╡ b7d95e93-12f7-4ac5-9c7f-01e67fa80af0
V₀ = Matrix([1 0 0; 0 1 0; 0 0 1])

# ╔═╡ 3ccafbb1-60c8-409d-a576-94821fec3ffb
Λ₀ = diagm([E₁,E₁,E₂])

# ╔═╡ d7f2e4df-cd5f-4b66-b097-3b12a8342e0b
md"
We could use eqn (11) from Bamieh to get first order perturbation to eigenvalues.

$\Lambda_1 = dg(V_0^* A_1 V_0)$
"

# ╔═╡ 42145ac1-715e-4f6e-86bc-87d1964ae165
Λ₁ = Diagonal(V₀' * A₁ * V₀)


# ╔═╡ 152ea087-984b-4e2e-a62f-45300256911c
md"
Unfortunately, we get nothing from the first order perturbation.
Let's go to the next level where we also allow the eigenstate to vary a little.

"

# ╔═╡ 09e5ac1e-60e6-48b6-b7fb-7e7015cd2dd4
Πdagcirc = Matrix([0 0 1/(E₁-E₂); 0 0 1/(E₁-E₂); 1/(E₂-E₁) 1/(E₂-E₁) 0])

# ╔═╡ 788a1d67-cf73-4d0e-bf09-28be4debcdd3
V₁ = -V₀ * (Πdagcirc .* (V₀ * A₁ * V₀))


# ╔═╡ bd97e7fa-0979-4218-9454-f9930b39e154
Λ₂ = simplify.(Diagonal(V₀ * A₁ * V₁))

# ╔═╡ 30775d34-c171-4639-8d28-e56014b00cd6
md"
Comparing this to our exact result above, we can see something is wrong.
For correction of first and second eigenstate, we got it wrong.
The energy correction term should be

1. $\delta\lambda_1 = 0$
2. $\delta\lambda_2 = -\frac{ \left(  a_i^{2} + a_r^{2} + b_i^{2} + b_r^{2} \right)}{(E_2 - E_1)}$
3. $\delta\lambda_3 = \frac{ \left(  a_i^{2} + a_r^{2} + b_i^{2} + b_r^{2} \right)}{(E_2 - E_1)}$
"

# ╔═╡ 8443f122-db68-4807-a21d-5970335009c3


# ╔═╡ 69f6c0fa-b5e4-45e1-8e5c-5e5b5cfd5bec


# ╔═╡ cdbde18f-c4d4-494e-b953-851fa07edc74


# ╔═╡ e0ae1294-6efb-4b4a-9e20-292064c940cf


# ╔═╡ Cell order:
# ╠═b07f1119-5968-4368-9da4-d549555901cc
# ╟─cf105f34-7086-11ee-1031-01e4aa101f91
# ╠═01d7dcb5-dfef-4647-83e9-3a5bdba57637
# ╟─831e47f6-5c36-4ed9-b0be-7a4e2fa6605c
# ╠═e2d00203-32d2-4cc3-95c3-2ec62ff4610d
# ╟─212bfe72-1510-4c0a-b3fa-4d13774787cd
# ╠═5a03812e-5a26-4f37-9f0f-5a264b31d3a3
# ╠═17c9e2a7-78b1-4079-a518-c533c7bff10a
# ╠═962a75c1-2ac5-42df-acf5-596c498bc92c
# ╠═4ae212d5-00bc-4f8a-acf8-1dbbd9a45f98
# ╠═b786c77f-deaa-49ec-b3f8-95553d1f01cb
# ╠═85f56594-4c27-4b3c-bcf0-e4915a6a6704
# ╠═85b5922d-979f-40ec-b6c4-db179e3f1526
# ╟─a9a93cd0-f21a-4382-93a1-c64c385cb36c
# ╠═af4d6281-d0a8-4d6f-9421-62644c5a14d3
# ╟─d83e1e50-da8c-4702-aaa7-8e7cce8bdc90
# ╠═b7d95e93-12f7-4ac5-9c7f-01e67fa80af0
# ╠═3ccafbb1-60c8-409d-a576-94821fec3ffb
# ╟─d7f2e4df-cd5f-4b66-b097-3b12a8342e0b
# ╠═42145ac1-715e-4f6e-86bc-87d1964ae165
# ╟─152ea087-984b-4e2e-a62f-45300256911c
# ╠═09e5ac1e-60e6-48b6-b7fb-7e7015cd2dd4
# ╠═788a1d67-cf73-4d0e-bf09-28be4debcdd3
# ╠═bd97e7fa-0979-4218-9454-f9930b39e154
# ╟─30775d34-c171-4639-8d28-e56014b00cd6
# ╠═8443f122-db68-4807-a21d-5970335009c3
# ╠═69f6c0fa-b5e4-45e1-8e5c-5e5b5cfd5bec
# ╠═cdbde18f-c4d4-494e-b953-851fa07edc74
# ╠═e0ae1294-6efb-4b4a-9e20-292064c940cf
