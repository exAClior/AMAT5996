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

Let's try to apply degenerate energy perturbation theory.
We notice for vector spanned by \{ [1 0 0] , [0 1 0]\},
we have energy degeneracy for $A_0$.

According to the paper, we need to find a special pair of linear combinations of
the above mentioned basis vectors so that we diagonalize $A_1$ in this subspace.
However, we could not as $A_1$ is already diagnoalized and degenerate.
"

# ╔═╡ 524b736b-fbb7-4425-a55a-ee7a73bba776
A₁[1:2,1:2]

# ╔═╡ 8443f122-db68-4807-a21d-5970335009c3
md"
In this case, the paper does not tell us what to do.
I found [Gottfried Chapter 3.7](https://books.google.co.jp/books/about/Quantum_Mechanics_Fundamentals.html) to be extremely helpful.
Using (273) directly gives us the effective Hamiltonian of a lower dimension (the dimension of degenerate subspace).
Reading off the (272) you could immediately see that the eigenvalue of this effective hamiltonian is the energy change from original Hamiltonian.
"

# ╔═╡ 69f6c0fa-b5e4-45e1-8e5c-5e5b5cfd5bec
Heff = Matrix(simplify.([ A₁[i,3]*A₁[3,j]/(E₁ - E₂) for i in 1:2, j in 1:2 ]))

# ╔═╡ cdbde18f-c4d4-494e-b953-851fa07edc74
simplify(det(Heff - λ * Matrix(I,2,2)))

# ╔═╡ e0ae1294-6efb-4b4a-9e20-292064c940cf
md"
It's obvious that one solution is $λ = 0$ and the other is $λ = \frac{|a|^2 + |b|^2}{E_2 - E_1}$
This matches our exact calculuation above.
"

# ╔═╡ b97425a6-6538-4fc5-82de-26475275a1f8
md"
# Problem 2: Scaling Laws for Hydrogen-like Atoms

Recall a hydrogenlike atom is defined as a cation with only one valence shell electron.

For the frist few properties, we rely on the effectiveness of Bohr model.
For details of derivation, please refer to [Libre Text](https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Stellar_Atmospheres_(Tatum)/07%3A_Atomic_Spectroscopy/7.04%3A_The_Bohr_Model_of_Hydrogen-like_Atoms).
"

# ╔═╡ a4f9ac12-a595-46d7-b8cc-870c27c4fa9e
@variables n;

# ╔═╡ ea7a825f-a545-4527-975b-ef01dc92b5c5
function r_scaling(n::Num)
    return n^2
end

# ╔═╡ 0cc39cdf-70d3-4a09-83e1-2296bc3a62ab
function E_scaling(n::Num)
    return -1/n^2
end

# ╔═╡ d3daec1c-6972-4024-b3e0-8a1bdc384a1d
function ΔE_scaling(n::Num)
    return simplify(E_scaling(n+1) - E_scaling(n))
end

# ╔═╡ 0767c090-9177-42e5-bf61-3738b293bcfb
# ~ 1/n^3
ΔE_scaling(n)

# ╔═╡ 04e8d2fa-fe02-4baf-ac0d-c6c28bfd8aee
md"
The transition dipole moment of ground state to $\ket{nl}$ is
"

# ╔═╡ 68afa6ad-21d4-48fc-81da-5047c208bc71


# ╔═╡ 5c82c343-5d9d-4399-99f5-a9a5c9bfc9bb


# ╔═╡ bd865d0b-79d0-433c-879d-a819dd8f0b69


# ╔═╡ 32b0b550-a4d6-495e-b6ae-d00b2e527d5b


# ╔═╡ 2482a924-f3d7-405f-983e-8b4a2c378ae2


# ╔═╡ c4a98565-d14f-4548-b4a3-852047bbdd1d


# ╔═╡ 4439df44-d75a-473d-957d-298d549d5854


# ╔═╡ 7fb1c334-76f5-4419-9db9-e57ae7e26fa1


# ╔═╡ eb7e93de-291a-41ce-92f9-3e7e68369802


# ╔═╡ 7e3ab2f8-75f2-44c0-a3d1-12d0142005a0


# ╔═╡ d41d54fd-f268-42ac-8b87-48c753bbbd81


# ╔═╡ 0085976f-3f40-42bb-98f0-5e11264aa93a


# ╔═╡ 8e514b0f-1cee-4be6-81c4-32afe4bcdaba


# ╔═╡ f80d2269-cac2-4231-bacc-788f2a61008f


# ╔═╡ cee00905-abb6-4dcc-993e-4c396f4a1c6b


# ╔═╡ 9459403d-b5d9-41ae-9761-fd4d2f485d81


# ╔═╡ e416e726-ca0a-4b21-af27-24059803dac4


# ╔═╡ 5d68dccc-e0a1-4e83-8319-437473c98916


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
# ╠═524b736b-fbb7-4425-a55a-ee7a73bba776
# ╟─8443f122-db68-4807-a21d-5970335009c3
# ╠═69f6c0fa-b5e4-45e1-8e5c-5e5b5cfd5bec
# ╠═cdbde18f-c4d4-494e-b953-851fa07edc74
# ╟─e0ae1294-6efb-4b4a-9e20-292064c940cf
# ╟─b97425a6-6538-4fc5-82de-26475275a1f8
# ╠═a4f9ac12-a595-46d7-b8cc-870c27c4fa9e
# ╠═ea7a825f-a545-4527-975b-ef01dc92b5c5
# ╠═0cc39cdf-70d3-4a09-83e1-2296bc3a62ab
# ╠═d3daec1c-6972-4024-b3e0-8a1bdc384a1d
# ╠═0767c090-9177-42e5-bf61-3738b293bcfb
# ╠═04e8d2fa-fe02-4baf-ac0d-c6c28bfd8aee
# ╠═68afa6ad-21d4-48fc-81da-5047c208bc71
# ╠═5c82c343-5d9d-4399-99f5-a9a5c9bfc9bb
# ╠═bd865d0b-79d0-433c-879d-a819dd8f0b69
# ╠═32b0b550-a4d6-495e-b6ae-d00b2e527d5b
# ╠═2482a924-f3d7-405f-983e-8b4a2c378ae2
# ╠═c4a98565-d14f-4548-b4a3-852047bbdd1d
# ╠═4439df44-d75a-473d-957d-298d549d5854
# ╠═7fb1c334-76f5-4419-9db9-e57ae7e26fa1
# ╠═eb7e93de-291a-41ce-92f9-3e7e68369802
# ╠═7e3ab2f8-75f2-44c0-a3d1-12d0142005a0
# ╠═d41d54fd-f268-42ac-8b87-48c753bbbd81
# ╠═0085976f-3f40-42bb-98f0-5e11264aa93a
# ╠═8e514b0f-1cee-4be6-81c4-32afe4bcdaba
# ╠═f80d2269-cac2-4231-bacc-788f2a61008f
# ╠═cee00905-abb6-4dcc-993e-4c396f4a1c6b
# ╠═9459403d-b5d9-41ae-9761-fd4d2f485d81
# ╠═e416e726-ca0a-4b21-af27-24059803dac4
# ╠═5d68dccc-e0a1-4e83-8319-437473c98916
