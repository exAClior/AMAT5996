using Yao, Symbolics, LinearAlgebra, Plots


# Problem 1
## Define variables
@variables E₁ E₂ aᵣ aᵢ bᵣ bᵢ λ

a = aᵣ + im * aᵢ
b = bᵣ + im * bᵢ

## Define Hamiltonian
origin_ham = Matrix([E₁ 0 0; 0 E₁ 0; 0 0 E₂])
perturbation = Matrix([0 0 b; 0 0 a; conj(b) conj(a) 0])
hamiltonian = origin_ham + perturbation

### Find eigenvalues exactly
chara_eqn = det(hamiltonian - λ .* Matrix(I, 3, 3))

# Examine the form of charasteristic equation, we can see that E₁ is a root of this equation
# We then divie $λ - E₁$ from the equation to get a lower order polynomial which
# could be solved by quadratic formula

substitute(chara_eqn, Dict(λ => E₁))
chara_eqn = simplify(chara_eqn / (λ - E₁))
chara_eqn *= -1
const_term = substitute(chara_eqn, Dict(λ => 0))
linear_term = substitute(simplify((chara_eqn - const_term) / λ), Dict(λ => 0))
quadratic_term = simplify((chara_eqn - const_term - linear_term * λ) / λ^2)

b2m4ac = real(simplify(linear_term^2 - 4 * quadratic_term * const_term))

root2 = simplify((-linear_term + sqrt(b2m4ac)) / (2 * quadratic_term))
println(root2)
root3 = simplify((-linear_term - sqrt(b2m4ac)) / (2 * quadratic_term))
println(root3)

# Now let's use first order perturbation
# we know the Hamiltonian before perturbation has eigenvalues E₁, E₁, E₂
# and eigenvectors [1, 0, 0], [0, 1, 0], [0, 0, 1] respectively
# using (11) in paper: A Tutorial on Matrix Perturbation Theory we have
V0 = Matrix(1.0 * I, 3, 3)
Λ₁ = Diagonal(V0 * perturbation * V0)
# nothing is here

# Now let's use second order perturbation
Πdag = Matrix([0 0 1/(E₁-E₂); 0 0 1/(E₁-E₂); 1/(E₂-E₁) 1/(E₂-E₁) 0])
V1 = -V0 * (Πdag .* (V0 * perturbation * V0))
Λ₂ = Diagonal(V0 * perturbation * V1)

simplify(Λ₂[1, 1])
simplify(Λ₂[2, 2])
simplify(Λ₂[3, 3])

lamb1np = simplify(E₁ + Λ₂[1, 1])
lamb2np = simplify(E₁ + Λ₂[2, 2])
lamb3np = simplify(E₂ + Λ₂[3, 3])

# they don't agree


# Problem 2
# Hydrogen-like atom: is an atom in which all but one of the electrons have been removed.
#
#
#
#
#
