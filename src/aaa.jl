"""
viscous_solution_aaa(u, x, t, X, Y)
viscous_solution_aaa(N, L, T, t, μ::Float64, X, Y)
viscous_solution_aaa(N, L, T, t, μ::Vector{Float64}, X, Y)

Given some numerical data `(u, x, t)` for the solution to viscous Burgers' equation, 
computes its analytic continuation in the complex plane specified by the points `(X, Y)` 
using the AAA algorithm. The latter two methods start by obtaining this data 
using a finite difference solution.
"""
function viscous_solution_aaa end
function viscous_solution_aaa(u, x, t, X, Y)
    U = zeros(ComplexF64, length(X), length(Y), length(t))
    for k in 1:length(t)
        r = aaa(x, u[:, k])
        for (i, ξ) in enumerate(X)
            for (j, ζ) in enumerate(Y)
                U[i, j, k] = r(ξ + im * ζ)
            end
        end
    end
    return U
end
function viscous_solution_aaa(N, L, T, t, μ::Float64, X, Y)
    u, x, t = viscous_solution_finite_diff(N, L, T, μ, t)
    return viscous_solution_aaa(u, x, t, X, Y)
end
function viscous_solution_aaa(N, L, T, t, μ::Vector{Float64}, X, Y)
    U = zeros(ComplexF64, length(X), length(Y), length(t), length(μ))
    for (i, m) in enumerate(μ)
        U[:, :, :, i] .= viscous_solution_aaa(N, L, T, t, m, X, Y)
    end
    return U
end

