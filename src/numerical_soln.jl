"""
    heateqs!(dφ, φ, p, t)

Computes the system of ordinary differential equations for solving 
the heat equation with homogeneous Neumann boundary conditions.
"""
function heateqs!(dφ, φ, p, t)
    Δx, N, μ = p
    dφ[1] = 2μ * (φ[2] - φ[1]) / Δx^2
    for i in 2:N-1
        dφ[i] = μ * (φ[i+1] - 2φ[i] + φ[i-1]) / Δx^2
    end
    dφ[N] = -2μ * (φ[N] - φ[N-1]) / Δx^2
    return nothing
end

"""
    heat_equation_finite_diff(N, L, T, μ)

Solves the heat equation with homogeneous Neumann boundary conditions 
with `N` grid points on the domain `[-L, L]` with diffusivity `μ` up to 
a time `T`.
"""
function heat_equation_finite_diff(N, L, T, μ)
    Δx = 2L / (N - 1)
    x = [-L + (i - 1) * Δx for i in 1:N]
    φ₀ = @. exp(-atan(x) / (2μ))
    p = (Δx, N, μ)
    prob = ODEProblem(heateqs!, φ₀, (0, T), p)
    sol = solve(prob)
    return sol
end

"""
    viscous_solution_finite_diff(N, L, T, μ, t)
    viscous_solution_finite_diff(N, L, T, μ::Vector{Float64}, t)

Computes the numerical solution to Burgers' equation with 
viscosity `μ` on the domain `[-L, L]` with homogeneous Dirichlet 
boundary conditions. The solution is given at the times specified 
in `t`. If `μ` is a vector, solutions are returned for each `μ[i]`.
"""
function viscous_solution_finite_diff end
function viscous_solution_finite_diff(N, L, T, μ::Float64, t)
    φ_sol = heat_equation_finite_diff(N, L, T, μ)
    Δx = 2L / (N - 1)
    M = length(t)
    x = [-L + (i - 1) * Δx for i in 1:N]
    u = zeros(N, M)
    φ_temp = zeros(N)
    for (j, τ) in enumerate(t)
        if τ > 0
            φ_temp .= φ_sol(τ)
            for i in 2:N-1
                u[i, j] = μ * (φ_temp[i-1] - φ_temp[i+1]) / (φ_temp[i] * Δx)
            end
        else
            for i in 2:N-1
                u[i, j] = 1 / (1 + x[i]^2)
            end
        end
    end
    return u, x, t
end
function viscous_solution_finite_diff(N, L, T, μ::Vector{Float64}, t)
    Δx = 2L / (N - 1)
    M = length(t)
    x = [-L + (i - 1) * Δx for i in 1:N]
    u = zeros(N, M, length(μ))
    for (i, m) in enumerate(μ)
        uu, _, _ = viscous_solution_finite_diff(N, L, T, m, t)
        u[:, :, i] .= uu
    end
    return u, x, t
end