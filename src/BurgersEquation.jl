module BurgersEquation

using LinearAlgebra
using Makie
using FastGaussQuadrature
using HypergeometricFunctions
using ForwardDiff
using SpecialFunctions
using ApproxFun
using Colors
using ComplexPortraits

###################################################################################################
## Rootfinding functions 
###################################################################################################
"""
    Newton_Method(F::Function, x::Vector{Float64}; τₐ = 1e-8, τᵣ = 1e-8, maxIters = 100)

Uses Newton's method with automatic differentiation to solve the equation F(x) = 0.

# Arguments 
- `F`: The function to find the root of.
- `x`: The initial guess. This will get updated in-place with the solution.

# Keyword Arguments
- `τₐ = 1e-8`: The absolute tolerance.
- `τᵣ = 1e-8`: The relative tolerance.
- `maxIters = 100`: The maximum allowable number of Newton iterations.

# Outputs 
- `x`: This value is updated in-place from the input `x` with the approximate solution to `F(x) = 0`.
"""
function Newton_Method(F::Function, x::Vector{Float64}; τₐ=1e-8, τᵣ=1e-8, maxIters=100)
    result = DiffResults.JacobianResult(x)
    ForwardDiff.jacobian!(result, F, x)
    Fx = DiffResults.value(result)
    Jx = DiffResults.jacobian(result)
    normF₀ = norm(Fx)
    res = missing
    for _ = 1:maxIters
        x .-= Jx \ Fx
        ForwardDiff.jacobian!(result, F, x)
        Fx .= DiffResults.value(result)
        Jx .= DiffResults.jacobian(result)
        normF = norm(Fx)
        if normF ≤ τᵣ * normF₀ + τₐ
            res = normF
            break
        end
    end
    if ismissing(res)
        error("Failed to converge.")
    end
    return nothing
end

"""
    Newton_Method(F::Function, x::Float64, τₐ = 1e-8, τᵣ = 1e-8, maxIters = 100)

Method for calling [`Newton_Method`](@ref) when `x` is a number rather than a vector.
"""
function Newton_Method(F::Function, x::Float64; τₐ=1e-8, τᵣ=1e-8, maxIters=500)
    Fx = F(x)
    Jx = ForwardDiff.derivative(F, x)
    normF₀ = norm(Fx)
    res = missing
    for _ = 1:maxIters
        x -= Jx \ Fx
        Fx = F(x)
        Jx = ForwardDiff.derivative(F, x)
        normF = norm(Fx)
        if normF ≤ τᵣ * normF₀ + τₐ
            res = normF
            break
        end
    end
    if ismissing(res)
        error("Failed to converge.")
    end
    return nothing
end

"""
    Twopoint_LineSearch(F::Function, x::Vector{Float64}; <keyword arguments>)

Uses a two-point linesearch to solve `F(x) = 0` with initial guess `x`.

# Arguments 
- `F`: Function to find root of.
- `x`: Initial guess for the solution to `F(x) = 0`.

# Keyword Arguments 
- `τₐ = 1e-8`: Absolute tolerance.
- `τᵣ = 1e-8`: Relative tolernace.
- `maxIters = 500`: Maximum allowable number of iterations.
- `σ₀ = 0.1`: Left end-point of interval used for safeguarding the λ value when determining the search direction.
- `σ₁ = 0.5`: Right end-pint of interval used for safeguarding the λ value when determining the search direction.
- `α = 1e-4`: The value of α to use in Armijo's rule.

# Outputs 
The found root is over-written into `x`.
"""
function Twopoint_LineSearch(F::Function, x::Vector{Float64}; τₐ=1e-8, τᵣ=1e-8, maxIters=500,
    σ₀=0.1, σ₁=0.5, α=1e-4)
    result = DiffResults.JacobianResult(x)
    ForwardDiff.jacobian!(result, F, x)
    Fx = DiffResults.value(result)
    Jx = DiffResults.jacobian(result)
    normF = norm(Fx)
    res = missing
    τ = τᵣ * normF + τₐ
    for _ in 1:maxIters
        normx = norm(x)
        if normx == 0
            ε = sqrt(eps(Float64))
        else
            ε = sqrt(eps(Float64)) * normx
        end
        δx = -Jx \ Fx
        λ = 1.0
        x⁺ = x + λ * δx
        g₀ = normF^2
        normFx⁺ = norm(F(x⁺))
        normFx⁺² = normFx⁺^2
        while normFx⁺ ≥ (1.0 - α * λ) * normF
            λ₁ = copy(λ)
            num = g₀ * λ₁^2
            den = normFx⁺² + (2λ₁ - 1.0) * g₀
            λ = num / den
            if λ < σ₀ * λ₁
                λ = σ₀ * λ₁
            elseif λ > σ₁ * λ₁
                λ = σ₁ * λ₁
            end
            x⁺ = x + λ * δx
            normFx⁺ = norm(F(x⁺))
            normFx⁺² = normFx⁺^2
        end
        x .= x⁺
        ForwardDiff.jacobian!(result, F, x)
        Fx .= DiffResults.value(result)
        Jx .= DiffResults.jacobian(result)
        normF = norm(Fx)
        if normF ≤ τ
            res = normF
            break
        end
    end
    if ismissing(res)
        error("Failed to converge.")
    end
    return nothing
end

###################################################################################################
## Inviscid solution functions
###################################################################################################
"""
    Inviscid_Solution(f::Function, x::AbstractVector, t::Float64; Δt = t/100)

Returns the solution to `∂ₜu  + u∂ₓu = 0` with initial condition `u(x, 0) = f(x)`, using the implicit solution 
`u(x, t) = f(x - ut)`. The solution at time `t` is found by using Newton's method, solving at times 0, Δt, 2Δt, …, t.

# Arguments 
- `f`: Defines the initial condition `u(x, 0) = f(x)`.
- `x`: The value of `x` on the real line to solve at.
- `t`: The time to give the solution at.

# Keyword Arguments
- `Δt = t / 100`: The spacing to use for breaking down the times between `0` and `t`.

# Output 
- `u`: The solution `u` such that `u[i]` is the solution at `x[i]`.
"""
function Inviscid_Solution(f::Function, x::AbstractVector, t::Float64; Δt=t / 100)
    u = f.(x)  # Initial guess
    if t > 0
        t_vals = Δt:Δt:t
        for τ in t_vals
            F = u -> u .- f.(x .- u * τ) # Function to solve 
            Newton_Method(F, u)
        end
    end
    return u
end

"""
    Inviscid_Solution(f::Function, x::AbstractVector, t::AbstractVector; Δt = maximum(t) / 100)

Method for calling `Inviscid_Solution` with a vector of times. 
"""
function Inviscid_Solution(f::Function, x::AbstractVector, t::AbstractVector; Δt=maximum(t) / 100)
    solns = zeros(length(x), length(t))
    for (j, T) in enumerate(t)
        solns[:, j] .= Inviscid_Solution(f, x, T; Δt=Δt)
    end
    return solns
end

###################################################################################################
## Functions for quadrature
###################################################################################################
"""
    Gauss_Legendre_Finite(f, a, b, nodes, weights)

Approximates the integral `∫ₐᵇf(s)ds` using Gauss-Legendre quadrature, mapping 
`[a, b]` into `[-1, 1]` by taking `s = (b-a)ξ/2 + (a+b)/2`.

# Arguments 
- `f`: The integrand.
- `a`: The left endpoint.
- `b`: The right endpoint.
- `nodes`: Gauss-Legendre nodes from `FastGaussQuadrature.gausslegendre`.
- `weights`: Gauss-Legendre weights from `FastGaussQuadrature.gausslegendre`.

# Outputs 
- `I`: An approximation to `∫ₐᵇf(s)ds`.
"""
function Gauss_Legendre_Finite(f, a, b, nodes, weights)
    ℓ = (b - a) / 2
    m = (a + b) / 2
    f_new = s -> f(ℓ * s + m) * ℓ
    return dot(weights, f_new.(nodes))
end

"""
    Gauss_Legendre_Right_Infinite(f, a, nodes, weights)

Computes the integral `∫f(s)ds` over the interval `[a, ∞)` using Gauss-Legendre quadrature, mapping 
`[a, ∞)` into `[-1, 1]` by taking `s = - cot[(π - 2arctan(a))(ξ - 1)/4]`.

# Arguments 
- `f`: The integrand.
- `a`: The left endpoint.
- `nodes`: Gauss-Legendre nodes from `FastGaussQuadrature.gausslegendre`.
- `weights`: Gauss-Legendre weights from `FastGaussQuadrature.gausslegendre`.

# Outputs 
- `I`: An approximation to `∫f(s)ds` over the interval `[a, ∞)`.
"""
function Gauss_Legendre_Right_Infinite(f, a, nodes, weights)
    g = s -> f(-cot((π - 2atan(a)) * (s - 1) / 4)) * csc(1 / 4 * (s - 1) * (π - 2atan(a)))^2
    return 1 / 4 * (π - 2atan(a)) * dot(weights, g.(nodes))
end

"""
    Gauss_Legendre_Left_Infinite(f, a, nodes, weights)

Computes the integral ∫f(s)ds over the interval `(-∞, a]` using Gauss-Legendre quadrature, mapping 
`(-∞, a]` into `[-1, 1]` by taking `s = -cot[(π + 2arctan(a))(ξ + 1)/4]`.

# Arguments 
- `f`: The integrand.
- `a`: The right endpoint.
- `nodes`: Gauss-Legendre nodes from `FastGaussQuadrature.gausslegendre`.
- `weights`: Gauss-Legendre weights from `FastGaussQuadrature.gausslegendre`.

# Outputs 
- `I`: An approximation to `∫f(s)ds` over the interval `(-∞, a]`.
"""
function Gauss_Legendre_Left_Infinite(f, a, nodes, weights)
    g = s -> f(-cot((π + 2atan(a)) * (s + 1) / 4)) * csc(1 / 4 * (s + 1) * (π + 2atan(a)))^2
    return 1 / 4 * (π + 2atan(a)) * dot(weights, g.(nodes))
end

"""
    Gauss_Hermite(f, nodes, weights)

Computes the integral `∫f(s)ds` over `ℝ` using Gauss-Hermite quadrature.

# Arguments 
- `f`: The integrand.
- `nodes`: Gauss-Hermite nodes from `FastGaussQuadrature.gausshermite`.
- `weights`: Gauss-Hermite weights from `FastGaussQuadrature.gausshermite`.

# Outputs 
- `I`: An approximation to `∫f(s)ds` over `ℝ`.
"""
function Gauss_Hermite(f, nodes, weights)
    return dot(weights, f.(nodes))
end

###################################################################################################
## Functions for computing the viscous solution
###################################################################################################
"""
    Viscous_Solution(x::Float64, y::Float64, t::Float64, μ::Float64, nodes, weights, glnodes, glweights)

Evaluates the viscous solution to Burgers' equation with initial condition `f(x) = 1/(1+x^2)` in the complex plane.

# Arguments 
- `x::Float64`: The real part of `(x, y)` to find the solution at.
- `y::Float64`: The imaginary part of `(x, y)` to find the solution at.
- `t::Float64`: The time to give the solution at.
- `μ::Float64`: The viscosity.
- `nodes`: The Gauss-Hermite nodes from `FastGaussQuadrature.gausshermite`.
-` weights`: The Gauss-Hermite weights from `FastGaussQuadrature.gausshermite`.
- `glnodes`: The Gauss-Legendre nodes from `FastGaussQuadrature.gausslegendre`.
- `glweights`: The Gauss-Legendre nodes from `FastGaussQuadrature.gausslegendre`.

# Outputs 
- `I`: An approximation to the viscous solution at `(x, y)` at time `t` with viscosity `μ`.
"""
function Viscous_Solution(x::Float64, y::Float64, t::Float64, μ::Float64, nodes, weights, glnodes, glweights)
    # Evaluate numerator  
    if abs(y) < 1 # Do we need to worry about the branch cut at all?
        f = s -> s * exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s))
        I₁ = Gauss_Hermite(f, nodes, weights)
    elseif y > 1 # Do we need to deform around the branch cut in the upper-half plane?
        # Γ₁: To the left of the branch point
        f₁₅ = s -> s * exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₁¹ = Gauss_Legendre_Left_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅: To the right of the branch point
        I₁⁵ = Gauss_Legendre_Right_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₂, Γ₄: The vertical walls going around the branch point 
        prefactor = 0.5 / (im * μ * t) * sinh(π / (4μ)) * exp(-0.25x^2 / (μ * t) - 0.5im * x * (y - 1) / (μ * t))
        f₂₅ = τ -> (-x - im * (y - 1 - τ)) * exp(0.25 * (y - 1 - τ)^2 / (μ * t)) * exp(-0.25im / μ * log((τ + 2) / τ)) * exp(0.5im * x * τ / (μ * t))
        I₁²⁴ = prefactor * Gauss_Legendre_Finite(f₂₅, 0.0, y - 1, glnodes, glweights)
        I₁ = I₁¹ + I₁⁵ + I₁²⁴
    elseif y < -1 # Do we need to deform around the branch cut in the lower-half plane?
        # Γ₁: To the left of the branch point 
        f₁₅ = s -> s * exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₁¹ = Gauss_Legendre_Left_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅: To the right of the branch point 
        I₁⁵ = Gauss_Legendre_Right_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₂, Γ₄: The vertical walls going around the branch point
        prefactor = 0.5im / (μ * t) * sinh(π / (4μ)) * exp(-0.25x^2 / (μ * t) - 0.5im * x * (y + 1) / (μ * t))
        f₂₅ = τ -> (-x - im * (y + 1 - τ)) * exp(0.25 * (y + 1 - τ)^2 / (μ * t)) * exp(-0.25im / μ * log(τ / (τ - 2))) * exp(0.5im * x * τ / (μ * t))
        I₁²⁴ = prefactor * Gauss_Legendre_Finite(f₂₅, y + 1, 0.0, glnodes, glweights)
        I₁ = I₁¹ + I₁⁵ + I₁²⁴
    elseif abs(y) == 1.0 # Special case where we are on the line that directly runs into the branch point
        # Γ₁
        f₁₅ = s -> s * exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₁¹ = Gauss_Legendre_Left_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅
        I₁⁵ = Gauss_Legendre_Right_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        I₁ = I₁¹ + I₁⁵
    end
    # Evaluate denominator. See above for an understanding of each if-else condition.
    if abs(y) < 1
        f = s -> exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s))
        I₂ = Gauss_Hermite(f, nodes, weights)
    elseif y > 1
        # Γ₁
        f₁₅ = s -> exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₂¹ = Gauss_Legendre_Left_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅
        I₂⁵ = Gauss_Legendre_Right_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₂, Γ₄
        prefactor = sinh(π / (4μ)) / (im * sqrt(μ * t)) * exp(-0.25x^2 / (μ * t) - 0.5im * x * (y - 1) / (μ * t))
        f₂₅ = τ -> exp(0.25 * (y - 1 - τ)^2 / (μ * t)) * exp(-0.25im / μ * log((τ + 2) / τ)) * exp(0.5im * x * τ / (μ * t))
        I₂²⁴ = prefactor * Gauss_Legendre_Finite(f₂₅, 0.0, y - 1, glnodes, glweights)
        I₂ = I₂¹ + I₂⁵ + I₂²⁴
    elseif y < -1
        # Γ₁
        f₁₅ = s -> exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₂¹ = Gauss_Legendre_Left_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅
        I₂⁵ = Gauss_Legendre_Right_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₂, Γ₄
        prefactor = im / sqrt(μ * t) * sinh(π / (4μ)) * exp(-0.25x^2 / (μ * t) - 0.5im * x * (y + 1) / (μ * t))
        f₂₅ = τ -> exp(0.25 * (y + 1 - τ)^2 / (μ * t)) * exp(-0.25im / μ * log(τ / (τ - 2))) * exp(0.5im * x * τ / (μ * t))
        I₂²⁴ = prefactor * Gauss_Legendre_Finite(f₂₅, y + 1, 0.0, glnodes, glweights)
        I₂ = I₂¹ + I₂⁵ + I₂²⁴
    elseif abs(y) == 1.0
        # Γ₁
        f₁₅ = s -> exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₂¹ = Gauss_Legendre_Left_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅
        I₂⁵ = Gauss_Legendre_Right_Infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        I₂ = I₂¹ + I₂⁵
    end
    return -2sqrt(μ / t) * I₁ / I₂
end

"""
    Viscous_Solution(x::Float64, t::Float64, μ::Float64, nodes, weights)

Evaluates the viscous solution to Burgers' equation with initial condition `f(x) = 1/(1+x^2)` on the real line.
Gauss-Hermite quadrature is used.

# Arguments 
- `x::Float64`: The point on the real line to give the solution at.
- `t::Float64`: The time to give the solution at.
- `μ::Float64`: The viscosity.
- `nodes`: The Gauss-Hermite nodes from `FastGaussQuadrature.gausshermite`.
-` weights`: The Gauss-Hermite weights from `FastGaussQuadrature.gausshermite`.

# Outputs 
- `I`: An approximation to the viscous solution at `x` at time `t` with viscosity `μ`.
"""
function Viscous_Solution(x::Float64, t::Float64, μ::Float64, nodes, weights)
    # Evaluate numerator  
    f = s -> s * exp(-0.5 / μ * atan(x + 2sqrt(μ * t) * s))
    I₁ = Gauss_Hermite(f, nodes, weights)
    # Evaluate denominator
    f = s -> exp(-0.5 / μ * atan(x + 2sqrt(μ * t) * s))
    I₂ = Gauss_Hermite(f, nodes, weights)
    return -2sqrt(μ / t) * I₁ / I₂
end

"""
    Viscous_Solution(x::AbstractVector, y::AbstractVector, t::AbstractVector, μ::AbstractVector)

Method for calling [`Viscous_Solution`](@ref) for a vector of points, times, and viscosities, in the complex plane.
The output is such that `u[i, j, k, ℓ]` is the solution at `x = x[i]`, `y = y[j]`, `t = t[k]`, and `μ = μ[ℓ]`.
"""
function Viscous_Solution(x::AbstractVector, y::AbstractVector, t::AbstractVector, μ::AbstractVector)
    num_nodes = 250
    nodes, weights = gausshermite(num_nodes)
    glnodes, glweights = gausslegendre(num_nodes)
    u = Array{ComplexF64}(zeros(length(y), length(x), length(t), length(μ)))
    for (ℓ, m) in enumerate(μ)
        for (k, T) in enumerate(t)
            for (i, X) in enumerate(x)
                for (j, Y) in enumerate(y)
                    if T > 0
                        u[i, j, k, ℓ] = Viscous_Solution(X, Y, T, m, nodes, weights, glnodes, glweights)
                    else
                        u[i, j, k, ℓ] = 1 / (1 + (X + im * Y)^2)
                    end
                end
            end
        end
    end
    return u
end

"""
    Viscous_Solution(x::AbstractVector, t::AbstractVector, μ::AbstractVector)

Method for calling [`Viscous_Solution`](@ref) for a vector of points, times, and viscosities, all on the real line. 
The output is such that `u[i, k, ℓ]` is the solution at `x = x[i]`, `t = t[k]`, and `μ = μ[ℓ]`.
"""
function Viscous_Solution(x::AbstractVector, t::AbstractVector, μ::AbstractVector)
    num_nodes = 250
    nodes, weights = gausshermite(num_nodes)
    u = zeros(length(x), length(t), length(μ))
    for (ℓ, m) in enumerate(μ)
        for (k, T) in enumerate(t)
            for (i, X) in enumerate(x)
                if T > 0
                    u[i, k, ℓ] = Viscous_Solution(X, T, m, nodes, weights)
                else
                    u[i, k, ℓ] = 1 / (1 + X^2)
                end
            end
        end
    end
    return u
end

###################################################################################################
## Functions for plotting phase portraits and analytical landscapes
###################################################################################################
"""
    portrait!(fig, x::AbstractVector, y::AbstractVector, Z::Matrix{ComplexF64}, i, j; 
    xlabel=L"x", ylabel=L"y", pltkwargs...)
   
Plots the phase portrait for the data in `Z`, defined over a structured grid.

# Arguments 
- `fig`: A Makie figure object.
- `x::AbstractVector`: The real part of the points that `Z` is computed at.
- `y::AbstractVector`: The imaginary part of the points that `Z` is computed at.
- `Z::Complex{Float64}`: The matrix of values to plot, defined such that `Z[i, j]` is at `(x[i], y[j])`.
- `i`: The row number to create the axis in `fig`.
- `j`: The column number to create the axis in `fig`.

# Keyword Arguments 
- `xlabel=L"x"`: Label for the x-axis.
- `ylabel=L"y"`: Label for the y-axis.
- `nist = true`: Whether to use the NIST color scheme. If false, the HSV color scheme is used.
- `pltkwargs...`: Extra keyword arguments for `Axis`.

# Output 
There is no output, but `fig` is updated with the figure in `fig[i, j]`.
"""
function portrait!(fig, x::AbstractVector, y::AbstractVector, Z::Matrix{ComplexF64}, i, j;
    xlabel=L"x", ylabel=L"y", nist=true, pltkwargs...)
    # Phase
    θ = angle.(-Z)
    # Obtain the NIST colormap http://dlmf.nist.gov/help/vrml/aboutcolor) or other
    nc = max(length(x), length(y))
    cm = nist ? ComplexPortraits.nist_colors(nc) : ComplexPortraits.hsv_colors(nc)
    # Make the plot 
    ax = Axis(fig[i, j], xlabel=xlabel, ylabel=ylabel; pltkwargs...)
    heatmap!(ax, x, y, θ, colorrange=(-π, π), colormap=cm)
    return nothing
end

"""
    landscape!(fig, x::AbstractVector, y::AbstractVector, Z::Matrix{ComplexF64}, i, j; 
    xlabel=L"x", ylabel=L"y", zlabel=L"|u(z, t)|", pltkwargs...)
   
Plots the analytical landscape for the data in `Z`, defined over a structured grid.

# Arguments 
- `fig`: A Makie figure object.
- `x::AbstractVector`: The real part of the points that `Z` is computed at.
- `y::AbstractVector`: The imaginary part of the points that `Z` is computed at.
- `Z::Complex{Float64}`: The matrix of values to plot, defined such that `Z[i, j]` is at `(x[i], y[j])`.
- `i`: The row number to create the axis in `fig`.
- `j`: The column number to create the axis in `fig`.

# Keyword Arguments 
- `xlabel=L"x"`: Label for the x-axis.
- `ylabel=L"y"`: Label for the y-axis.
- `zlabel=L"|u(z, t)|"`: Label for the z-axis.
- `nist = true`: Whether to use the NIST color scheme. If false, the HSV color scheme is used.
- `pltkwargs...`: Extra keyword arguments for `Axis3`.

# Output 
There is no output, but `fig` is updated with the figure in `fig[i, j]`.
"""
function landscape!(fig, x::AbstractVector, y::AbstractVector, Z::Matrix{ComplexF64}, i, j;
    xlabel=L"x", ylabel=L"y", zlabel="|u(z, t)|", nist=true, zlims=(0, maximum(abs.(Z[isfinite.(Z)]))), pltkwargs...)
    # Modulus and phase 
    R = abs.(Z)
    R[R.>zlims[2]] .= zlims[2]
    θ = angle.(-Z)
    # Obtain the NIST colormap http://dlmf.nist.gov/help/vrml/aboutcolor)
    nc = max(length(x), length(y))
    cm = nist ? ComplexPortraits.nist_colors(nc) : ComplexPortraits.hsv_colors(nc)
    # Make the plot 
    ax = Axis3(fig[i, j], xlabel=xlabel, ylabel=ylabel, zlabel=zlabel; pltkwargs...)
    surface!(ax, x, y, R, color=θ, colorrange=(-π, π), colormap=cm, shading=false)
    zlims!(ax, zlims[1], zlims[2])
    return nothing
end

###################################################################################################
## Functions for plotting the similarity solution
###################################################################################################
"""
    parabolicU(a, z)

Computes the parabolic cylinder function `U(a, z)`. We use Eq. 27 of https://mathworld.wolfram.com/ParabolicCylinderFunction.html.
"""
function parabolicU(a, z)
    Y₁ = 1 / sqrt(π) * exp(loggamma(0.25 - 0.5a)) / 2^(0.5a + 0.25) * exp(-0.25z^2.0) * HypergeometricFunctions.M(0.5a + 0.25, 0.5, 0.5z^2.0)
    Y₂ = 1 / sqrt(π) * exp(loggamma(0.75 - 0.5a)) / 2^(0.5a - 0.25) * z * exp(-0.25z^2.0) * HypergeometricFunctions.M(0.5a + 0.75, 1.5, 0.5z^2.0)
    U = cos(π * (0.25 + 0.5a)) * Y₁ - sin(π * (0.25 + 0.5a)) * Y₂
    return U
end

"""
    Φ₀(ξ, μ; type = "parabolic")

Computes the first-order small-time solution, given by `Φ₀ = [1/2√(2μ)]U(1/2 - i/4μ, iξ/√(2μ))/U(-1/2 - i/4μ, iξ/√(2μ))`.
"""
function Φ₀(ξ::ComplexF64, μ::Float64; type="parabolic")
    if type == "kummer"
        a₂ = -im / (8μ)
        a₁ = 1 + a₂
        b₂ = 0.5
        b₁ = 1 + b₂
        c = -ξ^2 / (4μ)
        U₁ = HypergeometricFunctions.U(a₁, b₁, c)
        U₂ = HypergeometricFunctions.U(a₂, b₂, c)
        return (ξ * im / (8μ)) * U₁ / U₂
    elseif type == "parabolic"
        a₁ = 0.5 - im / (4μ)
        a₂ = a₁ - 1.0
        b = im * ξ / sqrt(2μ)
        U₁ = parabolicU(a₁, b)
        U₂ = parabolicU(a₂, b)
        scale = 2sqrt(2μ)
        return U₁ / (scale * U₂)
    end
end

"""
    Φ₀(ℜξ::AbstractVector{Float64}, ℑξ::AbstractVector{Float64}, μ::AbstractVector{Float64}; type = "parabolic")

Computes the first-order small-time solution `Φ₀ = [1/2√(2μ)]U(1/2 - i/4μ, iξ/√(2μ))/U(-1/2 - i/4μ, iξ/√(2μ))`
for an array of values. The result is a matrix `A` such that `A[i, j, k]` is the solution at `ξ = ℜξ[i] + imℑξ[j]` and `μ = μ[k]`.
"""
function Φ₀(ℜξ::AbstractVector{Float64}, ℑξ::AbstractVector{Float64}, μ::AbstractVector{Float64}; type="parabolic")
    Φ₀_vals = Array{ComplexF64}(undef, length(ℜξ), length(ℑξ), length(μ))
    for (k, m) in enumerate(μ)
        for (j, y) in enumerate(ℑξ)
            for (i, x) in enumerate(ℜξ)
                ξ = complex(x, y)
                Φ₀_vals[i, j, k] = Φ₀(ξ, m; type=type)
            end
        end
    end
    return Φ₀_vals
end

end
