"""
    viscous_solution_numerator(x, y, t, μ, nodes, weights, glnodes, glweights)
    viscous_solution_denominator(x, y, t, μ, nodes, weights, glnodes, glweights)

Evaluates the numerator/denominator for the viscous solution to Burgers' equation with initial condition `f(x) = 1/(1+x^2)` in the complex plane.

# Arguments 
- `x`: The real part of `(x, y)` to find the solution at.
- `y`: The imaginary part of `(x, y)` to find the solution at.
- `T`: The time to give the solution at.
- `μ`: The viscosity.
- `nodes`: The Gauss-Hermite nodes from [`FastGaussQuadrature.gausshermite`](@ref).
- `weights`: The Gauss-Hermite weights from [`FastGaussQuadrature.gausshermite`](@ref).
- `glnodes`: The Gauss-Legendre nodes from [`FastGaussQuadrature.gausslegendre`](@ref).
- `glweights`: The Gauss-Legendre weights from [`FastGaussQuadrature.gausslegendre`](@ref).

# Outputs 
- `I`: An approximation to the numerator/denominator for the viscous solution at `(x, y)` at time `t` with viscosity `μ`.
"""
function viscous_solution_numerator(x, y, t, μ, nodes, weights, glnodes, glweights)
    # Evaluate numerator  
    if abs(y) < 1 # Do we need to worry about the branch cut at all?
        f = s -> s * exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s))
        I₁ = gauss_hermite(f, nodes, weights)
    elseif y > 1 # Do we need to deform around the branch cut in the upper-half plane?
        # Γ₁: To the left of the branch point
        f₁₅ = s -> s * exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₁¹ = gauss_legendre_left_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅: To the right of the branch point
        I₁⁵ = gauss_legendre_right_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₂, Γ₄: The vertical walls going around the branch point 
        prefactor = 0.5 / (im * μ * t) * sinh(π / (4μ)) * exp(-0.25x^2 / (μ * t) - 0.5im * x * (y - 1) / (μ * t))
        f₂₅ = τ -> (-x - im * (y - 1 - τ)) * exp(0.25 * (y - 1 - τ)^2 / (μ * t)) * exp(-0.25im / μ * log((τ + 2) / τ)) * exp(0.5im * x * τ / (μ * t))
        I₁²⁴ = prefactor * gauss_legendre_finite(f₂₅, 0.0, y - 1, glnodes, glweights)
        I₁ = I₁¹ + I₁⁵ + I₁²⁴
    elseif y < -1 # Do we need to deform around the branch cut in the lower-half plane?
        # Γ₁: To the left of the branch point 
        f₁₅ = s -> s * exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₁¹ = gauss_legendre_left_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅: To the right of the branch point 
        I₁⁵ = gauss_legendre_right_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₂, Γ₄: The vertical walls going around the branch point
        prefactor = 0.5im / (μ * t) * sinh(π / (4μ)) * exp(-0.25x^2 / (μ * t) - 0.5im * x * (y + 1) / (μ * t))
        f₂₅ = τ -> (-x - im * (y + 1 - τ)) * exp(0.25 * (y + 1 - τ)^2 / (μ * t)) * exp(-0.25im / μ * log(τ / (τ - 2))) * exp(0.5im * x * τ / (μ * t))
        I₁²⁴ = prefactor * gauss_legendre_finite(f₂₅, y + 1, 0.0, glnodes, glweights)
        I₁ = I₁¹ + I₁⁵ + I₁²⁴
    elseif abs(y) == 1.0 # Special case where we are on the line that directly runs into the branch point
        # Γ₁
        f₁₅ = s -> s * exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₁¹ = gauss_legendre_left_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅
        I₁⁵ = gauss_legendre_right_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        I₁ = I₁¹ + I₁⁵
    end
    return I₁
end
@doc (@doc viscous_solution_numerator) function viscous_solution_denominator(x, y, t, μ, nodes, weights, glnodes, glweights)
    if abs(y) < 1
        f = s -> exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s))
        I₂ = gauss_hermite(f, nodes, weights)
        return I₂
    elseif y > 1
        # Γ₁
        f₁₅ = s -> exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₂¹ = gauss_legendre_left_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅
        I₂⁵ = gauss_legendre_right_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₂, Γ₄
        prefactor = sinh(π / (4μ)) / (im * sqrt(μ * t)) * exp(-0.25x^2 / (μ * t) - 0.5im * x * (y - 1) / (μ * t))
        f₂₅ = τ -> exp(0.25 * (y - 1 - τ)^2 / (μ * t)) * exp(-0.25im / μ * log((τ + 2) / τ)) * exp(0.5im * x * τ / (μ * t))
        I₂²⁴ = prefactor * gauss_legendre_finite(f₂₅, 0.0, y - 1, glnodes, glweights)
        I₂ = I₂¹ + I₂⁵ + I₂²⁴
        return I₂
    elseif y < -1
        # Γ₁
        f₁₅ = s -> exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₂¹ = gauss_legendre_left_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅
        I₂⁵ = gauss_legendre_right_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₂, Γ₄
        prefactor = im / sqrt(μ * t) * sinh(π / (4μ)) * exp(-0.25x^2 / (μ * t) - 0.5im * x * (y + 1) / (μ * t))
        f₂₅ = τ -> exp(0.25 * (y + 1 - τ)^2 / (μ * t)) * exp(-0.25im / μ * log(τ / (τ - 2))) * exp(0.5im * x * τ / (μ * t))
        I₂²⁴ = prefactor * gauss_legendre_finite(f₂₅, y + 1, 0.0, glnodes, glweights)
        I₂ = I₂¹ + I₂⁵ + I₂²⁴
        return I₂
    elseif abs(y) == 1.0
        # Γ₁
        f₁₅ = s -> exp(-0.5 / μ * atan(x + im * y + 2sqrt(μ * t) * s)) * exp(-s^2)
        I₂¹ = gauss_legendre_left_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        # Γ₅
        I₂⁵ = gauss_legendre_right_infinite(f₁₅, -0.5x / sqrt(μ * t), glnodes, glweights)
        I₂ = I₂¹ + I₂⁵
        return I₂
    end
    return NaN + im * NaN
end

"""
    viscous_solution(x::T, y::T, t::T, μ::T, nodes, weights, glnodes, glweights) where {T}
    viscous_solution(x, t::T, μ::T, nodes, weights) where {T}
    viscous_solution(x::AbstractVector, y::AbstractVector, t::AbstractVector, μ::AbstractVector; num_nodes = 250)
    viscous_solution(x::AbstractVector, y::AbstractVector, t::AbstractVector, μ::Float64)
    viscous_solution(x::AbstractVector, y::AbstractVector, t::Float64, μ::Float64)
    viscous_solution(x::AbstractVector, t::AbstractVector, μ::AbstractVector)
    viscous_solution(x::AbstractVector, t::AbstractVector, μ::Float64)

Computes the viscous solution to Burgers' equation with initial condition `u(x, 0) = 1/(1+x^2)`. For the first two methods,
the result is simply a scalar. For the latter methods, the solution is such that `u[i, j, k, ℓ]` is the solution at 
`z = x[i] + im y[j]`, `t[k]`, `μ[ℓ]`, with dimensions dropped appropriately depending on the inputs given. Additionally,

- `nodes`: The Gauss-Hermite nodes from [`FastGaussQuadrature.gausshermite`](@ref).
- `weights`: The Gauss-Hermite weights from [`FastGaussQuadrature.gausshermite`](@ref).
- `glnodes`: The Gauss-Legendre nodes from [`FastGaussQuadrature.gausslegendre`](@ref).
- `glweights`: The Gauss-Legendre weights from [`FastGaussQuadrature.gausslegendre`](@ref).
"""
function viscous_solution end 
# Yes, there's a better way to do what we do below with a better usage of multiple dispatch.
function viscous_solution(x::T, y::T, t::T, μ::T, nodes, weights, glnodes, glweights) where {T}
    I₁ = viscous_solution_numerator(x, y, t, μ, nodes, weights, glnodes, glweights)
    I₂ = viscous_solution_denominator(x, y, t, μ, nodes, weights, glnodes, glweights)
    return -2sqrt(μ / t) * I₁ / I₂
end
function viscous_solution(x, t::T, μ::T, nodes, weights) where {T}
    if t > 0
        # Evaluate numerator  
        f = s -> s * exp(-0.5 / μ * atan(x + 2sqrt(μ * t) * s))
        I₁ = gauss_hermite(f, nodes, weights)
        # Evaluate denominator
        f = s -> exp(-0.5 / μ * atan(x + 2sqrt(μ * t) * s))
        I₂ = gauss_hermite(f, nodes, weights)
        return -2sqrt(μ / t) * I₁ / I₂
    else
        return 1 / (1 + x^2)
    end
end
function viscous_solution(x::AbstractVector, y::AbstractVector, t::AbstractVector, μ::AbstractVector; num_nodes = 250)
    nodes, weights = gausshermite(num_nodes)
    glnodes, glweights = gausslegendre(num_nodes)
    u = Array{ComplexF64}(zeros(length(x), length(y), length(t), length(μ)))
    for (ℓ, m) in enumerate(μ)
        for (k, T) in enumerate(t)
            for (i, X) in enumerate(x)
                for (j, Y) in enumerate(y)
                    if T > 0
                        u[i, j, k, ℓ] = viscous_solution(X, Y, T, m, nodes, weights, glnodes, glweights)
                    else
                        u[i, j, k, ℓ] = 1 / (1 + (X + im * Y)^2)
                    end
                end
            end
        end
    end
    return u
end
function viscous_solution(x::AbstractVector, y::AbstractVector, t::AbstractVector, μ::Float64)
    num_nodes = 250
    nodes, weights = gausshermite(num_nodes)
    glnodes, glweights = gausslegendre(num_nodes)
    u = Array{ComplexF64}(zeros(length(x), length(y), length(t)))
    for (k, T) in enumerate(t)
        for (i, X) in enumerate(x)
            for (j, Y) in enumerate(y)
                if T > 0
                    u[i, j, k] = viscous_solution(X, Y, T, μ, nodes, weights, glnodes, glweights)
                else
                    u[i, j, k] = 1 / (1 + (X + im * Y)^2)
                end
            end
        end
    end
    return u
end
function viscous_solution(x::AbstractVector, y::AbstractVector, t::Float64, μ::Float64)
    num_nodes = 250
    nodes, weights = gausshermite(num_nodes)
    glnodes, glweights = gausslegendre(num_nodes)
    u = Array{ComplexF64}(zeros(length(x), length(y)))
    for (i, X) in enumerate(x)
        for (j, Y) in enumerate(y)
            if t > 0
                u[i, j] = viscous_solution(X, Y, t, μ, nodes, weights, glnodes, glweights)
            else
                u[i, j] = 1 / (1 + (X + im * Y)^2)
            end
        end
    end
    return u
end
function viscous_solution(x::AbstractVector, t::AbstractVector, μ::AbstractVector)
    num_nodes = 250
    nodes, weights = gausshermite(num_nodes)
    u = zeros(length(x), length(t), length(μ))
    for (ℓ, m) in enumerate(μ)
        for (k, T) in enumerate(t)
            for (i, X) in enumerate(x)
                if T > 0
                    u[i, k, ℓ] = viscous_solution(X, T, m, nodes, weights)
                else
                    u[i, k, ℓ] = 1 / (1 + X^2)
                end
            end
        end
    end
    return u
end
function viscous_solution(x::AbstractVector, t::AbstractVector, μ::Float64)
    num_nodes = 250
    nodes, weights = gausshermite(num_nodes)
    u = zeros(length(x), length(t))
    for (j, T) in enumerate(t)
        for (i, X) in enumerate(x)
            if T > 0
                u[i, j] = Viscous_Solution(X, T, μ, nodes, weights)
            else
                u[i, j] = 1 / (1 + X^2)
            end
        end
    end
    return u
end


