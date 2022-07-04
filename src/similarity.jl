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
    Φ₀(ξ, μ; type="parabolic")
    Φ₀(x, y, μ; type="parabolic")
    Φ₀(x::AbstractVector, y::AbstractVector, μ::AbstractVector; type="parabolic")
    Φ₀(x::AbstractVector{Float64}, y::AbstractVector{Float64}, μ::Float64; type="parabolic")

Computes the first-order small-time solution, given by `Φ₀ = [1/2√(2μ)]parabolicU(1/2 - i/4μ, iξ/√(2μ))/parabolicU(-1/2 - i/4μ, iξ/√(2μ))`. In 
the second method, `ξ = x + im y`, and similarly in the third method. The third method returns a matrix `A` such that `A[i, j, k]` is the value of 
`Φ₀` at `z = x[i] + im y[j]` and `μ[k]`.  The last method returns `A[i, j]`, which is the same as before except without `μ[k]` (as it is a scalar).
"""
function Φ₀ end
function Φ₀(ξ, μ; type="parabolic")
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
function Φ₀(x, y, μ; type="parabolic")
    return Φ₀(complex(x, y), μ; type=type)
end
function Φ₀(x::AbstractVector, y::AbstractVector, μ::AbstractVector; type="parabolic")
    Φ₀_vals = Array{ComplexF64}(undef, x, y, length(μ))
    for (k, m) in enumerate(μ)
        for (j, Y) in enumerate(y)
            for (i, X) in enumerate(x)
                Φ₀_vals[i, j, k] = Φ₀(X, Y, m; type=type)
            end
        end
    end
    return Φ₀_vals
end
function Φ₀(x::AbstractVector{Float64}, y::AbstractVector{Float64}, μ::Float64; type="parabolic")
    Φ₀_vals = Array{ComplexF64}(undef, x, y)
    for (j, Y) in enumerate(y)
        for (i, X) in enumerate(x)
            ξ = complex(X, Y)
            Φ₀_vals[i, j] = Φ₀(ξ, μ; type=type)
        end
    end
    return Φ₀_vals
end