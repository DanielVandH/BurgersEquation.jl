"""
    viscous_solution_large_time_Ψ(η, μ)
    viscous_solution_large_time_Ψ(x::AbstractVector, y::AbstractVector, μ::AbstractVector)
    viscous_solution_large_time_Ψ(x::AbstractVector, μ::AbstractVector)
    viscous_solution_large_time_Ψ(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Number}
    viscous_solution_large_time(x, t::Float64, μ::Float64k, t₀ = 0.0)
    viscous_solution_large_time(x::AbstractVector{T}, t::AbstractVector{Float64}, μ::AbstractVector{Float64}, t₀ = 0.0) where {T<:Number}
    viscous_solution_large_time(x::AbstractVector{T}, y::AbstractVector{T}, t::AbstractVector{Float64}, μ::AbstractVector{Float64}, t₀ = 0.0) where {T<:Number}
    viscous_solution_large_time(x::AbstractVector{T}, t, μ::AbstractVector{Float64}, t₀ = 0.0) where {T<:Number}
    
Function for computing the large-time limit of the exact solution to the viscous problem.
"""
function viscous_solution_large_time end
@doc (@doc viscous_solution_large_time) function viscous_solution_large_time_Ψ(η, μ)
    2 / sqrt(π) * exp(-η^2 / 4) / (coth(π / (4μ)) - erf(η / 2))
end
@doc (@doc viscous_solution_large_time) function viscous_solution_large_time_Ψ(x::AbstractVector, y::AbstractVector, μ::AbstractVector)
    sol = zeros(ComplexF64, length(x), length(y), length(μ))
    for (k, m) in enumerate(μ)
        for (j, yy) in enumerate(y)
            for (i, xx) in enumerate(x)
                sol[i, j, k] = viscous_solution_large_time_Ψ(complex(xx, yy), m)
            end
        end
    end
    sol
end
@doc (@doc viscous_solution_large_time) function viscous_solution_large_time_Ψ(x::AbstractVector, μ::AbstractVector)
    sol = zeros(Float64, length(x), length(μ))
    for (j, m) in enumerate(μ)
        for (i, xx) in enumerate(x)
            sol[i, j] = viscous_solution_large_time_Ψ(xx, m)
        end
    end
    sol
end
function viscous_solution_large_time(x, t::Float64, μ::Float64, t₀ = 0.0)
    η = x/sqrt(μ*(t-t₀))
    Ψ = viscous_solution_large_time_Ψ(η, μ)
    return sqrt(μ/(t-t₀)) * Ψ
end
function viscous_solution_large_time(x::AbstractVector{T}, t::AbstractVector{Float64}, μ::AbstractVector{Float64}, t₀ = 0.0) where {T<:Number}
    vals = zeros(T, length(x), length(t), length(μ))
    for (k, μμ) in enumerate(μ)
        for (j, tt) in enumerate(t)
            for (i, xx) in enumerate(x)
                vals[i, j, k] = viscous_solution_large_time(xx, tt, μμ, t₀)
            end
        end
    end
    return vals
end
function viscous_solution_large_time(x::AbstractVector{T}, t, μ::AbstractVector{Float64}, t₀ = 0.0) where {T<:Number}
    viscous_solution_large_time(x, [t], μ, t₀)[:, 1, :]
end
function viscous_solution_large_time(x::AbstractVector{T}, y::AbstractVector{T}, t::AbstractVector{Float64}, μ::AbstractVector{Float64}, t₀ = 0.0) where {T<:Number}
    vals = zeros(ComplexF64, length(x), length(y), length(t), length(μ))
    for (ℓ, μμ) in enumerate(μ)
        for (k, tt) in enumerate(t)
            for (j, yy) in enumerate(y)
                for (i, xx) in enumerate(x)
                    z = xx + im * yy
                    vals[i, j, k, ℓ] = viscous_solution_large_time(z, tt, μμ, t₀)
                end
            end
        end
    end
    return vals
end