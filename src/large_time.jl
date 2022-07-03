"""
viscous_solution_large_time(x, μ::Float64, t::Float64)
viscous_solution_large_time(x::Vector{T}, μ::Vector{Float64}, t::Vector{Float64}) where {T<:Number}
viscous_solution_large_time(x::Vector{T}, y::Vector{T}, μ::Vector{Float64}, t::Vector{Float64}) where {T<:Number}
viscous_solution_large_time_Ψ(x::Vector{T}, y::Vector{T}, μ::Vector{Float64}) where {T<:Number}

Function for computing the large-time limit of the exact solution to the viscous problem. The last method computes 
`Ψ` rather than `u`.
"""
function viscous_solution_large_time end 
function viscous_solution_large_time(x, t::Float64, μ::Float64)
    return 2sqrt(μ / (π * t)) * exp(-x^2 / (4 * μ * t)) / (coth(π / (4μ)) - erf(x / (2sqrt(μ * t))))
end
function viscous_solution_large_time(x::Vector{T}, t::Vector{Float64}, μ::Vector{Float64}) where {T<:Number}
    vals = zeros(T, length(x), length(t), length(μ))
    for (k, μμ) in enumerate(μ)
        for (j, tt) in enumerate(t)
            for (i, xx) in enumerate(x)
                vals[i, j, k] = viscous_solution_large_time(xx, tt, μμ)
            end
        end
    end
    return vals
end
function viscous_solution_large_time(x::Vector{T}, y::Vector{T}, t::Vector{Float64}, μ::Vector{Float64}) where {T<:Number}
    vals = zeros(ComplexF64, length(x), length(y), length(t), length(μ))
    for (ℓ, μμ) in enumerate(μ)
        for (k, tt) in enumerate(t)
            for (j, yy) in enumerate(y)
                for (i, xx) in enumerate(x)
                    z = xx + im * yy
                    vals[i, j, k, ℓ] = viscous_solution_large_time(z, tt, μμ)
                end
            end
        end
    end
    return vals
end
@doc (@doc viscous_solution_large_time) function viscous_solution_large_time_Ψ(x::Vector{T}, y::Vector{T}, μ::Vector{Float64}) where {T<:Number}
    vals = viscous_solution_large_time(x, y, μ, [1.0])
    vals = vals[:, :, 1, :]
    return vals
end