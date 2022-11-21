"""
    maximum_slope(t::Float64, μ::Float64, nodes, weights, x₀; ic = 1)
    maximum_slope(t::AbstractVector, μ::AbstractVector)

For a given `μ` and `t`, finds the maximum slope for the viscous solution on the real line. The latter method does this 
for all pairs of `t` and `μ`,

# Arguments 
- `T`: Time to solve at.
- `μ`: Viscosity.
- `nodes`: Gauss-Hermite nodes from `FastGaussQuadrature.gausshermite`.
- `weights`: Gauss-Hermite weights from `FastGaussQuadrature.gausshermite`.
- `x₀`: An initial guess for the position of the maximum slope.

Use `ic` to control the initial condition; see [`viscous`](@ref).

# Outputs 
- The maximum absolute slope.
"""
function maximum_slope end
function maximum_slope(t::Float64, μ::Float64, nodes, weights, x₀; ic=1)
    if μ > 0
        ∇u = x -> abs(ForwardDiff.derivative(x -> viscous_solution(x[1], t, μ, nodes, weights; ic=ic), x[1]))
    else
        if ic == 1
            ∇u = x -> abs(ForwardDiff.derivative(x -> inviscid_solution(y -> 1 / (1 + y^2), x[1], t), x[1]))
        elseif ic == 2
            ∇u = x -> abs(ForwardDiff.derivative(x -> inviscid_solution(y -> 1 / (1 + y^2)^2, x[1], t), x[1]))
        elseif ic == 3
            ∇u = x -> abs(ForwardDiff.derivative(x -> inviscid_solution(y -> 1 / (1 + y^2)^(1 / 2), x[1], t), x[1]))
        end
    end
    sol = Optim.optimize(x -> -∇u(x), [x₀], LBFGS(); autodiff=:forward)
    max_slope = -sol.minimum
    return sol.minimizer[1], max_slope
end
function maximum_slope(t::AbstractVector, μ::AbstractVector; ic=1)
    num_nodes = 250
    nodes, weights = gausshermite(num_nodes)
    slopes = zeros(length(t), length(μ))
    positions = zeros(length(t), length(μ))
    @inbounds for (j, m) in enumerate(μ)
        for (i, τ) in enumerate(t)
            positions[i, j], slopes[i, j] = maximum_slope(τ, m, nodes, weights, i == 1 ? 0.0 : positions[i-1, j]; ic=ic)
        end
    end
    return slopes
end

"""
    findmaxima(slopes::Matrix{Float64})

Returns the (1) indices and (2) maxima for each column of `slopes`.
If no maxima can be found, the index returned is `-1` and the peak is `NaN`.
"""
function Peaks.findmaxima(slopes::Matrix{Float64})
    indices = Vector{Int64}(undef, size(slopes, 2))
    vals = Vector{Float64}(undef, size(slopes, 2))
    for (i, x) in enumerate(eachcol(slopes))
        try
            idx, pk = findmaxima(x)
            indices[i], vals[i] = idx[1], pk[1]
        catch
            indices[i], vals[i] = -1, NaN
        end
    end
    return indices, vals
end