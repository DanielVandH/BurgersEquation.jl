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
function maximum_slope(t::Float64, μ::Float64, nodes, weights, x₀; ic=1, alg=LBFGS())
    if μ > 0
        ∇u = x -> abs(ForwardDiff.derivative(x -> viscous_solution(x[1], t, μ, nodes, weights; ic=ic), x[1]))
    else
        if ic == 1
            f′ = y -> -2y / (1 + y^2)^2
            ∇u = x -> begin
                x = x[1]
                u = inviscid_solution(y -> 1 / (1 + y^2), x, t)
                res = abs(f′(x - u * t) / (1 + t * f′(x - u * t)))
                if res > 1e6
                    return NaN
                else
                    return res
                end
            end
        elseif ic == 2
            f′ = y -> -4y / (1 + y^2)^3
            ∇u = x -> begin
                x = x[1]
                u = inviscid_solution(y -> 1 / (1 + y^2)^2, x, t)
                res = abs(f′(x - u * t) / (1 + t * f′(x - u * t)))
                if res > 1e6
                    return NaN
                else
                    return res
                end
            end
        elseif ic == 3
            f′ = y -> -y / (1 + y^2)^(3 / 2)
            ∇u = x -> begin
                x = x[1]
                u = inviscid_solution(y -> 1 / (1 + y^2)^(1 / 2), x, t)
                res = abs(f′(x - u * t) / (1 + t * f′(x - u * t)))
                if res > 1e6
                    return NaN
                else
                    return res
                end
            end
        end
    end
    sol = Optim.optimize(x -> -∇u(x), [x₀], alg; autodiff=:forward)
    max_slope = -sol.minimum
    return sol.minimizer[1], max_slope
end
function maximum_slope(t::AbstractVector, μ::AbstractVector; ic=1, alg=LBFGS(), init=0.0, num_nodes=250)
    nodes, weights = gausshermite(num_nodes)
    slopes = zeros(length(t), length(μ))
    positions = zeros(length(t), length(μ))
    @inbounds for (j, m) in enumerate(μ)
        for (i, τ) in enumerate(t)
            @show i / length(t), j / length(μ)
            if i == 1 || positions[i-1, j] ≠ NaN
                try
                    positions[i, j], slopes[i, j] = maximum_slope(τ, m, nodes, weights, i == 1 ? init : positions[i-1, j]; ic=ic, alg=alg)
                catch
                end
                if i > 1
                    if abs(positions[i, j] - positions[i-1, j]) > 0.01
                        positions[i:end, j] .= NaN
                        slopes[i:end, j] .= NaN
                    end
                end
            end
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