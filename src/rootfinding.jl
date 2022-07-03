"""
    newton_method(F, x::AbstractVector; τₐ=1e-8, τᵣ=1e-8, maxIters=500)
    newton_method(F, x; τₐ=1e-8, τᵣ=1e-8, maxIters=500)
    newton_method(F, F′, x; τₐ=1e-8, τᵣ=1e-8, maxIters=500)

Solve the equation `F(x) = 0` using Newton's method. The first method is for vectors, 
the second for scalars, and the third is for scalar problems with a provided derivative. 

# Arguments 
- `F`: Function to find the zero of.
- `F′`: Derivative of `F`.
- `x`: Initial guess for the solution to `F(x) = 0`.

# Keyword Arguments 
- `τₐ = 1e-8`: Absolute tolerance.
- `τᵣ = 1e-8`: Relative tolerance.
- `maxIters = 500`: Maximum allwoable number of Newton iterations.

# Outputs 
In the first method, the value is updated in-place into `x`. The value is returned for the 
latter two methods.
"""
function newton_method end
function newton_method(F, x::AbstractVector; τₐ=1e-8, τᵣ=1e-8, maxIters=500)
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
        x .= NaN
    end
    return nothing
end
function newton_method(F, x; τₐ=1e-8, τᵣ=1e-8, maxIters=500)
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
    return x
end
function newton_method(F, F′, x; τₐ=1e-8, τᵣ=1e-8, maxIters=500)
    Fx = F(x)
    Jx = F′(x)
    normF₀ = norm(Fx)
    res = missing
    for _ = 1:maxIters
        x -= Jx \ Fx
        Fx = F(x)
        Jx = F′(x)
        normF = norm(Fx)
        if normF ≤ τᵣ * normF₀ + τₐ
            res = normF
            break
        end
    end
    if ismissing(res)
        error("Failed to converge.")
    end
    return x
end

"""
    twopoint_linesearch(F, x; <keyword arguments>)

Solve `F(x) = 0` with initial guess `x` using a two-point linesearch.

# Arguments 
- `F`: Function to find the zero of.
- `x`: Initial guess for the solution to `F(x) = 0`.

# Keyword Arguments 
- `τₐ = 1e-8`: Absolute tolerance.
- `τᵣ = 1e-8`: Relative tolerance.
- `maxIters = 500`: Maximum allowable number of iterations.
- `σ₀ = 0.1`: Left end-point of interval used for safeguarding the λ value when determining the search direction.
- `σ₁ = 0.5`: Right end-pint of interval used for safeguarding the λ value when determining the search direction.
- `α = 1e-4`: The value of α to use in Armijo's rule.

# Outputs 
The found root is updated in-place into `x`.

"""
function twopoint_linesearch(F, x; τₐ=1e-8, τᵣ=1e-8, maxIters=500, σ₀=0.1, σ₁=0.5, α=1e-4)
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
        x .= NaN
    end
end