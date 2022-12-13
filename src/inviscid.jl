"""
    inviscid_solution(f, x::AbstractVector, t::Number; Δt=t / 100)
    inviscid_solution(f, x::AbstractVector, t::AbstractVector; Δt=maximum(t) / 100)
    inviscid_solution(f, x::Number, t::Number; Δt=maximum(t) / 100)
    inviscid_solution(f, f′, x, y, t::Number; Δt=t / 100)

Computes the solution to the inviscid Burgers equation `uₜ + uuₓ = 0` with initial 
condition `u(x, 0) = f(x)` using the implicit solution `u(x, t) = f(x - ut)`. The 
solution at time `t` is found by using Newton's method, solving at times `0:Δt:t`.

# Arguments 
- `f`: The initial condition `u(x, 0) = f(x)`.
- `f′`: If computing the solution in the complex plane, you need to provide the derivative `f′` of `f`.
- `x`: The values of `x` on the real line to compute the solution at.
- `y`: If this argument is computed, then the solutions are given for `complex(x, y)`.
- `t`: The time(s) to compute the solution at.

# Keyword Arguments 
- `Δt = maximum(t) / 100`: The spacing to use for breaking down the times between `0` and `t`.

# Output 
- `u`: The solution `u` such that `u[i]` is the solution at `x[i]` if `t` is a scalar, and `u[i, j]` is the solution at `x[i]`, `t[j]` if `t` is a vector. 

If instead `y` is provided, then the solutions are for `u[i, j]` for the solution at `complex(x[i], y[j])` when `t` is a scalar, 
and `u[i, j, k]` would be at `complex(x[i], y[j])` when `t = t[k]`.

Yes, there are better ways to do everything here via multiple dispatch.
"""
function inviscid_solution end
function inviscid_solution(f, x::AbstractVector, t::Number; Δt=t / 100)
    u = f.(x)
    if t > 0
        t_vals = Δt:Δt:t
        for τ in t_vals
            F = u -> u .- f.(x .- u * τ)
            newton_method(F, u; maxIters=100)
        end
    end
    return u
end
function inviscid_solution(f, x::AbstractVector, t::AbstractVector; Δt=maximum(t) / 100)
    solns = zeros(eltype(x), length(x), length(t))
    for (j, T) in enumerate(t)
        solns[:, j] .= inviscid_solution(f, x, T; Δt)
    end
    return solns
end
function inviscid_solution(f, x::Number, t::Number; Δt=maximum(t) / 100)
    solns = inviscid_solution(f, [x], t; Δt)
    return solns[1]
end
function inviscid_solution(f, f′, x, y, t::Number; Δt=t / 100)
    z = complex(x, y)
    u = f(z)
    if t > 0
        t_vals = Δt:Δt:t
        for τ in t_vals
            F = u -> u - f(z - u * τ)
            F′ = u -> 1 + τ * f′(z - u * τ)
            u = newton_method(F, F′, u; maxIters=100)
        end
    end
    return u
end
function inviscid_solution(f, f′, x::Number, t::Number; Δt=t / 100)
    u = f(x)
    if t > 0
        t_vals = Δt:Δt:t
        for τ in t_vals
            F = u -> u - f(x - u * τ)
            F′ = u -> 1 + τ * f′(x - u * τ)
            u = newton_method(F, F′, u; maxIters=100)
        end
    end
    return u
end
function inviscid_solution(f, f′, x::AbstractVector, y::AbstractVector, t::AbstractVector; Δt=t / 100)
    u = zeros(ComplexF64, length(x), length(y), length(t))
    for (i, X) in enumerate(x)
        for (j, Y) in enumerate(y)
            for (k, T) in enumerate(t)
                u[i, j, k] = inviscid_solution(f, f′, X, Y, T; Δt=Δt[k])
            end
        end
    end
    return u
end
function inviscid_solution(f::Function, f′::Function, x::AbstractVector, t::AbstractVector)
    u = zeros(length(x), length(t))
    u[:, 1] .= f.(x)
    for j in 2:lastindex(t)
        τ = t[j]
        prev_u = @views u[:, j-1]
        for i in eachindex(x)
            F = u -> u - f(x[i] - u * τ)
            F′ = u -> 1.0 + τ * f′(x[i] - u .* τ)
            try
                u[i, j] = newton_method(F, F′, prev_u[i])
            catch
                u[i, j] = NaN
            end
        end
    end
    return u
end

"""
    inviscid_singularities(t)

Computes the locations of the branch points in the inviscid solution at time `t`. Returns the 
solution values and the locations, ordered by first, second, third, and fourth quadrant.
"""
function inviscid_singularities(t)
    a = 4t^2
    b = -4t^2
    c = d = 0.0
    e = 1.0
    p = (8 * a * c - 3 * b^2) / (8 * a^2)
    q = (b^3 - 4 * a * b * c + 8 * a^2 * d) / (8 * a^3)
    Δ₀ = c^2 - 3 * b * d + 12 * a * e
    Δ₁ = 2 * c^2 - 9 * b * c * d + 27 * b^2 * e + 27 * a * d^2 - 72 * a * c * e
    Δ = 256 * a^3 * e^3 - 192 * a^2 * c^2 * e^2 + 144 * a^2 * c * d^2 * e - 27 * a^2 * d^4 + 144 * a * b^2 * c * e^2 - 6 * a * b^2 * d^2 * e - 80 * a * b * c^2 * d * e + 18 * a * b * c * d^3 + 16 * a * c^4 * e - 4 * a * c^3 * d^2 - 27 * b^4 * e^2 + 18 * b^3 * c * d * e - 4 * b^3 * d^3 - 4 * b^2 * c^3 * e + b^2 * c^2 * d^2
    if Δ₀ ≠ 0.0 || (Δ == 0.0 && Δ₀ == 0.0)
        Q = ((Δ₁ + sqrt(Complex(Δ₁^2 - 4Δ₀^3))) / 2)^(1 / 3)
    else
        Q = cbrt(Δ₁)
    end
    S = 1 / 2 * sqrt(-2 / 3 * p + 1 / (3 * a) * (Q + Δ₀ / Q))
    if S == 0.0
        Q = Q * exp(2π * im / 3)
        S = 1 / 2 * sqrt(-2 / 3 * p + 1 / (3 * a) * (Q + Δ₀ / Q))
        if S == 0.0
            Q = Q * exp(2π * im / 3)
            S = 1 / 2 * sqrt(-2 / 3 * p + 1 / (3 * a) * (Q + Δ₀ / Q))
        end
    end
    ξ = (-1.0, -1.0, 1.0, 1.0)
    ζ = (1.0, -1.0, 1.0, -1.0)
    u = @. -b / (4 * a) + ξ * S + 1 / 2 * ζ * sqrt(-4 * S^2 - 2 * p - ξ * q / S)
    z = @. u * t + ξ * sqrt(1 / u - 1)
    u = [u for u in u]
    z = [z for z in z]
    i_idx = Vector{Int64}([])
    for i in 1:4
        if real(z[i]) ≥ 0.0 && imag(z[i]) ≥ 0.0
            push!(i_idx, 1)
        elseif real(z[i]) ≤ 0.0 && imag(z[i]) ≥ 0.0
            push!(i_idx, 2)
        elseif real(z[i]) ≤ 0.0 && imag(z[i]) ≤ 0.0
            push!(i_idx, 3)
        elseif real(z[i]) ≥ 0.0 && imag(z[i]) ≤ 0.0
            push!(i_idx, 4)
        end
    end
    new_idx = sortperm(i_idx)
    return u[new_idx], z[new_idx]
end

"""
    exact_inviscid(x, t)

Computes the inviscid solution with the initial condition `1/(1+x^2)`, returning the three 
roots to `u = f(x - ut)`.
"""
function exact_inviscid(x, t)
    a = t^2
    b = -2 * t * x
    c = 1 + x^2
    d = -1.0
    Δ₀ = b^2 - 3 * a * c
    Δ₁ = 2 * b^3 - 9 * a * b * c + 27 * a^2 * d
    C = cbrt((Δ₁ + sqrt(Δ₁^2 - 4Δ₀^3)) / 2)
    if C == 0
        C = cbrt((Δ₁ - sqrt(Δ₁^2 - 4Δ₀^3)) / 2)
    end
    if C ≠ 0
        ξ₁ = -1 / (3a) * (b + C + Δ₀ / C)
        ζ = (-1 + sqrt(3) * im) / 2
        ξ₂ = -1 / (3a) * (b + ζ * C + Δ₀ / (ζ * C))
        ξ₃ = -1 / (3a) * (b + ζ^2 * C + Δ₀ / (ζ^2 * C))
        ξr = ξ₁, ξ₂, ξ₃
    else
        ξr = -b / (3a), -b / (3a), -b / (3a)
    end
    ξr_reals = real.(ξr)
    i = argmin(ξr_reals)
    return real(ξr[i])
end

