"""
    cubic_saddle_points(x, t)

Computes the roots to the cubic polynomial `s³ - s^2x + s + t - x`.
"""
function cubic_saddle_points(x, t)
    A = 2^(2 / 3) * (18.0x - 27.0t + 3.0 * sqrt(3.0) * sqrt(4.0 + 27.0t^2 - 36.0t * x + 8.0x^2 - 4.0t * x^3 + 4.0x^4) + 2.0x^3)^(1 / 3) / 6
    B = (3 - x^2) / (9A)
    s₁ = x / 3 + A - B
    s₂ = x / 3 + (-1.0 - im * sqrt(3.0)) * A / 2.0 - (-1.0 + im * sqrt(3.0)) * B / 2
    s₃ = x / 3 + (-1.0 + im * sqrt(3.0)) * A / 2.0 - (-1.0 - im * sqrt(3.0)) * B / 2
    return s₁, s₂, s₃
end

"""
    cubic_discriminant(x, t)

Computes the discriminant of the cubic polynomial `s³ - s^2x + s + t - x`.
"""
function cubic_discriminant(x, t)
    return -4.0 - 27.0t^2 + 36.0t * x - 8.0x^2 + 4.0t * x^3 - 4.0x^4
end

""" 
    cubic_caustic(x, t)

Plots the values of `cubic_discriminant(x, t)` that correspond to `Δ = 0`, for a given `x`.
"""
function cubic_caustic(x)
    t₁ = 2x / 3 + 2x^3 / 27 - 2sqrt((x^2 - 3)^3) / 27
    t₂ = 2x / 3 + 2x^3 / 27 + 2sqrt((x^2 - 3)^3) / 27
    return t₁, t₂
end

"""
    saddle_point_approximation_μ(z, t, μ)

Computes the small-μ saddle point approximation, `D(μ) ∼ sqrt(2μπ)[D₁ + D₂]` at time `t`,
for the denominator of the solution to the viscous Burgers' equation with initial 
condition `u(x, 0) = 1/(1+x^2)`, where `Dⱼ = 1/sqrt(|h′′(sⱼ)|)exp(iθⱼ + h(sⱼ)/μ)`, 
assuming that `z` is a pole of the exact solution of `u`.
"""
function saddle_point_approximation_μ(z, t, μ)
    # Define the phase function 
    h(s) = -0.5atan(s) - (z - s)^2 / (4t)
    h′(s) = 0.5z / t - 0.5s / t - 0.5 / (s^2 + 1.0)
    h′′(s) = s / (s^2 + 1.0)^2 - 0.5 / t
    # Select the saddle points 
    s₁, s₂, s₃ = cubic_saddle_points(z, t)
    Diff₁₂ = abs(real(h(s₁)) - real(h(s₂)))
    Diff₁₃ = abs(real(h(s₁)) - real(h(s₃)))
    Diff₂₃ = abs(real(h(s₂)) - real(h(s₃)))
    Diffs = [Diff₁₂, Diff₁₃, Diff₂₃]
    idx = findmin(Diffs)[2]
    if idx == 1
        s₁, s₂ = s₁, s₂
    elseif idx == 2
        s₁, s₂ = s₁, s₃
    elseif idx == 3
        s₁, s₂ = s₂, s₃
    end
    # Compute θ
    α₁ = angle(h′′(s₁))
    α₂ = angle(h′′(s₂))
    θ₁¹ = -0.5α₁ + 0.5π
    θ₁² = -0.5α₁ + 1.5π
    if cos(θ₁¹) ≥ 0
        θ₁ = θ₁¹
    else
        θ₁ = θ₁²
    end
    θ₂¹ = -0.5α₂ + 0.5π
    θ₂² = -0.5α₂ + 1.5π
    if cos(θ₂¹) ≥ 0
        θ₂ = θ₂¹
    else
        θ₂ = θ₂²
    end
    # Compute the approximations
    D₁ = abs(h′′(s₁))^(-0.5) * exp(im * θ₁ + h(s₁) / μ)
    D₂ = abs(h′′(s₂))^(-0.5) * exp(im * θ₂ + h(s₂) / μ)
    D̃ = sqrt(2μ * π) * (D₁ + D₂)
    return D̃
end