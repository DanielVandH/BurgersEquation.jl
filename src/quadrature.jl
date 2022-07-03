"""
    gauss_legendre_finite(f, a, b, nodes, weights):         Int[f(s), a, b]
    gauss_legendre_right_infinite(f, a, nodes, weights):    Int[f(s), a, ∞]
    gauss_legendre_left_infinite(f, a, nodes, weights):     Int[f(s), -∞, a]
    gauss_hermite(f, nodes, weights):                       Int[f(s), -∞, ∞]

Methods for Gauss-Legendre / Gauss-Hermite quadrature. In all cases, `(nodes, weights)` give the 
appropriate quadrature `nodes` and `weights` (see [`FastGaussQuadrature.gausslegendre`](@ref) 
or [`FastGaussQuadrature.gausshermite`](@ref)). 
"""
function gauss_legendre end
@doc (@doc gauss_legendre) function gauss_legendre_finite(f, a, b, nodes, weights)
    ℓ = (b-a)/2 
    m = (a+b)/2 
    f_new = s -> f(ℓ * s + m) * ℓ 
    return dot(weights, f_new.(nodes))
end
@doc (@doc gauss_legendre) function gauss_legendre_right_infinite(f, a, nodes, weights)
    g = s -> f(-cot((π - 2atan(a)) * (s - 1) / 4)) * csc(1 / 4 * (s - 1) * (π - 2atan(a)))^2
    return 1 / 4 * (π - 2atan(a)) * dot(weights, g.(nodes))
end
@doc (@doc gauss_legendre) function gauss_legendre_left_infinite(f, a, nodes, weights)
    g = s -> f(-cot((π + 2atan(a)) * (s + 1) / 4)) * csc(1 / 4 * (s + 1) * (π + 2atan(a)))^2
    return 1 / 4 * (π + 2atan(a)) * dot(weights, g.(nodes))
end
@doc (@doc gauss_legendre) function gauss_hermite(f, nodes, weights)
    return dot(weights, f.(nodes))
end