"""
    complex_split_denominator(z, t, μ, nodes, weights, glnodes, glweights)

Computes the denominator of the exact solution `u` to the viscous Burgers' equation 
with initial condition `u(x, 0) = 1/(1+x^2)` and splits it into real and imaginary parts.    
"""
function complex_split_denominator(z, t, μ, nodes, weights, glnodes, glweights)
    D = viscous_solution_denominator(z..., t, μ, nodes, weights, glnodes, glweights)
    return [real(D), imag(D)]
end

"""
    saddle_point_approximation_μ_split(z, t, μ)

Computes the small-μ saddle point approximation and splits it into real 
and imaginary parts, assuming `z` is a pole of the exact solution.
"""
function saddle_point_approximation_μ_split(z, t, μ)
    D̃ = saddle_point_approximation_μ(complex(z...), t, μ)
    return [real(D̃), imag(D̃)]
end

"""
    Φ₀_split(z, μ)

Returns the real and imaginary parts, respectively, of `Φ₀(z, μ)`.
"""
function Φ₀_split(z, μ)
    U = Φ₀(complex(z...), μ)
    return [real(U), imag(U)]
end
