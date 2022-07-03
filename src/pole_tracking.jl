"""
    locate_pole(z₀, t, μ, nodes, weights, glnodes, glweights; newton = true, kwargs...)
    locate_pole(z₀, t, μ)

Estimates the location of a pole in the denominator at time `t` and viscosity `μ` 
of the exact solution based on an initial guess `z₀`, using numerical quadrature for the computation.
If `newton = true`, Newton's method is used, otherwise a two-point linesearch is used. The additional
`kwargs` are for the nonlinear solvers. The latter method uses the saddle point approximation to the 
exact solution.
"""
function locate_pole end
function locate_pole(z₀, t, μ, nodes, weights, glnodes, glweights; newton=true, kwargs...)
    D = z -> complex_split_denominator(z, t, μ, nodes, weights, glnodes, glweights)
    z = [real(z₀), imag(z₀)]
    newton ? newton_method(D, z; kwargs...) : twopoint_linesearch(D, z; kwargs...)
    if all(isnan, z)
        return NaN + im * NaN
    end
    return complex(z...)
end
function locate_pole(z₀, t, μ)
    D = z -> saddle_point_approximation_μ_split(z, t, μ)
    z = [real(z₀), imag(z₀)]
    twopoint_lineSearch(D, z)
    return complex(z...)
end

"""
    Φ₀_pole(z, z₀, μ)

Returns an approximate solution `z` for the solution to `Φ₀(z, μ) = ∞`, using the initial guess 
`z = z₀`.
"""
function Φ₀_pole(z₀, μ; kwargs...)
    try
        U = z -> 1.0 ./ Φ₀_split(z, μ)
        z = [real(z₀), imag(z₀)]
        Newton_Method(U, z; kwargs...)
        return complex(z...)
    catch
        return NaN + im * NaN
    end
end

"""
    Ψ_pole(μ, z₀)

Search for poles to the large-limit function `Ψ` by solving `coth(π/4μ) = erf(η/2√μ)`, 
starting at `μ[end]` with initial guess `z₀`, and then backtracking in `μ`.
"""
function Ψ_pole(μ, z₀)
    roots = zeros(ComplexF64, length(μ))
    F = η -> coth(π / (4μ[end])) - erf(η / (2sqrt(μ[end])))
    F′ = η -> -exp(-η^2 / (4μ[end])) / sqrt(μ[end] * π)
    roots[end] = newton_method(F, F′, z₀)
    for j in (length(μ)-1):-1:1
        F = η -> coth(π / (4μ[j])) - erf(η / (2sqrt(μ[j])))
        F′ = η -> -exp(-η^2 / (4μ[j])) / sqrt(μ[j] * π)
        try
            roots[j] = newton_method(F, F′, roots[j+1])
        catch
            roots[j] = roots[j+1]
        end
    end
    return roots
end

"""
    tracking_poles_saddle(z, t, μ)
    tracking_poles_exact(z, t, μ)
    tracking_poles_saddle(z::Vector{ComplexF64}, t::Vector{Float64}, μ::Vector{Float64})
    tracking_poles_exact(z::Vector{ComplexF64}, t::Vector{Float64}, μ::Vector{Float64})
    tracking_poles_exact(μ::AbstractVector{Float64}, z₀::ComplexF64; t_max=10.0, t_min=1e-6)
    tracking_poles_Φ₀(μ, z₀; Δμ = 1e-5, μ_max = 1.0)
    tracking_poles_aaa(μ::Float64, z₀::ComplexF64; t_max=10.0, t_min=1e-6, L=30, xN=2000)

Tracks the closest pole to the real line in the complex plane for the exact solution `u` to 
Burgers' equation with initial condition `u(x, 0) = 1/(1+x^2)`, using the saddle point approximation or the exact
solution (evaluated with quadrature). An initial estimate `z` for the pole location at a time `t` must be provided. 
It is assumed that `μ` is small if the saddle point approximation is being used. 
    
For the methods which take in vectors, results are given for each pair `(z[j], t[j], μ[j])`. 

For the latter `exact` method, poles are tracked backwards in time, starting at `t_max` and ending at `t_min`.
The argument `z₀` is an initial guess for the pole location at `t_max` and `μ = μ[end]`. 

The method for `Φ₀` does so for the similarity solution `Φ₀`, starting from an initial guess of the pole location `z₀` at `μ` and returning 
the locations over `μ:Δμ:μ_max`.

The AAA method tracks poles backwards in time using the AAA algorithm, starting at `t_max` and ending at `t_min`.
The argument `z₀` is an initial guess for the pole location at `t_max` and `μ = μ[end]`.
"""
function tracking_poles end 
@doc (@doc tracking_poles) function tracking_poles_saddle(z, t, μ)
    # Setup
    N = 10000
    Δt = t / N
    t_vals = 0:Δt:t
    t_idx = N + 1
    pole_locs = Vector{ComplexF64}(zeros(N + 1))
    pole_locs[t_idx] = locate_pole(z, t, μ)
    # Track poles from t, t-Δt, t-2Δt, …, 0
    for j = t_idx:-1:2
        try
            pole = locate_pole(pole_locs[j], t_vals[j-1], μ)
            if isnan(pole)
                pole_locs[j-1] = pole_locs[j]
            else
                pole_locs[j-1] = pole
            end
        catch
            pole_locs[j-1] = pole_locs[j]#.= NaN + im * NaN
            #break
        end
    end
    return t_vals, pole_locs
end
@doc (@doc tracking_poles) function tracking_poles_exact(z, t, μ)
    # Setup
    N = 20000
    Δt = t / N
    t_vals = 0:Δt:t
    t_idx = N + 1
    pole_locs = Vector{ComplexF64}(zeros(N + 1))
    num_nodes = 250
    nodes, weights = gausshermite(num_nodes)
    glnodes, glweights = gausslegendre(num_nodes)
    pole_locs[t_idx] = locate_pole(z, t, μ, nodes, weights, glnodes, glweights; newton=false)
    # Track poles from t, t-Δt, t-2Δt, …, 0
    for j = t_idx:-1:2
        try
            pole = locate_pole(pole_locs[j], t_vals[j-1], μ, nodes, weights, glnodes, glweights; newton=false)
            if isnan(pole)
                pole_locs[j-1] = pole_locs[j]
            else
                pole_locs[j-1] = pole
            end
        catch
            pole_locs[j-1] = pole_locs[j]#.= NaN + im * NaN
            #break
        end
        if abs(pole_locs[j-1] - pole_locs[j]) > 0.1
            pole_locs[j-1] = pole_locs[j]
        end
    end
    return t_vals, pole_locs
end
@doc (@doc tracking_poles) function tracking_poles_saddle(z::Vector{ComplexF64}, t::Vector{Float64}, μ::Vector{Float64})
    t_vals = Vector{Vector{Float64}}(undef, length(z))
    pole_locs = Vector{Vector{ComplexF64}}(undef, length(z))
    for j in 1:length(z)
        t_vals[j], pole_locs[j] = tracking_poles_saddle(z[j], t[j], μ[j])
    end
    return t_vals, pole_locs
end
@doc (@doc tracking_poles) function tracking_poles_exact(z::Vector{ComplexF64}, t::Vector{Float64}, μ::Vector{Float64})
    t_vals = Vector{Vector{Float64}}(undef, length(z))
    pole_locs = Vector{Vector{ComplexF64}}(undef, length(z))
    for j in 1:length(z)
        t_vals[j], pole_locs[j] = tracking_poles_exact(z[j], t[j], μ[j])
    end
    return t_vals, pole_locs
end
@doc (@doc tracking_poles) function tracking_poles_exact(μ::AbstractVector{Float64}, z₀::ComplexF64; t_max=10.0, t_min=1e-6)
    # Setup
    N = 1001
    Δt = (t_max - t_min) / (N - 1)
    t_vals = t_min:Δt:t_max
    num_nodes = 500
    ghnodes, ghweights = gausshermite(num_nodes)
    glnodes, glweights = gausslegendre(num_nodes)
    pole_locs = Array{ComplexF64}(undef, length(μ), N) # pole_locs[i, j] is the estimated pole location at μ[i], t_vals[j]
    # Now refine the initial guess 
    pole_locs[end, N] = locate_pole(z₀, t_max, μ[end], ghnodes, ghweights, glnodes, glweights; newton=false)
    # Now track the poles, one μ at a time. Need to start with μ[end]
    for j in (N-1):-1:1
        z₀_new = pole_locs[end, j+1]
        try
            pole_locs[end, j] = locate_pole(z₀_new, t_vals[j], μ[end], ghnodes, ghweights, glnodes, glweights; newton=true, maxIters=500, τₐ=1e-2, τᵣ=1e-2)
        catch
            #pole_locs[end, j:-1:1] .= NaN + im * NaN
            pole_locs[end, j] = z₀_new
            #break
        end
        if abs(pole_locs[end, j] - z₀_new) > 0.1 # Sometimes we find a pole, but it's just so obviously wrong 
            pole_locs[end, j] = z₀_new
        end
    end
    # Now that we have these poles, we can move onto the next μ.
    for i in (length(μ)-1):-1:1
        pole_locs[i, N] = locate_pole(pole_locs[i+1, N], t_max, μ[i], ghnodes, ghweights, glnodes, glweights; newton=true, maxIters=500, τₐ=1e-2, τᵣ=1e-2)
        for j in (N-1):-1:1
            z₀_new = pole_locs[i, j+1]
            try
                pole_locs[i, j] = locate_pole(z₀_new, t_vals[j], μ[i], ghnodes, ghweights, glnodes, glweights; newton=true, maxIters=500, τₐ=1e-2, τᵣ=1e-2)
                if isnan(pole_locs[i, j])
                    #pole_locs[i, j:-1:1] .= NaN + im * NaN
                    #break
                    pole_locs[i, j] = z₀_new
                end
            catch
                #pole_locs[i, j:-1:1] .= NaN + im * NaN
                #break
                pole_locs[i, j] = z₀_new
            end
            if abs(pole_locs[i, j] - z₀_new) > 0.1 # Sometimes we find a pole, but it's just so obviously wrong 
                pole_locs[i, j] = z₀_new
            end
        end
    end
    return t_vals, pole_locs
end
@doc (@doc tracking_poles) function tracking_poles_Φ₀(μ, z₀; Δμ=1e-5, μ_max=1.0, kwargs...)
    μ_vals = μ:Δμ:μ_max
    num_μ = length(μ_vals)
    pole_locs = Vector{ComplexF64}(zeros(num_μ))
    pole_locs[1] = Φ₀_pole(z₀, μ; kwargs...)
    if isnan(pole_locs[1])
        pole_locs[1] = z₀
    end
    for j in 2:num_μ
        pole_locs[j] = Φ₀_pole(pole_locs[j-1], μ_vals[j]; kwargs...)
        if isnan(pole_locs[j])
            pole_locs[j] = pole_locs[j-1]
        end
    end
    return pole_locs, μ_vals
end
@doc (@doc tracking_poles) function tracking_poles_aaa(μ::Float64, z₀::ComplexF64; t_max=10.0, t_min=1e-6, L=30, xN=2000)
    ## Computes the solutions on the real line 
    M = 2500
    Δt = (t_max - t_min) / (M - 1)
    t_vals = t_min:Δt:t_max
    U, X, _ = Viscous_Solution_Finite_Diff(xN, L, t_max, μ, t_vals)
    ## Compare using chebpts (closest points to, at least)
    Xc = 5chebpts(250)
    idx = Vector{Int64}(zeros(250))
    for i in 1:250
        j = argmin(abs.(X .- Xc[i]))
        Xc[i] = X[j]
        idx[i] = j
    end
    unique!(idx)
    unique!(Xc)
    ## Preallocate pole locations 
    pole_locs = zeros(ComplexF64, M)

    ## Refine the initial guess 
    cleanup_fnc = (pol, res, r) -> begin
        cond₁ = @. abs(res) < 1e-13         # Residues too small
        cond₂ = @. abs(imag(res)) < 1e-12   # Unlikely to have real residues 
        cond₃ = @. abs(r(pol)) < 1e3        # Poles should have large function values 
        cond₄ = @. abs(imag(pol)) < 1e-12   # Not going to have real poles in this problem
        findall(cond₁ .|| cond₂ .|| cond₃ .|| cond₄)
    end
    postcleanup_fnc! = (pol, res, r) -> begin
        cond₁ = @. abs(res) < 1e-4     
        ii = findall(cond₁)
        deleteat!(pol, ii)
        return nothing
    end
    @views r = aaa(Xc, U[idx, end, end]; cleanup_fnc)
    pp, rr, _ = prz(r)
    postcleanup_fnc!(pp, rr, r)
    ii = argmin(@. abs(pp - z₀))
    #pole_locs[end, end] = pp[ii]
    pole_locs[end] = pp[ii]
    if abs(pole_locs[end] - z₀) > 0.1
        pole_locs[end] = z₀
    end

    ## Now, for given given μ, backtrack 
    for j in (M-1):-1:1
        try
            @views r = aaa(Xc, U[idx, j, end]; cleanup_fnc)
            pp, _, _ = prz(r)
            #ii = argmin(@. abs(pp - pole_locs[end, j+1]))
            ii = argmin(@. abs(pp - pole_locs[j+1]))
            #pole_locs[end, j] = pp[ii]
            pole_locs[j] = pp[ii]
        catch
            pole_locs[j] = pole_locs[j+1]
        end
        if abs(pole_locs[j] - pole_locs[j+1]) > 0.1
            pole_locs[j] = pole_locs[j+1]
        end
    end
    return t_vals, pole_locs
end