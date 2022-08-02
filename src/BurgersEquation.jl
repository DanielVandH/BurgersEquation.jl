module BurgersEquation 

using LinearAlgebra
using Makie 
using CairoMakie
using FastGaussQuadrature
using HypergeometricFunctions
using ForwardDiff
using SpecialFunctions
using Colors 
using BaryRational
using ComplexPortraits
using LaTeXStrings
using StatsBase 
using Contour 
using ColorSchemes 
using JLD2 
using Peaks 
using Optim 
using OrderedCollections
using DifferentialEquations
using ArbNumerics

include("rootfinding.jl");          export newton_method, twopoint_linesearch
include("inviscid.jl");             export inviscid_solution, inviscid_singularities
include("quadrature.jl");           export gauss_legendre_finite, gauss_legendre_left_infinite, gauss_legendre_right_infinite, gauss_hermite
include("viscous.jl");              export viscous_solution
include("plotting.jl");             export portrait!, landscape!, LINSPECER_12_J, LINSPECER_250_J
include("similarity.jl");           export parabolicU, Φ₀, large_ξ_roots_θ
include("large_time.jl");           export viscous_solution_large_time, viscous_solution_large_time_Ψ, viscous_solution_large_time_Ψ_roots_θ, viscous_solution_large_time_Ψ_roots_ρ
include("steepest_descents.jl");    export saddle_point_approximation_μ, cubic_caustic, cubic_discriminant, cubic_saddle_points
include("split_functions.jl");      export complex_split_denominator, saddle_point_approximation_μ_split, Φ₀_split
include("pole_tracking.jl");        export tracking_poles_exact, tracking_poles_saddle, locate_pole, tracking_poles_Φ₀, Φ₀_pole, tracking_poles_aaa, tracking_poles_Ψ
include("slope_analysis.jl");       export maximum_slope
include("numerical_soln.jl");       export viscous_solution_finite_diff
include("aaa.jl");                  export viscous_solution_aaa

const f = x -> 1 / (1 + x^2)
const ∫f = s -> atan(s)
const tₛ = 8sqrt(3) / 9
export f, ∫f, tₛ

end