using BurgersEquation
using CairoMakie
using LinearAlgebra
using FastGaussQuadrature
using Contour
using OrderedCollections
using Peaks
using DelimitedFiles
using ArbNumerics
using JLD2
using ForwardDiff
using Optim
using HypergeometricFunctions
using DataFrames
using CSV

include("constants.jl")
include("introduction.jl")
include("exact_solution.jl")
include("colour_wheel.jl")
include("parabolic_cylinder_function.jl")
include("steepest_descent_preliminaries.jl")
include("tracking_poles.jl")
include("balancing_advection_diffusion.jl")
include("large_time_solution.jl")
include("large_time_solution_roots.jl")
include("aaa.jl")
include("aaa_closest_poles.jl")
# include("aaa_tracking_poles.jl") < --- See tracking_poles.jl
include("small_to_large_transition.jl")

## Some extra figures and results added after review
# Other initial conditions
include("other_initial_conditions.jl")

# Enstrophy computations 
include("enstrophy.jl")