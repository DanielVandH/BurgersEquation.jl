# BurgersEquation.jl

This package contains code for reproducing the numerical results in our paper (...). Installation instructions are given below, and a description of the contents follow.

# Installation 
You can install the package into Julia using 

```
using Pkg 
Pkg.add(url = "https://github.com/DanielVandH/BurgersEquation.jl")
using BurgersEquation
```

although this just gives you access to the functions in `src`, which may not be so useful on their own. Another option is to clone the repository into your `dev` environment and use the package locally, allowing you access to the package and to the code in the `paper` directory of the repository. To do this, you can do the following (tested in VS Code on a Windows 10 machine)
```
julia> ] develop https://github.com/DanielVandH/BurgersEquation.jl
```
This code will create a package in the `dev` directory of your Julia environment. On Windows, this folder is located at `C:\Users\USER\.julia\dev\BurgersEquation`. You can then go into VS Code and change your folder to this new directory, and you should have complete ability to work with the code in this repository. If it's worked correctly, then your environment should be `BurgersEquation`; your current active environment can be seen by running
```
julia> Base.active_project()
```
On my machine, this gives
```
julia> Base.active_project()
"c:\\Users\\USER\\.julia\\dev\\BurgersEquation\\Project.toml"
```
It is important that you have this environment activated so that you have access to the same Julia packages used in the paper rather than new, possibly breaking, ones.

There may be some issues with phase portraits and LaTeX tick labels in 3D plots. If this is the case, then try running
```
using Pkg
Pkg.add(url = "https://github.com/luchr/ComplexPortraits.jl", rev = "master")
Pkg.add(url = "https://github.com/DanielVandH/Makie.jl")
```
and then restart your Julia environment and try again. This second call installs my fork of `Makie.jl`, implementing a fix for the tick label issues. If even that fails, perhaps try
```
using Pkg 
Pkg.instantiate()
```
You may also need to install my `BaryRational.jl` fork,
```
using Pkg
Pkg.add(url = "https://github.com/DanielVandH/BaryRational.jl")
end
```
although none of this code is actually called directly (but you can call it if you want to, see the `/src/` description below).

If you want to run any of the MATLAB code (none of which is called directly), you also need Chebfun.

# Contents 
The two main folders in the repository are `paper` and `src`.

- `/src/`

  This folder contains all the `using` and `include` commands that define the functions that do all the heavy lifting in the paper. Most of the functions used are documented with docstrings. For example, if you want clearer documentation on our functions for plotting portraits, you can write
  ```
  julia> ?portrait!
  search: portrait!

  portrait!(fig, x, y, Z, i, j; xlabel=L"x", ylabel=L"y", nist=true, pltkwargs...)

  Plots the phase portrait for the data in Z, defined over a structured grid.

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •  fig: A Makie figure object.

    •  x: The real part of the points that Z is computed at.

    •  y: The imaginary part of the points that Z is computed at.

    •  Z: The matrix of values to plot, defined such that Z[i, j] is at (x[i], y[j]).

    •  i: The row number to create the axis in fig.

    •  j: The column number to create the axis in fig.

  Keyword Arguments
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    •  xlabel=L"x": Label for the x-axis.

    •  ylabel=L"y": Label for the y-axis.

    •  nist = true: Whether to use the NIST color scheme. If false, the HSV color scheme is used.

    •  pltkwargs...: Extra keyword arguments for Axis.

  Output
  ≡≡≡≡≡≡≡≡

  There is no output, but fig is updated with the figure in fig[i, j].
  ```
  A full list of functions exported from the package (many are not exported, simply performing intermediate operations) can be seen using `names`:
  ```
  julia> names(BurgersEquation)
  40-element Vector{Symbol}:
  :BurgersEquation
  :LINSPECER_12_J 
  :LINSPECER_250_J
  :complex_split_denominator
  :cubic_caustic
  :cubic_discriminant
  :cubic_saddle_points
  :f
  :gauss_hermite
  ⋮
  :viscous_solution_large_time
  :viscous_solution_large_time_Ψ
  :viscous_solution_large_time_Ψ_roots_θ
  :viscous_solution_large_time_Ψ_roots_ρ
  :Φ₀
  :Φ₀_pole
  :Φ₀_split
  :∫f
  ```
  Not all of these functions are used in the paper, for example `numerical_solution.jl` is not used in the paper, but contains some methods for solving Burgers' equation numerically with finite differences by transforming to the heat equation first with the Cole-Hopf transformation. Similarly, `aaa.jl` is not used in the paper but you can use it instead for producing the AAA results (with less accuracy, unfortunately) with Julia rather than with Chebfun in MATLAB (described below). How we would do this is given in the commented code in the AAA files of the `/paper/` code, described below.
  
- `/paper/`

  This folder contains the code that actually produces the results in the paper. The main script in this filder is `paper_code.jl` which actually runs all the other scripts. In order of how we include the files in `paper_code.jl`, a description of each script is as follows:
    - `constants.jl`: This is a script which defines some constants used in the code. The first constants are just setting up values for saving figures. The second set of constants comes from the MATLAB script `burger_aaa.m` which we use to make some AAA plots. The constants are given here so that you do not need MATLAB to look at the plots. Corresponding Julia code is given in `aaa.jl`, but commented out. The third set of constants comes from the MATLAB script `burger_tracking_poles.m` which we use to obtain the data for tracking poles with the AAA algorithm. A Julia version for obtaining this data is given in `aaa_tracking_poles.jl`, but commented out. All the data from these files is saved in `/paper/data/`. The MATLAB script `fix_i.m` is used to change `i` to `im` in the data files.
    - `introduction.jl`: This contains code for producing the figure in our introduction. In particular, a figure is produced which plots the inviscid solution at the shock time and the viscous solution at this shock time and other times.
    - `exact_solution.jl`: This contains code for plotting the exact solution.
    - `colour_wheel.jl`: This plots the colour wheel that we use for visualising solutions in the complex plane.
    - `parabolic_cylinder_function.jl`: This contains code for plotting the small-time similarity solution.
    - `steepest_descent_preliminaries.jl`: This contains code for two figures. First, we make the plot for the values of the discriminant, showing the nature of the saddle points in different regions. We then make a plot which shouts the contour that we integrate over when applying the saddle point method, along with the surface to explain the issue with branch cuts.
    - `tracking_poles.jl`: This contains code for two figures. The first figure shows results for tracking poles with the saddle point method for small $\mu$, which we also compare to our method of tracking poles from the exact solution. We then also write `include("aaa_tracking_poles.jl")` which calls the code in `aaa_tracking_poles.jl` that has results for tracking these poles with the AAA algorithm, making use of the data from `constants.jl`. The second figure produced shows results for tracking poles for larger values of $\mu$, although these are not shown in the paper.
    - `balancing_advection_diffusion.jl`: This contains code for our experiments on how we can relate the steepest of the solution on the real line to the proximity of a pole in the complex plane to the real line. In particular, the code is used to try and find the value $\mu$ that balances the effect of diffusion and advection. The first figure produced from this experiment directly compares these slopes and proximities, and the value in `breakdown_μ` shows the value of $\mu$ that balances these effects according to the slopes, and `μb` is the value according to the proximities. This value of `μb` is not reliably computed, though, since we aren't using specialised methods for highly oscillatory integrals. The second figure produced here shows how we obtain the value of $\mu$ above which there are no longer any poles in the lower-half plane in the similarity solution, with the value given in `μˢ`. The second figure produced here is not shown in the paper.
    - `large_time_solution.jl`: This script contains code for plotting the large-time solution, along with the asymptotic results for the poles.
    - `large_time_solution_roots.jl`: This script contains code for tracking the roots in the large-time solution, although none of this is shown in the paper.
    - `aaa.jl`: This script contains code for plotting the AAA approximants, using the data from `constants.jl` obtained from the MATLAB script `burger_aaa.m`. The corresponding Julia code for producing these results is commented out at the top if you are interested.
    - `aaa_closest_poles.jl`: This script produces a figure comparing the AAA approximant and the exact solution, making the point that the closest poles are quite accurate. The data used for these figures comes from `constants.jl` obtained from the MATLAB script `burger_aaa.m`, but the corresponding Julia code to make the figures is commented out at the top.
    - `aaa_tracking_poles.jl`: This script contains code for plotting the poles tracked from the AAA approximants, using the data from `constants.jl` obtained from the MATLAB script `burger_tracking_poles.m`. This is included when running `tracking_poles.jl`. The corresponding Julia code that could be used to make these figures is commented out at the top.
    - `small_to_large_transition.jl`: This script contains code for looking at how the solution evolves for large time, noting that the pole-zero pairs get pushed apart to a separation distance of infinity.

Provided you have installed the package correctly, there should be no errors if you simply write `include("paper_code.jl")` in the REPL. If you uncomment out the Julia code for the AAA results, you need to install `BaryRational.jl`.
