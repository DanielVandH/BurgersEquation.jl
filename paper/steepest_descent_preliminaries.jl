## Discriminant plot 
x = -5:0.01:5
t = 0:0.01:5
X = vec(x' .* ones(length(t)))
T = vec(ones(length(x))' .* t)
discriminant_vals = vec(cubic_discriminant.(X, T))
xc = 1.733:0.01:5
caustic_vals1 = zeros(length(xc))
caustic_vals2 = zeros(length(xc))
for (i, x) in enumerate(xc)
    caustic_vals1[i], caustic_vals2[i] = cubic_caustic(x)
end

fig = Figure(fontsize=26, resolution=(800, 400))
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"t",
    xticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]),
    yticks=([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [L"0", L"1", L"2", L"3", L"4", L"5"]))
scatter!(ax, X[discriminant_vals.<0.0], T[discriminant_vals.<0.0], color=:blue, markersize=4)
scatter!(ax, X[discriminant_vals.>0.0], T[discriminant_vals.>0.0], color=:red, markersize=4)
lines!(ax, xc, caustic_vals1, color=:black, linewidth=6)
lines!(ax, xc, caustic_vals2, color=:black, linewidth=6)
xlims!(ax, -5, 5)
ylims!(ax, 0, 5)
text!(ax, L"\Delta < 0", position=(-2.0, 2.0), color=:white, textsize=36.4)
text!(ax, L"\Delta > 0", position=(3.0, 4.0), color=:white, textsize=36.4)
scatter!(ax, [sqrt(3)], [tₛ], color=:black, markersize=9)
scatter!(ax, [sqrt(3)], [tₛ], color=:white, markersize=6)

save("$FIGURES/discriminant_plot.$EXTENSION", fig)

## Contour 
x = 2.0713 + 0.48208im
t = 2.0
s₁, s₂, s₃ = cubic_saddle_points(x, t)

h(s) = -0.5atan(s) - (x - s)^2 / (4t)
h(u, v) = h(u + im * v)
ℜh(u, v) = real(h(u, v))
ℑh(u, v) = imag(h(u, v))
ℑh(u, v, k) = ℑh(u, v) - k

# Saddle points and other parameters 
S = [s₁, s₂, s₃]
k = imag.(h.(S))
PTS = [-1.0 0.0; 1.0 0.0]
sadd = [real.(S) imag.(S)]
domain = [-2.0 2.0; -2.0 2.0]
colors = cgrad(LINSPECER_250_J, LinRange(-1, 1, 100), categorical=false)

fig = Figure(fontsize=26, resolution=(800, 400))
ga = fig[1, 1:2] = GridLayout()

# Contour plot 
ax = Axis(ga[1, 2], xlabel=L"\mathrm{Re}(s)", ylabel=L"\mathrm{Im}(s)",
    title=L"(b)$:$ Contours", titlealign=:left,
    xticks=([-2.0, -1.0, 0.0, 1.0, 2.0], [L"-2", L"-1", L"0", L"1", L"2"]),
    yticks=([-2.0, -1.0, 0.0, 1.0, 2.0], [L"-2", L"-1", L"0", L"1", L"2"]))
u = domain[1, 1]:0.01:domain[1, 2]
v = domain[2, 1]:0.01:domain[2, 2]
z = [ℜh(u, v) for u in u, v in v]
contourf!(ax, u, v, z, colorrange=(-1, 1), colormap=colors)
z₂ = [ℑh(u, v) for u in u, v in v]
contour!(ax, u, v, z₂, levels=sort(k), linewidth=3, color=[:blue, :magenta, :red])
scatter!(ax, real.(S), imag.(S), markersize=12, color=:black)
scatter!(ax, real.(S), imag.(S), markersize=8, color=:white)
xlims!(ax, -2, 2)
ylims!(ax, -2, 2)
X = [-2.36842 -1.9378; -1.33971 -0.909091; -0.263158 0.0239234; 0.31105 0.263158; 0.76555 0.956938; 1.41148 1.9378];
Y = [0.287081 0.215311; 0.191388 0.167464; 0.239234 0.358852; 0.669856 1.29187; 1.10048 0.980861; 0.837321 0.669856];
arrows!(ax, X[:, 1], Y[:, 1], X[:, 2] .- X[:, 1], Y[:, 2] .- Y[:, 1], color=:white, arrowsize=13.0, linewidth=3.5)

# Surface plot 
ax = Axis3(ga[1, 1], xlabel=L"\mathrm{Re}(s)", ylabel=L"\mathrm{Im}(s)",
    zlabel=L"\mathrm{Re}(h)", title=L"(a)$:$ Surface", titlealign=:left,
    azimuth=12.0, elevation=0.65,
    xticks=([-2.0, -1.0, 0.0, 1.0, 2.0], [L"-2", L"-1", L"0", L"1", L"2"]),
    yticks=([-2.0, -1.0, 0.0, 1.0, 2.0], [L"-2", L"-1", L"0", L"1.", L"2"]),
    zticks=([-1.0, 0.0, 1.0], [L"-1", L"0", L"1"]))
surface!(ax, u, v, z, color=z, colorrange=(-1, 1), colormap=colors, overdraw=true)
C = contours(u, v, z₂, sort(k))
contour_colors = [:blue, :magenta, :red]
for (k, cl) in enumerate(levels(C))
    for line in Contour.lines(cl)
        us, vs = coordinates(line)
        lines!(ax, us, vs, ℜh.(us, vs), color=contour_colors[k], overdraw=false)
    end
end
colgap!(ga, 90)
save("$FIGURES/contour_surface_plot.$EXTENSION", fig)

## Is the pair of saddle points the same for all poles z(t)?

# Need to start by computing the poles 
t = [2.0, 2.0, 2.0]
μ = [0.05, 0.1468, 0.5]
z = [2.0713428071028632 + 0.48208285060903094im,
    2.125956790502127 + 1.1121784646136845im,
    2.3122739258833174 + 2.711309573638537im]
t_vals_exact, pole_locs_exact = tracking_poles_exact(z, t, μ) # See tracking_poles.jl 
outer_size = length(t_vals_exact[1]) # same for 2, 3 
spacing = 50
_t_vals_exact = [t_vals_exact[i][j] for j in 1:spacing:outer_size, i in eachindex(μ)]
_pole_locs_exact = [pole_locs_exact[i][j] for j in 1:spacing:outer_size, i in eachindex(μ)]
outer_size = size(_pole_locs_exact, 1) # same for 2, 3

# Now, for all the poles, we need to compute the saddle points (z₁, z₂, z₃)
s₁ = deepcopy(_pole_locs_exact)
s₂ = deepcopy(_pole_locs_exact)
s₃ = deepcopy(_pole_locs_exact)
for i in eachindex(μ)
    for j in 1:outer_size
        if j == 1
            s₁[j, i], s₂[j, i], s₃[j, i] = ComplexF64(NaN, NaN), ComplexF64(NaN, NaN), ComplexF64(NaN, NaN)
        else
            _s₁, _s₂, _s₃ = cubic_saddle_points(_pole_locs_exact[j, i], _t_vals_exact[j, i])
            h(s) = -0.5atan(s) - (_pole_locs_exact[j, i] - s)^2 / (4.0 * _t_vals_exact[j, i])
            h(u, v) = h(u + im * v)
            Diff₁₂ = abs(real(h(_s₁)) - real(h(_s₂)))
            Diff₁₃ = abs(real(h(_s₁)) - real(h(_s₃)))
            Diff₂₃ = abs(real(h(_s₂)) - real(h(_s₃)))
            Diffs = [Diff₁₂, Diff₁₃, Diff₂₃]
            idx = findmin(Diffs)[2]
            if idx == 1
                s₁[j, i], s₂[j, i], s₃[j, i] = _s₁, _s₂, _s₃
            elseif idx == 2
                s₁[j, i], s₂[j, i], s₃[j, i] = _s₁, _s₃, _s₂
            elseif idx == 3
                s₁[j, i], s₂[j, i], s₃[j, i] = _s₂, _s₃, _s₁
            end
        end
    end
end

# Now, for all the poles, we need to compute the contour data 
domain = [-4.0 4.0; -4.0 4.0]
colors = cgrad(LINSPECER_250_J, LinRange(-1, 1, 100), categorical=false)
u = domain[1, 1]:0.01:domain[1, 2]
v = domain[2, 1]:0.01:domain[2, 2]
contour_real = [zeros(length(u), length(v)) for _ in 1:outer_size, _ in eachindex(μ)]
contour_imag = [zeros(length(u), length(v)) for _ in 1:outer_size, _ in eachindex(μ)]
for i in eachindex(μ)
    for j in 1:outer_size
        h(s) = -0.5atan(s) - (_pole_locs_exact[j, i] - s)^2 / (4.0 * _t_vals_exact[j, i])
        h(u, v) = h(u + im * v)
        ℜh(u, v) = j == 1 ? NaN : real(h(u, v))
        ℑh(u, v) = j == 1 ? NaN : imag(h(u, v))
        contour_real[j, i] .= [ℜh(u, v) for u in u, v in v]
        contour_imag[j, i] .= [ℑh(u, v) for u in u, v in v]
    end
end

# Next, compute all the levels k 
k = [Float64[] for _ in axes(_pole_locs_exact, 1), _ in axes(_pole_locs_exact, 2)]
for i in eachindex(μ)
    for j in 1:outer_size
        h(s) = -0.5atan(s) - (_pole_locs_exact[j, i] - s)^2 / (4.0 * _t_vals_exact[j, i])
        _s₁, _s₂, _s₃ = s₁[j, i], s₂[j, i], s₃[j, i]
        push!(k[j, i], [imag(h(_s₁)), imag(h(_s₂)), imag(h(_s₃))]...)
    end
end

# Prepare the axes 
j = Observable(2)
fig = Figure(fontsize=48, resolution=(2270.2766f0, 779.6094f0))

ax1 = Axis(fig[1, 1], xlabel=L"\mathrm{Re}(s)", ylabel=L"\mathrm{Im}(s)",
    title=Makie.lift(_j -> L"(a): $\mu = %$(μ[1])$, $t = %$(rpad(round(_t_vals_exact[_j, 1], digits = 3), 5, '0'))$", j),
    titlealign=:left,
    xticks=(-4:2:4, [L"%$s" for s in -4:2:4]),
    yticks=(-4:2:4, [L"%$s" for s in -4:2:4]),
    width=600,
    height=600,
    aspect=1)
ax2 = Axis(fig[1, 2], xlabel=L"\mathrm{Re}(s)", ylabel=L"\mathrm{Im}(s)",
    title=Makie.lift(_j -> L"(b): $\mu = %$(μ[2])$, $t = %$(rpad(round(_t_vals_exact[_j, 2], digits = 3), 5, '0'))$", j),
    titlealign=:left,
    xticks=(-4:2:4, [L"%$s" for s in -4:2:4]),
    yticks=(-4:2:4, [L"%$s" for s in -4:2:4]),
    width=600,
    height=600,
    aspect=1)
ax3 = Axis(fig[1, 3], xlabel=L"\mathrm{Re}(s)", ylabel=L"\mathrm{Im}(s)",
    title=Makie.lift(_j -> L"(c): $\mu = %$(μ[3])$, $t = %$(rpad(round(_t_vals_exact[_j, 3], digits = 3), 5, '0'))$", j),
    titlealign=:left,
    xticks=(-4:2:4, [L"%$s" for s in -4:2:4]),
    yticks=(-4:2:4, [L"%$s" for s in -4:2:4]),
    width=600,
    height=600,
    aspect=1)

contourf!(ax1, u, v, @lift(contour_real[$j, 1]), colorrange=(-1, 1), colormap=colors)
contourf!(ax2, u, v, @lift(contour_real[$j, 2]), colorrange=(-1, 1), colormap=colors)
contourf!(ax3, u, v, @lift(contour_real[$j, 3]), colorrange=(-1, 1), colormap=colors)

for i in 1:3
    contour!(ax1, u, v, @lift(contour_imag[$j, 1]), levels=@lift([k[$j, 1][i]]), linewidth=3, color=[:magenta, :red, :blue][i])
    contour!(ax2, u, v, @lift(contour_imag[$j, 2]), levels=@lift([k[$j, 2][i]]), linewidth=3, color=[:magenta, :red, :blue][i])
    contour!(ax3, u, v, @lift(contour_imag[$j, 3]), levels=@lift([k[$j, 3][i]]), linewidth=3, color=[:magenta, :red, :blue][i])
end

scatter!(ax1, @lift([real(s₃[$j, 1])]), @lift([imag(s₃[$j, 1])]), markersize=12, color=:black)
scatter!(ax1, @lift([real(s₃[$j, 1])]), @lift([imag(s₃[$j, 1])]), markersize=8, color=:white)
scatter!(ax1, @lift([real(s₁[$j, 1]), real(s₂[$j, 1])]), @lift([imag(s₁[$j, 1]), imag(s₂[$j, 1])]), markersize=12, color=:black)
scatter!(ax1, @lift([real(s₁[$j, 1]), real(s₂[$j, 1])]), @lift([imag(s₁[$j, 1]), imag(s₂[$j, 1])]), markersize=8, color=:red)
scatter!(ax2, @lift([real(s₃[$j, 2])]), @lift([imag(s₃[$j, 2])]), markersize=12, color=:black)
scatter!(ax2, @lift([real(s₃[$j, 2])]), @lift([imag(s₃[$j, 2])]), markersize=8, color=:white)
scatter!(ax2, @lift([real(s₁[$j, 2]), real(s₂[$j, 2])]), @lift([imag(s₁[$j, 2]), imag(s₂[$j, 2])]), markersize=12, color=:black)
scatter!(ax2, @lift([real(s₁[$j, 2]), real(s₂[$j, 2])]), @lift([imag(s₁[$j, 2]), imag(s₂[$j, 2])]), markersize=8, color=:red)
scatter!(ax3, @lift([real(s₃[$j, 3])]), @lift([imag(s₃[$j, 3])]), markersize=12, color=:black)
scatter!(ax3, @lift([real(s₃[$j, 3])]), @lift([imag(s₃[$j, 3])]), markersize=8, color=:white)
scatter!(ax3, @lift([real(s₁[$j, 3]), real(s₂[$j, 3])]), @lift([imag(s₁[$j, 3]), imag(s₂[$j, 3])]), markersize=12, color=:black)
scatter!(ax3, @lift([real(s₁[$j, 3]), real(s₂[$j, 3])]), @lift([imag(s₁[$j, 3]), imag(s₂[$j, 3])]), markersize=8, color=:red)

# Animate 
j_rng = 2:outer_size
secs = 15
framerate = floor(Int64, outer_size / secs)
record(fig, "$FIGURES/saddle_point_contour_animation.mp4", j_rng; framerate=framerate) do _j
    @show _j
    j[] = _j
end
