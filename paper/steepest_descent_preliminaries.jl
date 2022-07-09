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

