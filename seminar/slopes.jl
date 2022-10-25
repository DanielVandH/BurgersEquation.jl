using CairoMakie
using .BurgersEquation

μ = [0.0, 0.05, 0.1, 1.0]
t = LinRange(1e-6, 5.0, 1000)
#slopes = maximum_slope(t, μ)
indices, vals = findmaxima(slopes)
fig = Figure(fontsize=38)
ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"Maximum$ $ absolute slope",
    xticks=([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [L"0", L"1", L"2", L"3", L"4", L"5"]), height=400, width=600,
    yticks=([0.0, 2.0, 4.0, 6.0, 8.0], [L"0", L"2", L"4", L"6", L"8"]))
lines!(ax, t[1:5:283], slopes[1:5:283, 1], color=:black, linestyle=:dash)
lines!(ax, t, slopes[:, 2], color=:blue, linestyle=:solid)
lines!(ax, t, slopes[:, 3], color=:red, linestyle=:solid)
lines!(ax, t, slopes[:, 4], color=:darkgreen, linestyle=:solid)
scatter!(ax, [t[indices[2]]], [vals[2]], color=:black)
scatter!(ax, [t[indices[3]]], [vals[3]], color=:black)
xlims!(ax, 0, 5.0)
ylims!(ax, 0, 3.5)
text!(ax, L"\mu = 0.0", position=(1.02, 1.3), rotation=1.0, color=:black, textsize=36.4)
text!(ax, L"\mu = 0.05", position=(3.0, 1.3), rotation=0.0, color=:black, textsize=36.4)
text!(ax, L"\mu = 0.1", position=(2.2, 0.7), rotation=-0.05, color=:black, textsize=36.4)
text!(ax, L"\mu = 1.0", position=(1.5, 0.17), rotation=0.0, color=:black, textsize=36.4)
save("seminar/figures/slopes.pdf", fig)
save("seminar/figures/slopes.png", fig)
fig