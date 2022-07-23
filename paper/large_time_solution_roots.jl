μ = LinRange(0.01, 25.0, 2500)
z₁ = 4.0 + 22.0im
z₄ = 4.0 - 22.0im
z₃ = -10.5 - 25.0im
z₂ = -10.5 + 25.0im
roots₁ = tracking_poles_Ψ(z₁, μ)
roots₂ = tracking_poles_Ψ(z₂, μ)
roots₃ = tracking_poles_Ψ(z₃, μ)
roots₄ = tracking_poles_Ψ(z₄, μ)
z = [z₁, z₂, z₃, z₄]
rootsᵢ = [roots₁, roots₂, roots₃, roots₄]

fig = Figure(fontsize=36, resolution=(1600, 800))
ax = Axis(fig[1, 1], width=400, height=400,
    xlabel=L"\mathrm{Re}(\eta)",
    ylabel=L"\mathrm{Im}(\eta)",
    xticks=([-10, -7.5, -5.0, -2.5, 0, 2.5, 5.0], [L"-10", L"-7.5", L"-5", L"-2.5", L"0", L"2.5", L"5"]),
    yticks=([-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25], [L"-25", L"-20", L"-15", L"-10", L"-5", L"0", L"5", L"10", L"15", L"20", L"25"]),
    title=L"(a):$ $ Pole distance",
    titlealign=:left)
idx = 1:10:2500
colors = cgrad(LINSPECER_250_J, μ[idx] / 25, categorical=false)
[scatter!(ax, [real(rootsᵢ[j][idx[i]])], [imag(rootsᵢ[j][idx[i]])], color=colors[i], markersize=4) for j in 1:4, i in 1:length(idx)]

ax = Axis(fig[1, 2], width=400, height=400,
    xlabel=L"\mu",
    ylabel=L"\mathrm{Im}(\eta)",
    xticks=([0, 5, 10, 15, 20, 25], [L"0", L"5", L"10", L"15", L"20", L"25"]),
    yticks=([-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25], [L"-25", L"-20", L"-15", L"-10", L"-5", L"0", L"5", L"10", L"15", L"20", L"25"]),
    title=L"(b):$ $ Pole locations",
    titlealign=:left)
_colors = [:red, :black, :blue, :purple]
[scatter!(ax, [μ[idx[i]]], [imag(rootsᵢ[j][idx[i]])], color=_colors[j], markersize=4) for j in 1:4, i in 1:length(idx)]#not going to bother using lines! for this one
text!(ax, L"Quadrant $1$", position=(5, 13), color=:black, textsize=36.4, rotation = 0.3)
text!(ax, L"Quadrant $2$", position=(7, 4.5), color=:black, textsize=36.4, rotation = 0.3)
text!(ax, L"Quadrant $3$", position=(4, -17), color=:black, textsize=36.4, rotation = -0.33)
text!(ax, L"Quadrant $4$", position=(7, -9.7), color=:black, textsize=36.4, rotation = -0.33)

Colorbar(fig[1, 3], limits=(0.0, 25.0), colormap=colors, label=L"\mu", vertical=true,
    ticks=([0.0, 5.0, 10.0, 15.0, 20.0, 25.0], [L"0", L"5", L"10", L"15", L"20", L"25"]))
resize_to_layout!(fig)

save("$FIGURES/large_time_solution_plots_all_roots.$EXTENSION", fig)

###########################################
## Closest pole
########################################### 

x = LinRange(-5, 5, 500)
y = LinRange(-5, 5, 500)
μ = [0.1, 1.0]
vals = viscous_solution_large_time_Ψ(x, y, μ)

fig = Figure()
portrait!(fig, x, y, vals[:, :, 1], 1, 1; width=400,height=400,aspect=1)
scatter!(fig.content[end], [2.07,2.07,-4.14],[-0.9238,0.9238,0.0], color=:black, markersize=13)
portrait!(fig, x, y, vals[:, :, 2], 1, 2; width=400,height=400,aspect=1)
scatter!(fig.content[end], [2.201,2.201,-4.401],[1.59,-1.59,0.0], color=:black, markersize=13)
resize_to_layout!(fig)
fig
save("$FIGURES/large_time_solution_plots_all_roots_closest_cubic.$EXTENSION", fig)
