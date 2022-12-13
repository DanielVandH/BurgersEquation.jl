x = LinRange(-200.0, 200.0, 250)
y = LinRange(0, 200.0, 250)
t = [100.0, 250.0, 500.0, 1000.0]
μ = [1.0]
u_vals = viscous_solution(x, y, t, μ)

fig = Figure(fontsize=38, resolution=(2970.0f0, 551.9453f0))
portrait!(fig, x, y, u_vals[:, :, 1, 1], 1, 1;
    xlabel=L"x", ylabel=L"y", title=L"(a): $t = 100$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
arrows!(fig.content[end], [-25.0], [25.0], [40.0], [0.0], color=:white, linewidth=4)
text!(fig.content[end], [-10.0], [25.0]; text=L"s_0^{(1)}", textsize=48, color=:black)
arrows!(fig.content[end], [100.0], [15.0], [-40.0], [0.0], color=:black, linewidth=4)
portrait!(fig, x, y, u_vals[:, :, 2, 1], 1, 2;
    xlabel=L"x", ylabel=L"y", title=L"(b): $t = 250$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
arrows!(fig.content[end], [-17.0], [40.0], [40.0], [0.0], color=:white, linewidth=4)
text!(fig.content[end], [-30.0], [10.0]; text=L"s_0^{(1)}", textsize=48, color=:black)
arrows!(fig.content[end], [130.0], [25.0], [-40.0], [0.0], color=:black, linewidth=4)
portrait!(fig, x, y, u_vals[:, :, 3, 1], 1, 3;
    xlabel=L"x", ylabel=L"y", title=L"(c): $t = 500$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
arrows!(fig.content[end], [-5.0], [55.0], [40.0], [0.0], color=:white, linewidth=4)
text!(fig.content[end], [-10.0], [20.0]; text=L"s_0^{(1)}", textsize=48, color=:black)
arrows!(fig.content[end], [170.0], [33.0], [-40.0], [0.0], color=:black, linewidth=4)
portrait!(fig, x, y, u_vals[:, :, 4, 1], 1, 4;
    xlabel=L"x", ylabel=L"y", title=L"(d): $t = 1000$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
arrows!(fig.content[end], [16.0], [75.0], [40.0], [0.0], color=:white, linewidth=4)
text!(fig.content[end], [10.0], [40.0]; text=L"s_0^{(1)}", textsize=48, color=:black)
arrows!(fig.content[end], [120.0], [25.0], [50.0], [15.0], color=:black, linewidth=4)

fig

save("$FIGURES/small_to_large_transition.$EXTENSION", fig)
