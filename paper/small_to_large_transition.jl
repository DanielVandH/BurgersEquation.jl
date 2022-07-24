x = LinRange(-200.0, 200.0, 250)
y = LinRange(0, 200.0, 250)
t = [100.0, 250.0, 500.0, 1000.0]
μ = [1.0]
u_vals = viscous_solution(x, y, t, μ)

fig = Figure(fontsize=34, resolution=(3000.1049f0, 740.89844f0))
portrait!(fig, x, y, u_vals[:, :, 1, 1], 1, 1;
    xlabel=L"x", ylabel=L"y", title=L"(a): $t = 100$",
    width=600, height=600, aspect=1,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 2, 1], 1, 2;
    xlabel=L"x", ylabel=L"y", title=L"(b): $t = 250$",
    width=600, height=600, aspect=1,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 3, 1], 1, 3;
    xlabel=L"x", ylabel=L"y", title=L"(c): $t = 500$",
    width=600, height=600, aspect=1,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 4, 1], 1, 4;
    xlabel=L"x", ylabel=L"y", title=L"(d): $t = 1000$",
    width=600, height=600, aspect=1,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)

save("$FIGURES/small_to_large_transition.$EXTENSION", fig)
