x = -1:0.001:1
y = -1:0.001:1
z = [x + im * y for x in x, y in y]
z[abs.(z).>1] .= NaN .+ im * NaN
fig = Figure(fontsize=38, resolution=(600, 600))
portrait!(fig, x, y, z, 1, 1; nist=NIST,
    xticks=([-1.0, -0.5, 0.0, 0.5, 1.0], [L"-1", L"-0.5", L"0", L"0.5", L"1"]),
    yticks=([-1.0, -0.5, 0.0, 0.5, 1.0], [L"-1", L"-0.5", L"0", L"0.5", L"1"]))
save("$FIGURES/colour_wheel.$EXTENSION", fig)