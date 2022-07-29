x = [-6.0:0.01:6.0, -6.0:0.01:6.0, -6.0:0.01:6.0, -6:0.01:6]
y = [-6.0:0.01:6.0, -6.0:0.01:6.0, -6.0:0.01:6.0, -6:0.01:6]
μ = [2.0, 1.0, 0.5, 0.1]
Φ₀_vals = Φ₀.(x, y, μ)
t = 1e-6
z = Vector{Matrix{ComplexF64}}(undef, 4)
u = Vector{Matrix{ComplexF64}}(undef, 4)
num_nodes = 250
nodes, weights = gausshermite(num_nodes)
glnodes, glweights = gausslegendre(num_nodes)
for j in 1:4
    xx = x[j]
    yy = y[j]
    z[j] = [im .+ sqrt(t) * (xx + im * yy) for xx in xx, yy in yy]
    u[j] = viscous_solution.(real(z[j]), imag(z[j]), t, μ[j], Ref(nodes), Ref(weights), Ref(glnodes), Ref(glweights))
end

fig = Figure(fontsize=33, resolution=(2000, 1600))
plot_layout = [[(1, 1), (2, 1), (3, 1), (4, 1)], [(1, 2), (2, 2), (3, 2), (4, 2)],
    [(1, 3), (2, 3), (3, 3), (4, 3)], [(1, 4), (2, 4), (3, 4), (4, 4)]]
for j = 1:4
    portrait!(fig, x[j], y[j], Φ₀_vals[j], plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
        title=L"(%$(ALPHABET[j])): $\mu = %$(μ[j])$", titlealign=:left,
        xlabel=L"\mathrm{Re}(\xi)", ylabel=L"\mathrm{Im}(\xi)", width=450, height=450,
        xticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]),
        yticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]))

    portrait!(fig, real.(z[j]), imag.(z[j]), u[j], plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
        title=L"(%$(ALPHABET[j+8])): $\mu = %$(μ[j])$, $t = 10^{-6}$", titlealign=:left,
        xlabel=L"x", ylabel=L"y", width=450, height=450,
        xticks=([-0.005, 0.0, 0.005], [L"-0.005", L"0", L"0.005"]),
        yticks=([0.995, 1.0, 1.005], [L"0.995", L"1", L"1.005"]))
end
resize_to_layout!(fig)

save("$FIGURES/similarity_solutions.$EXTENSION", fig)
