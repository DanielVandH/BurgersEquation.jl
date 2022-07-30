bt=2048
x1 = collect(-6.0:0.01:6.0); x2 = collect(-6.0:0.01:6.0); x3 = collect(-6.0:0.01:6.0); x4 = ArbReal.(-6:0.05:6, bits=bt)
y1 = collect(-6.0:0.01:6.0); y2 = collect(-6.0:0.01:6.0); y3 = collect(-6.0:0.01:6.0); y4 = ArbReal.(-6:0.05:6, bits=bt)
μ1 = 2.0; μ2 = 1.0; μ3 = 0.5; μ4 = ArbReal(0.1, bits = bt)
Φ₀_vals = [Φ₀(x1, y1, μ1), Φ₀(x2, y2, μ2), Φ₀(x3, y3, μ3), Φ₀(x4, y4, μ4)]
t = 1e-6
z = Vector{Matrix{Union{ComplexF64,ArbComplex{P}} where {P}}}(undef, 4)
u = Vector{Matrix{Union{ComplexF64,ArbComplex{P}} where {P}}}(undef, 4)
num_nodes = 250
nodes, weights = gausshermite(num_nodes)
glnodes, glweights = gausslegendre(num_nodes)
x = [x1, x2, x3, x4]
y = [y1, y2, y3, y4]
μ = [μ1, μ2, μ3, μ4]
for j in 1:4
    xx = Float64.(x[j])
    yy = Float64.(y[j])
    z[j] = [im .+ sqrt(t) * (xx + im * yy) for xx in xx, yy in yy]
    u[j] = viscous_solution.(real(z[j]), imag(z[j]), t, Float64(μ[j]), Ref(nodes), Ref(weights), Ref(glnodes), Ref(glweights))
end

fig = Figure(fontsize=38, resolution=(2000, 1600))
plot_layout = [[(1, 1), (2, 1), (3, 1), (4, 1)], [(1, 2), (2, 2), (3, 2), (4, 2)],
    [(1, 3), (2, 3), (3, 3), (4, 3)], [(1, 4), (2, 4), (3, 4), (4, 4)]]
for j = 1:4
    portrait!(fig, x[j], y[j], Complex{Float64}.(Φ₀_vals[j]), plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
        title=L"(%$(ALPHABET[j])): $\mu = %$(Float64(μ[j]))$", titlealign=:left,
        xlabel=L"\mathrm{Re}(\xi)", ylabel=L"\mathrm{Im}(\xi)", width=450, height=450,
        xticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]),
        yticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]))

    portrait!(fig, real.(z[j]), imag.(z[j]), Complex{Float64}.(u[j]), plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
        title=L"(%$(ALPHABET[j+8])): $\mu = %$(Float64(μ[j]))$, $t = 10^{-6}$", titlealign=:left,
        xlabel=L"x", ylabel=L"y", width=450, height=450,
        xticks=([-0.005, 0.0, 0.005], [L"-0.005", L"0", L"0.005"]),
        yticks=([0.995, 1.0, 1.005], [L"0.995", L"1", L"1.005"]))
end
resize_to_layout!(fig)

fig

save("$FIGURES/similarity_solutions.$EXTENSION", fig)
