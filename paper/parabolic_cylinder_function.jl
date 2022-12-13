bt = 2048
x1 = collect(-6.0:0.01:6.0)
x2 = collect(-6.0:0.01:6.0)
x3 = collect(-6.0:0.01:6.0)
x4 = ArbReal.(-6:0.05:6, bits=bt)
y1 = collect(-6.0:0.01:6.0)
y2 = collect(-6.0:0.01:6.0)
y3 = collect(-6.0:0.01:6.0)
y4 = ArbReal.(-6:0.05:6, bits=bt)
μ1 = 2.0
μ2 = 1.0
μ3 = 0.5
μ4 = ArbReal(0.1, bits=bt)
#Φ₀_vals = [Φ₀(x1, y1, μ1), Φ₀(x2, y2, μ2), Φ₀(x3, y3, μ3), Φ₀(x4, y4, μ4)]
#Φ₀_vals64 = Matrix{Complex{Float64}}.(Φ₀_vals)
@load "paper/data/Phi064vals.jld2"
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
reverse!(plot_layout)
for j = 1:4
    portrait!(fig, x[j], y[j], Complex{Float64}.(Φ₀_vals64[j]), plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
        title=L"(%$(ALPHABET[4-j+1])): $\mu = %$(Float64(μ[j]))$", titlealign=:left,
        xlabel=L"\mathrm{Re}(\xi)", ylabel=L"\mathrm{Im}(\xi)", width=450, height=450,
        xticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]),
        yticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]))

    portrait!(fig, real.(z[j]), imag.(z[j]), Complex{Float64}.(u[j]), plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
        title=L"(%$(ALPHABET[8-j+1])): $\mu = %$(Float64(μ[j]))$, $t = 10^{-6}$", titlealign=:left,
        xlabel=L"x", ylabel=L"y", width=450, height=450,
        xticks=([-0.005, 0.0, 0.005], [L"-0.005", L"0", L"0.005"]),
        yticks=([0.995, 1.0, 1.005], [L"0.995", L"1", L"1.005"]))
end
ξ1 = [complex(1.9, 4.1), [sqrt(8n * Float64(μ1) * π) for n in 2:100]...]
ξ2 = [complex(1.75, 2.5), complex(3.75, 4.4), complex(5.0, 5.6),
    [sqrt(8n * Float64(μ2) * π) for n in 4:100]...]
ξ3 = [complex(1.5, 1.2), complex(2.7, 2.8), complex(3.6, 3.75),
    complex(4.4, 4.5), complex(5.0, 5.2), complex(5.6, 5.7),
    [sqrt(8n * Float64(μ3) * π) for n in 7:100]...]
ξ4 = [complex(1.2, -0.1), complex(1.5, 0.5), complex(1.8, 1.0), complex(2.2, 1.4),
    complex(2.4, 1.8),
    complex(2.6, 2.0), complex(2.8, 2.2),
    complex(3.0, 2.5), complex(3.2, 2.7),  complex(3.4, 3.0),
    complex(3.6, 3.1), complex(3.7, 3.3), complex(3.9, 3.5),
    complex(4.025, 3.65), [sqrt(8n * Float64(μ4) * π) for n in 13:100]...]
ξ = [ξ1, ξ2, ξ3, ξ4]
ξs = [large_ξ_roots_ρ.(abs.(ξ[i]), Float64(μ[i])) for i in 1:4]
ρ = 1.0:0.01:200.0
θs = [large_ξ_roots_θ.(ρ, Float64(μ[i]), 1) for i in 1:4]
ρs = [@. ρ * exp(im * θs[i]) for i in 1:4]
ξs[4][1]=complex(NaN,NaN)
ξs[4][7]=2.84+2.2825im
ξs[4][12] = 3.728+3.346im
push!(ξs[4], 5.230+4.935im)
[lines!(fig.content[i], real(ρs[j]), imag(ρs[j]), color=:black, linewidth=3) for (i, j) in zip([1, 3, 5, 7], 1:4)]
[scatter!(fig.content[i], real(ξs[j]), imag(ξs[j]), color=:white, markersize=j < 4 ? 8 : 6) for (i, j) in zip([1, 3, 5, 7], 1:4)]
[xlims!(fig.content[i], -6, 6) for i in [1, 3, 5, 7]]
[ylims!(fig.content[i], -6, 6) for i in [1, 3, 5, 7]]

ρ = 2.4:0.01:200.0
θs2 = [large_ξ_roots_θ.(ρ, Float64(μ[i]), 2) for i in 1:4]
ρs2 = [@. ρ * exp(im * θs2[i]) for i in 1:4]
[lines!(fig.content[i], real(ρs2[j]), imag(ρs2[j]), color=:black, linewidth=3) for (i, j) in zip([1, 3, 5, 7], 1:4)]

ξ12 = [sqrt(8n * Float64(μ1) * π) for n in 1:100]
ξ22 = [sqrt(8n * Float64(μ2) * π) for n in 1:100]
ξ32 = [sqrt(8n * Float64(μ3) * π) for n in 1:100]
ξ42 = [complex(-1.6, 1.6), complex(-2.0, 2.0), complex(-2.3, 2.3),
    complex(-2.6, 2.6), complex(-2.85, 2.85), complex(-3.1, 3.1),
    complex(-3.3, 3.30), complex(-3.5, 3.5), complex(-3.65, 3.65),
    [sqrt(8n * Float64(μ4) * π) for n in 12:100]...]
ξ2 = [ξ12, ξ22, ξ32, ξ42]
ξs2 = [large_ξ_roots_ρ.(abs.(ξ2[i]), Float64(μ[i]), 2) for i in 1:4]
ρ2 = 2.4:0.01:200.0
θs2 = [large_ξ_roots_θ.(ρ2, Float64(μ[i]), 2) for i in 1:4]
ρs2 = [@. ρ2 * exp(im * θs2[i]) for i in 1:4]
[lines!(fig.content[i], real(ρs2[j]), imag(ρs2[j]), color=:black, linewidth=3) for (i, j) in zip([1, 3, 5, 7], 1:4)]
[scatter!(fig.content[i], real(ξs2[j]), imag(ξs2[j]), color=:white, markersize=j < 4 ? 8 : 6) for (i, j) in zip([1, 3, 5, 7], 1:4)]
[xlims!(fig.content[i], -6, 6) for i in [1, 3, 5, 7]]
[ylims!(fig.content[i], -6, 6) for i in [1, 3, 5, 7]]

resize_to_layout!(fig)

fig

save("$FIGURES/similarity_solutions.$EXTENSION", fig)