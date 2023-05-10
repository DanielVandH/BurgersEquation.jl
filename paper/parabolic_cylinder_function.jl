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
    complex(3.0, 2.5), complex(3.2, 2.7), complex(3.4, 3.0),
    complex(3.6, 3.1), complex(3.7, 3.3), complex(3.9, 3.5),
    complex(4.025, 3.65), [sqrt(8n * Float64(μ4) * π) for n in 13:100]...]
ξ = [ξ1, ξ2, ξ3, ξ4]
ξs = [large_ξ_roots_ρ.(abs.(ξ[i]), Float64(μ[i])) for i in 1:4]
ρ = 1.0:0.01:200.0
θs = [large_ξ_roots_θ.(ρ, Float64(μ[i]), 1) for i in 1:4]
ρs = [@. ρ * exp(im * θs[i]) for i in 1:4]
ξs[4][1] = complex(NaN, NaN)
ξs[4][7] = 2.84 + 2.2825im
ξs[4][12] = 3.728 + 3.346im
push!(ξs[4], 5.230 + 4.935im)
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

#if RUN_EXTRAS
##### Just need to get some data 
## Getting some more results 
# First, for 0.1, 0.5, 1, 2 
μ1 = 2.0
μ2 = 1.0
μ3 = 0.5
μ4 = 0.1
μ = [μ1, μ2, μ3, μ4]
ξ1 = [complex(1.9, 4.1), [sqrt(8n * Float64(μ1) * π) for n in 2:100]...]
ξ2 = [complex(1.75, 2.5), complex(3.75, 4.4), complex(5.0, 5.6),
    [sqrt(8n * Float64(μ2) * π) for n in 4:100]...]
ξ3 = [complex(1.5, 1.2), complex(2.7, 2.8), complex(3.6, 3.75),
    complex(4.4, 4.5), complex(5.0, 5.2), complex(5.6, 5.7),
    [sqrt(8n * Float64(μ3) * π) for n in 7:100]...]
ξ4 = [complex(1.2, -0.1), complex(1.5, 0.5), complex(1.8, 1.0), complex(2.2, 1.4),
    complex(2.4, 1.8),
    complex(2.6, 2.0), complex(2.8, 2.2),
    complex(3.0, 2.5), complex(3.2, 2.7), complex(3.4, 3.0),
    complex(3.6, 3.1), complex(3.7, 3.3), complex(3.9, 3.5),
    complex(4.025, 3.65), [sqrt(8n * Float64(μ4) * π) for n in 13:100]...]
ξ = [ξ1, ξ2, ξ3, ξ4]
ξs = [large_ξ_roots_ρ.(abs.(ξ[i]), Float64(μ[i])) for i in 1:4]
ξs[4][1] = complex(NaN, NaN)
ξs[4][7] = 2.84 + 2.2825im
ξs[4][12] = 3.728 + 3.346im
push!(ξs[4], 5.230 + 4.935im)

ξ12 = [sqrt(8n * Float64(μ1) * π) for n in 1:100]
ξ22 = [sqrt(8n * Float64(μ2) * π) for n in 1:100]
ξ32 = [sqrt(8n * Float64(μ3) * π) for n in 1:100]
ξ42 = [complex(-1.6, 1.6), complex(-2.0, 2.0), complex(-2.3, 2.3),
    complex(-2.6, 2.6), complex(-2.85, 2.85), complex(-3.1, 3.1),
    complex(-3.3, 3.30), complex(-3.5, 3.5), complex(-3.65, 3.65),
    [sqrt(8n * Float64(μ4) * π) for n in 12:100]...]
ξ2 = [ξ12, ξ22, ξ32, ξ42]
ξs2 = [large_ξ_roots_ρ.(abs.(ξ2[i]), Float64(μ[i]), 2) for i in 1:4]

# Now consider μ = 0.25, 0.75, 1.5
bt = 2048
μ5 = ArbReal(0.25, bits=bt)
μ6 = 0.75
μ7 = 1.5
x5 = ArbReal.(-6.0:0.05:6.0, bits=bt)
x6 = collect(-6.0:0.01:6.0)
x7 = collect(-6.0:0.01:6.0)
y5 = ArbReal.(-6.0:0.05:6.0, bits=bt)
y6 = collect(-6.0:0.01:6.0)
y7 = collect(-6.0:0.01:6.0)
Φ₀_valsc = [Φ₀(x5, y5, μ5), Φ₀(x6, y6, μ6), Φ₀(x7, y7, μ7)]
tc = 1e-6
zc = Vector{Matrix{Union{ComplexF64,ArbComplex{P}} where {P}}}(undef, 3)
uc = Vector{Matrix{Union{ComplexF64,ArbComplex{P}} where {P}}}(undef, 3)
num_nodes = 250
nodes, weights = gausshermite(num_nodes)
glnodes, glweights = gausslegendre(num_nodes)
xc = [x5, x6, x7]
yc = [y5, y6, y7]
μc = [μ5, μ6, μ7]
for j in 1:3
    xx = Float64.(xc[j])
    yy = Float64.(yc[j])
    zc[j] = [im .+ sqrt(tc) * (xx + im * yy) for xx in xx, yy in yy]
    uc[j] = viscous_solution.(real(zc[j]), imag(zc[j]), tc, Float64(μc[j]), Ref(nodes), Ref(weights), Ref(glnodes), Ref(glweights))
end
fig = Figure(fontsize=38, resolution=(2000, 1600))
plot_layoutc = [[(1, 1), (2, 1)], [(1, 2), (2, 2)], [(1, 3), (2, 3)]]
for j = 1:3
    portrait!(fig, xc[j], yc[j], Complex{Float64}.(Φ₀_valsc[j]), plot_layoutc[j][1][1], plot_layoutc[j][1][2]; nist=NIST,
        title=L"(%$(ALPHABET[j])): $\mu = %$(Float64(μc[j]))$", titlealign=:left,
        xlabel=L"\mathrm{Re}(\xi)", ylabel=L"\mathrm{Im}(\xi)", width=450, height=450,
        xticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]),
        yticks=([-5.0, 0.0, 5.0], [L"-5", L"0", L"5"]))

    portrait!(fig, real.(zc[j]), imag.(zc[j]), Complex{Float64}.(uc[j]), plot_layoutc[j][2][1], plot_layoutc[j][2][2]; nist=NIST,
        title=L"(%$(ALPHABET[3+j])): $\mu = %$(Float64(μc[j]))$, $t = 10^{-6}$", titlealign=:left,
        xlabel=L"x", ylabel=L"y", width=450, height=450,
        xticks=([-0.005, 0.0, 0.005], [L"-0.005", L"0", L"0.005"]),
        yticks=([0.995, 1.0, 1.005], [L"0.995", L"1", L"1.005"]))
end

ξ5 = [complex(1.29770, 0.40120), complex(1.775, 1.5250), [sqrt(8n * Float64(μ5) * π) for n in 2:100]...]
ξ52 = [complex(-1.8765, 2.29575), complex(-2.6525, 2.9635), [sqrt(8n * Float64(μ5) * π) for n in 4:100]...]
ξ6 = [complex(1.52595, 1.78915), complex(3.1749, 3.58125), [sqrt(8n * Float64(μ6) * π) for n in 2:100]...]
ξ62 = [sqrt(8n * Float64(μ6) * π) for n in 1:100]
ξ7 = [complex(1.7593, 3.2794), complex(4.3502, 5.5392), [sqrt(8n * Float64(μ7) * π) for n in 3:100]...]
ξ72 = [complex(-3.44, 5.16), [sqrt(8n * Float64(μ7) * π) for n in 2:100]...]

ξsc = [large_ξ_roots_ρ.(abs.([ξ5, ξ6, ξ7][i]), Float64.([μ5, μ6, μ7][i])) for i in 1:3]
ξs2c = [large_ξ_roots_ρ.(abs.([ξ52, ξ62, ξ72][i]), Float64.([μ5, μ6, μ7][i]), 2) for i in 1:3]
ξsc[3][1] = complex(1.7593, 3.2794)
ξsc[3][2] = complex(4.3502, 5.5392)

tgt = 105
for _ξ in (ξs, ξs2, ξsc, ξs2c)
    [append!(_ξ[i], repeat([complex(NaN, NaN)], tgt - length(_ξ[i]))) for i in eachindex(_ξ)]
end

idx = [sortperm(abs.(ξs[i])) for i in eachindex(ξs)]
idx2 = [sortperm(abs.(ξs2[i])) for i in eachindex(ξs2)]
idxc = [sortperm(abs.(ξsc[i])) for i in eachindex(ξsc)]
idx2c = [sortperm(abs.(ξs2c[i])) for i in eachindex(ξs2c)]
[permute!(ξs[i], idx[i]) for i in eachindex(ξs)]
[permute!(ξs2[i], idx2[i]) for i in eachindex(ξs2)]
[permute!(ξsc[i], idxc[i]) for i in eachindex(ξsc)]
[permute!(ξs2c[i], idx2c[i]) for i in eachindex(ξs2c)]

[scatter!(fig.content[i], real(ξsc[j]), imag(ξsc[j]), color=:black, markersize=j < 3 ? 9 : 9) for (i, j) in zip([1, 3, 5], 1:3)]
[scatter!(fig.content[i], real(ξs2c[j]), imag(ξs2c[j]), color=:black, markersize=j < 3 ? 9 : 9) for (i, j) in zip([1, 3, 5], 1:3)]

[xlims!(fig.content[i], -9, 9) for i in [1, 3, 5]]
[ylims!(fig.content[i], -9, 9) for i in [1, 3, 5]]

A1 = hcat(reduce(hcat, ξs), reduce(hcat, ξsc))
A2 = hcat(reduce(hcat, ξs2), reduce(hcat, ξs2c))

D1 = DataFrame(A1, string.(Float64.([μ1, μ2, μ3, μ4, μ5, μ6, μ7])))
D2 = DataFrame(A2, string.(Float64.([μ1, μ2, μ3, μ4, μ5, μ6, μ7])))
CSV.write("$EXTRAS/first_quadrant.csv", D1)
CSV.write("$EXTRAS/second_quadrant.csv", D2)
#end