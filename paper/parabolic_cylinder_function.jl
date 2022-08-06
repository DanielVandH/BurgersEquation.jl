bt = 2048
x1 = collect(-6.0:0.01:6.0);
x2 = collect(-6.0:0.01:6.0);
x3 = collect(-6.0:0.01:6.0);
x4 = ArbReal.(-6:0.05:6, bits=bt);
y1 = collect(-6.0:0.01:6.0);
y2 = collect(-6.0:0.01:6.0);
y3 = collect(-6.0:0.01:6.0);
y4 = ArbReal.(-6:0.05:6, bits=bt);
μ1 = 2.0;
μ2 = 1.0;
μ3 = 0.5;
μ4 = ArbReal(0.1, bits=bt);
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

## Add the roots now 
#=
ρ = 0.01:0.01:10.0
θ_roots = zeros(ComplexF64, length(ρ), length(μ))
for j in eachindex(μ)
    for i in eachindex(ρ)
        θ_roots[i, j] = large_ξ_roots_θ(ρ[i], Float64(μ[j]), 1)
    end
end

lines!(fig.content[1], real(θ_roots[:, 1]), imag(θ_roots[:, 1]), color=:black, linewidth=6)
lines!(fig.content[3], real(θ_roots[:, 2]), imag(θ_roots[:, 2]), color=:black, linewidth=6)
lines!(fig.content[5], real(θ_roots[:, 3]), imag(θ_roots[:, 3]), color=:black, linewidth=6)
lines!(fig.content[7], real(θ_roots[:, 4]), imag(θ_roots[:, 4]), color=:black, linewidth=6)
=#

fig = Figure()
j = 1
xp = vec(real(z[j]))
yp = vec(imag(z[j]))
zp = vec(angle.(-Complex{Float64}.(Φ₀_vals64[j])))
ax = Axis(fig[1, 1])
scatter!(ax, vec(real(z[j])), vec(imag(z[j])), color=vec(angle.(-Complex{Float64}.(Φ₀_vals64[j]))), colormap=ComplexPortraits.nist_colors(241))
on(events(ax).mousebutton) do event
    if event.button == Mouse.left
        mp = events(ax).mouseposition[]
        @show mp
    end
end
fig

scene = Scene()
scatter!(scene, xp, yp, color=zp, colormap=ComplexPortraits.nist_colors(241))
on(events(scene).mousebutton) do event
    if event.button == Mouse.left
        mp = events(scene).mouseposition[]
        @show mp
    end
end

j = 1
fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, vec(real(z[j])), vec(imag(z[j])), color=vec(angle.(-Complex{Float64}.(Φ₀_vals64[j]))), colormap=ComplexPortraits.nist_colors(241))
fig

ξ1 = [complex(0.0020, 1.0041)]
ξ2 = [complex(1.650e-3, 1.0024), complex(0.00365, 1.0044), complex(0.00500, 1.0057)]
ξ3 = [complex(0.0014, 1.0012), complex(0.00265, 1.00270), complex(0.0036, 1.0037),
    complex(0.00435, 1.0045), complex(0.0051, 1.00525), complex(0.0056, 1.0058)]
ξ4 = [complex(0.0013, 0.9998), complex(0.0015, 1.0005), complex(0.0017, 1.00095),
    complex(0.0021, 1.0012), complex(0.00235, 1.0017), complex(0.00260, 1.0020),
    complex(0.00280, 1.0023), complex(0.0030, 1.0025), complex(0.00320, 1.0027),
    complex(0.00340, 1.0029)]

scatter!(ax, real(ξ2), imag(ξ2))
ξ = large_ξ_roots_ρ.(abs.(ξ4), 0.1)
scatter!(ax, real(ξ), imag(ξ))

save("$FIGURES/similarity_solutions.$EXTENSION", fig)
