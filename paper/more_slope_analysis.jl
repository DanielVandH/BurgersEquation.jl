## 1/(1+x^2)
μ = [0, 0.01, 0.02, 0.03, 0.04, 0.05]
tₛ₁ = 8sqrt(3) / 9
t = [0.1, 0.5, 1.0, tₛ₁]
x = collect(LinRange(-10, 10, 2000))
f = y -> 1 / (1 + y^2)
f′ = y -> -2y / (1 + y^2)^2

# u(x, t)
u_inviscid = inviscid_solution(f, f′, x, t)
u_viscous = viscous_solution(x, t, μ[2:end]; ic=1)

fig = Figure(fontsize=38, resolution=(2239.2961f0, 1120.5579f0))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"$u(x,t)$",
    title=L"(a): $\mu = 0$, $u(x, 0) = 1/(1+x^2)$", titlealign=:left,
    width=600, height=400)
series!(ax, x, u_inviscid', linewidth=4, color=[:red, :black, :blue, :darkgreen])

for i in axes(u_viscous, 3)
    u_viscous_μ = @views u_viscous[:, :, i]
    ax = Axis(fig[i > 2 ? 2 : 1, (i%3)+1],
        xlabel=L"x", ylabel=L"$u(x,t)$",
        title=L"(%$(ALPHABET[i+1])): $\mu = %$(μ[i+1])$, $u(x, 0) = 1/(1+x^2)$", titlealign=:left,
        width=600, height=400)
    series!(ax, x, u_viscous_μ', linewidth=4, color=[:red, :black, :blue, :darkgreen])
    ylims!(ax, 0, 1)
end

fig_runge = fig

# uₓ(x, t)
uₓ_inviscid = similar(u_inviscid)
for j in eachindex(t)
    for i in eachindex(x)
        uₓ_inviscid[i, j] = f′(x[i] - u_inviscid[i, j] * t[j]) / (1 + t[j] * f′(x[i] - u_inviscid[i, j] * t[j]))
    end
end
num_nodes = 250
nodes, weights = gausshermite(num_nodes)
∇u = (x, t, μ) -> ForwardDiff.derivative(x -> viscous_solution(x, t, μ, nodes, weights; ic=1), x)
uₓ_viscous = similar(u_viscous)
for (i, X) in enumerate(x)
    for (j, T) in enumerate(t)
        for (k, M) in enumerate(μ[2:end])
            uₓ_viscous[i, j, k] = ∇u(X, T, M)
        end
    end
end

fig = Figure(fontsize=38, resolution=(2239.2961f0, 1120.5579f0))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"$u_x(x, t)$",
    title=L"(a): $\mu = 0$, $u(x, 0) = 1/(1+x^2)$", titlealign=:left,
    width=600, height=400)
series!(ax, x, abs.(uₓ_inviscid'), linewidth=4, color=[:red, :black, :blue, :darkgreen])
ylims!(ax, 0, 4)

for i in axes(u_viscous, 3)
    uₓ_viscous_μ = @views uₓ_viscous[:, :, i]
    ax = Axis(fig[i > 2 ? 2 : 1, (i%3)+1],
        xlabel=L"x", ylabel=L"$u_x(x, t)$",
        title=L"(%$(ALPHABET[i+1])): $\mu = %$(μ[i+1])$, $u(x, 0) = 1/(1+x^2)$", titlealign=:left,
        width=600, height=400)
    series!(ax, x, abs.(uₓ_viscous_μ'), linewidth=4, color=[:red, :black, :blue, :darkgreen])
    ylims!(ax, 0, 4)
end

fig_runge_du = fig

# Alternative views 
fig = Figure(fontsize=38, resolution=(1461.5542f0, 1120.5579f0))
time_labels = ["0.1", "0.5", "1", "8/27^{1/2}"]
for j in eachindex(t)
    ax = Axis(fig[j > 2 ? 2 : 1, mod1(j, 2)],
        xlabel=L"x", ylabel=L"$u(x, t)",
        title=L"(%$(ALPHABET[j])): $t = %$(time_labels[j])$, $u(x, 0) = 1/(1+x^2)$", titlealign=:left,
        width=600, height=400)
    lines!(ax, x, u_inviscid[:, j], color=:red)
    lines!(ax, x, u_viscous[:, j, 1], color=:black)
    lines!(ax, x, u_viscous[:, j, 2], color=:blue)
    lines!(ax, x, u_viscous[:, j, 3], color=:darkgreen)
    lines!(ax, x, u_viscous[:, j, 4], color=:magenta)
    lines!(ax, x, u_viscous[:, j, 5], color=:orange)
    ylims!(ax, 0.8, 1.1)
end

fig_runge_t = fig

fig = Figure(fontsize=38, resolution=(1461.5542f0, 1120.5579f0))
time_labels = ["0.1", "0.5", "1", "8/27^{1/2}"]
for j in eachindex(t)
    ax = Axis(fig[j > 2 ? 2 : 1, mod1(j, 2)],
        xlabel=L"x", ylabel=L"$u_x(x, t)",
        title=L"(%$(ALPHABET[j])): $t = %$(time_labels[j])$, $u(x, 0) = 1/(1+x^2)$", titlealign=:left,
        width=600, height=400)
    lines!(ax, x, abs.(uₓ_inviscid[:, j]), color=:red)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 1]), color=:black)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 2]), color=:blue)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 3]), color=:darkgreen)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 4]), color=:magenta)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 5]), color=:orange)
    ylims!(ax, 0, 2)
end

fig_runge_t_du = fig

## 1/(1+x^2)^2
μ = [0, 0.01, 0.02, 0.03, 0.04, 0.05]
tₛ₂ = 54sqrt(5) / 125
t = [0.1, 0.5, 0.75, tₛ₂]
x = collect(LinRange(-10, 10, 2000))
f = z -> 1 / (1 + z^2)^2
f′ = z -> -4z / (1 + z^2)^3

# u(x, t)
u_inviscid = inviscid_solution(f, f′, x, t)
u_viscous = viscous_solution(x, t, μ[2:end]; ic=2)

fig = Figure(fontsize=38, resolution=(2239.2961f0, 1120.5579f0))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"$u(x,t)$",
    title=L"(a): $\mu = 0$, $u(x, 0) = 1/(1+x^2)^2$", titlealign=:left,
    width=600, height=400)
series!(ax, x, u_inviscid', linewidth=4, color=[:red, :black, :blue, :darkgreen])

for i in axes(u_viscous, 3)
    u_viscous_μ = @views u_viscous[:, :, i]
    ax = Axis(fig[i > 2 ? 2 : 1, (i%3)+1],
        xlabel=L"x", ylabel=L"$u(x,t)$",
        title=L"(%$(ALPHABET[i+1])): $\mu = %$(μ[i+1])$, $u(x, 0) = 1/(1+x^2)^2$", titlealign=:left,
        width=600, height=400)
    _ax[i > 2 ? 2 : 1, (i%3)+1] = ax
    series!(ax, x, u_viscous_μ', linewidth=4, color=[:red, :black, :blue, :darkgreen])
end

fig_square = fig

# uₓ(x, t)
uₓ_inviscid = similar(u_inviscid)
for j in eachindex(t)
    for i in eachindex(x)
        uₓ_inviscid[i, j] = f′(x[i] - u_inviscid[i, j] * t[j]) / (1 + t[j] * f′(x[i] - u_inviscid[i, j] * t[j]))
    end
end
num_nodes = 250
nodes, weights = gausshermite(num_nodes)
∇u = (x, t, μ) -> ForwardDiff.derivative(x -> viscous_solution(x, t, μ, nodes, weights; ic=2), x)
uₓ_viscous = similar(u_viscous)
for (i, X) in enumerate(x)
    for (j, T) in enumerate(t)
        for (k, M) in enumerate(μ[2:end])
            uₓ_viscous[i, j, k] = ∇u(X, T, M)
        end
    end
end

fig = Figure(fontsize=38, resolution=(2239.2961f0, 1120.5579f0))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"$u_x(x, t)$",
    title=L"(a): $\mu = 0$, $u(x, 0) = 1/(1+x^2)^2$", titlealign=:left,
    width=600, height=400)
series!(ax, x, abs.(uₓ_inviscid'), linewidth=4, color=[:red, :black, :blue, :darkgreen])
ylims!(ax, 0, 4)

for i in axes(u_viscous, 3)
    uₓ_viscous_μ = @views uₓ_viscous[:, :, i]
    ax = Axis(fig[i > 2 ? 2 : 1, (i%3)+1],
        xlabel=L"x", ylabel=L"$u_x(x, t)$",
        title=L"(%$(ALPHABET[i+1])): $\mu = %$(μ[i+1])$, $u(x, 0) = 1/(1+x^2)^2$", titlealign=:left,
        width=600, height=400)
    _ax[i > 2 ? 2 : 1, (i%3)+1] = ax
    series!(ax, x, abs.(uₓ_viscous_μ'), linewidth=4, color=[:red, :black, :blue, :darkgreen])
    ylims!(ax, 0, 4)
end

fig_square_du = fig

# Alternative views 
fig = Figure(fontsize=38, resolution=(1461.5542f0, 1120.5579f0))
time_labels = ["0.1", "0.5", "0.75", "54/5^{5/2}"]
for j in eachindex(t)
    ax = Axis(fig[j > 2 ? 2 : 1, mod1(j, 2)],
        xlabel=L"x", ylabel=L"$u(x, t)",
        title=L"(%$(ALPHABET[j])): $t = %$(time_labels[j])$, $u(x, 0) = 1/(1+x^2)^2$", titlealign=:left,
        width=600, height=400)
    lines!(ax, x, u_inviscid[:, j], color=:red)
    lines!(ax, x, u_viscous[:, j, 1], color=:black)
    lines!(ax, x, u_viscous[:, j, 2], color=:blue)
    lines!(ax, x, u_viscous[:, j, 3], color=:darkgreen)
    lines!(ax, x, u_viscous[:, j, 4], color=:magenta)
    lines!(ax, x, u_viscous[:, j, 5], color=:orange)
    ylims!(ax, 0.8, 1.1)
end

fig_square_t = fig

fig = Figure(fontsize=38, resolution=(1461.5542f0, 1120.5579f0))
time_labels = ["0.1", "0.5", "0.75", "54/5^{5/2}"]
for j in eachindex(t)
    ax = Axis(fig[j > 2 ? 2 : 1, mod1(j, 2)],
        xlabel=L"x", ylabel=L"$u_x(x, t)",
        title=L"(%$(ALPHABET[j])): $t = %$(time_labels[j])$, $u(x, 0) = 1/(1+x^2)^2$", titlealign=:left,
        width=600, height=400)
    lines!(ax, x, abs.(uₓ_inviscid[:, j]), color=:red)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 1]), color=:black)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 2]), color=:blue)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 3]), color=:darkgreen)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 4]), color=:magenta)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 5]), color=:orange)
    ylims!(ax, 0, 4)
end

fig_square_t_du = fig

## 1/(1+x^2)^(1/2)
μ = [0, 0.01, 0.02, 0.03, 0.04, 0.05]
tₛ₃ = 3sqrt(3) / 2
t = [0.1, 0.5, 1.0, tₛ₃]
x = collect(LinRange(-10, 10, 2000))
f = z -> 1 / (1 + z^2)^(1 / 2)
f′ = z -> -z / (1 + z^2)^(3 / 2)

# u(x, t)
u_inviscid = inviscid_solution(f, f′, x, t)
u_viscous = viscous_solution(x, t, μ[2:end]; ic=3)

fig = Figure(fontsize=38, resolution=(2239.2961f0, 1120.5579f0))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"$u(x,t)$",
    title=L"(a): $\mu = 0$, $u(x, 0) = 1/(1+x^2)^{1/2}$", titlealign=:left,
    width=600, height=400)
series!(ax, x, u_inviscid', linewidth=4, color=[:red, :black, :blue, :darkgreen])

for i in axes(u_viscous, 3)
    u_viscous_μ = @views u_viscous[:, :, i]
    ax = Axis(fig[i > 2 ? 2 : 1, (i%3)+1],
        xlabel=L"x", ylabel=L"$u(x,t)$",
        title=L"(%$(ALPHABET[i+1])): $\mu = %$(μ[i+1])$, $u(x, 0) = 1/(1+x^2)^{1/2}$", titlealign=:left,
        width=600, height=400)
    series!(ax, x, u_viscous_μ', linewidth=4, color=[:red, :black, :blue, :darkgreen])
end

fig_sqrt = fig

# uₓ(x, t)
uₓ_inviscid = similar(u_inviscid)
for j in eachindex(t)
    for i in eachindex(x)
        uₓ_inviscid[i, j] = f′(x[i] - u_inviscid[i, j] * t[j]) / (1 + t[j] * f′(x[i] - u_inviscid[i, j] * t[j]))
    end
end
num_nodes = 250
nodes, weights = gausshermite(num_nodes)
∇u = (x, t, μ) -> ForwardDiff.derivative(x -> viscous_solution(x, t, μ, nodes, weights; ic=3), x)
uₓ_viscous = similar(u_viscous)
for (i, X) in enumerate(x)
    for (j, T) in enumerate(t)
        for (k, M) in enumerate(μ[2:end])
            uₓ_viscous[i, j, k] = ∇u(X, T, M)
        end
    end
end

fig = Figure(fontsize=38, resolution=(2239.2961f0, 1120.5579f0))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"$u_x(x, t)$",
    title=L"(a): $\mu = 0$, $u(x, 0) = 1/(1+x^2)^{1/2}$", titlealign=:left,
    width=600, height=400)
series!(ax, x, abs.(uₓ_inviscid'), linewidth=4, color=[:red, :black, :blue, :darkgreen])
ylims!(ax, 0, 4)

for i in axes(u_viscous, 3)
    uₓ_viscous_μ = @views uₓ_viscous[:, :, i]
    ax = Axis(fig[i > 2 ? 2 : 1, (i%3)+1],
        xlabel=L"x", ylabel=L"$u_x(x, t)$",
        title=L"(%$(ALPHABET[i+1])): $\mu = %$(μ[i+1])$, $u(x, 0) = 1/(1+x^2)^{1/2}$", titlealign=:left,
        width=600, height=400)
    series!(ax, x, abs.(uₓ_viscous_μ'), linewidth=4, color=[:red, :black, :blue, :darkgreen])
    ylims!(ax, 0, 4)
end

fig_sqrt_du = fig

# Alternative views 
fig = Figure(fontsize=38, resolution=(1461.5542f0, 1120.5579f0))
time_labels = ["0.1", "0.5", "1", "27^{1/2}/2"]
for j in eachindex(t)
    ax = Axis(fig[j > 2 ? 2 : 1, mod1(j, 2)],
        xlabel=L"x", ylabel=L"$u(x, t)",
        title=L"(%$(ALPHABET[j])): $t = %$(time_labels[j])$, $u(x, 0) = 1/(1+x^2)^{1/2}$", titlealign=:left,
        width=600, height=400)
    lines!(ax, x, u_inviscid[:, j], color=:red)
    lines!(ax, x, u_viscous[:, j, 1], color=:black)
    lines!(ax, x, u_viscous[:, j, 2], color=:blue)
    lines!(ax, x, u_viscous[:, j, 3], color=:darkgreen)
    lines!(ax, x, u_viscous[:, j, 4], color=:magenta)
    lines!(ax, x, u_viscous[:, j, 5], color=:orange)
    ylims!(ax, 0.8, 1.1)
end

fig_sqrt_t = fig

fig = Figure(fontsize=38, resolution=(1461.5542f0, 1120.5579f0))
time_labels = ["0.1", "0.5", "1", "27^{1/2}/2"]
for j in eachindex(t)
    ax = Axis(fig[j > 2 ? 2 : 1, mod1(j, 2)],
        xlabel=L"x", ylabel=L"$u_x(x, t)",
        title=L"(%$(ALPHABET[j])): $t = %$(time_labels[j])$, $u(x, 0) = 1/(1+x^2)^{1/2}$", titlealign=:left,
        width=600, height=400)
    lines!(ax, x, abs.(uₓ_inviscid[:, j]), color=:red)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 1]), color=:black)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 2]), color=:blue)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 3]), color=:darkgreen)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 4]), color=:magenta)
    lines!(ax, x, abs.(uₓ_viscous[:, j, 5]), color=:orange)
    ylims!(ax, 0, 2)
end
fig_sqrt_t_du = fig
