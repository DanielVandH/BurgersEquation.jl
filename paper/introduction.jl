x = -12:0.1:12
t = [0.0, tₛ, 5.0, 10.0]
μ = [0.1, 1.0]
t₀ = 25.26

u_inviscid = inviscid_solution(x -> 1 / (1 + x^2), x, tₛ)
u_viscous = viscous_solution(x, t, μ)

fig = Figure(fontsize=38, resolution=(800, 400))
ab = ['a', 'b']
for k in 1:length(μ)
    ax = Axis(fig[1, k], xlabel=L"x", ylabel=L"u(x,t)",
        title=L"(%$(ab[k])): $\mu = %$(μ[k])$", titlealign=:left,
        width=600, height=400,
        xticks=([-20, -10, 0, 10, 20], [L"-20", L"-10", L"0", L"10", L"20"]),
        yticks=([0, 0.25, 0.50, 0.75, 1.00], [L"0", L"0.25", L"0.5", L"0.75", L"1"]))
    k == 1 && lines!(ax, x, u_inviscid, color=:red, linestyle=:dash, linewidth=2.3)
    for j in 1:length(t)
        lines!(ax, x, u_viscous[:, j, k], color=:black, linestyle=:solid, linewidth=2.3)
    end
end
resize_to_layout!(fig)

save("$FIGURES/inviscid_compared_to_viscous_solution.$EXTENSION", fig)

