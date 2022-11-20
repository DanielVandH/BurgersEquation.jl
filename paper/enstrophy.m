mu = [0.0, linspace(1e-2, 0.5, 250)];
L = 30;
a = -L; b = L;
d = [a b];
M = 501; N = 250;
t_max = 5; t_min = 0;
enstrophies = zeros(M, length(mu));
dt = (t_max - t_min) / (M-1); t = t_min:dt:t_max;
opts = pdeset('Eps', 1e-7);
bc.left = @(u) diff(u); bc.right = @(u) diff(u);
x = chebfun(@(x) x, d);
f = 1./(1+x.^2);
for k = length(mu):-1:2
    k
    pdefun = @(t, x, u) mu(k)*diff(u, 2) - diff(u.^2 / 2);
    [t, u] = pde15s(pdefun, t, f, bc, opts);
    enstrophies(:, k) = sum(diff(u)^2)*0.5;
end
%%
pdefun = @(t, x, u) - diff(u.^2 / 2);
bcc.left = @(u) diff(u);
[t, u] = pde23t(pdefun, t(1:151), f, bcc, pdeset('Eps', 1e-7));


%%
writematrix(t, 'data/T_ENSTROPHY.dat');
writematrix(mu', 'data/mu_ENSTROPHY.dat');
writematrix(enstrophies, 'data/ENSTROPHY.dat', 'Delimiter', ' ');