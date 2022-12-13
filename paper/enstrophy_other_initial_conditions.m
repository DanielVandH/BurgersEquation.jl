%% 1/(1+x^2)^2
mu = [0.0, linspace(1e-2, 0.5, 250)];
L = 30;
a = -L; b = L;
d = [a b];
M = 5001; N = 250;
t_max = 25; t_min = 0;
enstrophies = zeros(M, length(mu));
dt = (t_max - t_min) / (M-1); t = t_min:dt:t_max;
opts = pdeset('Eps', 1e-7);
bc.left = @(u) diff(u); bc.right = @(u) diff(u);
x = chebfun(@(x) x, d);
f = 1./(1+x.^2).^2;
for k = length(mu):-1:2
    k
    pdefun = @(t, x, u) mu(k)*diff(u, 2) - diff(u.^2 / 2);
    [t, u] = pde15s(pdefun, t, f, bc, opts);
    enstrophies(:, k) = sum(diff(u)^2)*0.5;
end
pdefun = @(t, x, u) - diff(u.^2 / 2);
bcc.left = @(u) diff(u);
[t, u] = pde15s(pdefun, t(t < 0.8), f, bcc, opts);
enstrophies(1:length(t), 1) = sum(diff(u)^2)*0.5;
enstrophies((length(t)+1):end, 1) = NaN;
t = t_min:dt:t_max;
writematrix(t', 'data/T_ENSTROPHY_1O1PX2A2.dat');
writematrix(mu', 'data/mu_ENSTROPHY_1O1PX2A2.dat');
writematrix(enstrophies, 'data/ENSTROPHY_1O1PX2A2.dat', 'Delimiter', ' ');

%% 1/(1+x^2)^(1/2)
mu = [0.0, linspace(1e-2, 0.5, 250)];
L = 30;
a = -L; b = L;
d = [a b];
M = 5001; N = 250;
t_max = 25; t_min = 0;
enstrophies = zeros(M, length(mu));
dt = (t_max - t_min) / (M-1); t = t_min:dt:t_max;
opts = pdeset('Eps', 1e-7);
bc.left = @(u) diff(u); bc.right = @(u) diff(u);
x = chebfun(@(x) x, d);
f = 1./(1+x.^2).^(1/2);
for k = length(mu):-1:2
    k
    pdefun = @(t, x, u) mu(k)*diff(u, 2) - diff(u.^2 / 2);
    [t, u] = pde15s(pdefun, t, f, bc, opts);
    enstrophies(:, k) = sum(diff(u)^2)*0.5;
end
pdefun = @(t, x, u) - diff(u.^2 / 2);
bcc.left = @(u) diff(u);
[t, u] = pde15s(pdefun, t(t < 1.2), f, bcc, opts);
enstrophies(1:length(t), 1) = sum(diff(u)^2)*0.5;
enstrophies((length(t)+1):end, 1) = NaN;
t = t_min:dt:t_max;
writematrix(t', 'data/T_ENSTROPHY_1OSQRT1PX2.dat');
writematrix(mu', 'data/mu_ENSTROPHY_1OSQRT1PX2.dat');
writematrix(enstrophies, 'data/ENSTROPHY_1OSQRT1PX2.dat', 'Delimiter', ' ');