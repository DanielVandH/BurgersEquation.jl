%% Start by finding all the solutions
mu = [0.05, 0.1468, 0.5];
L = 15; scaled_L = 5;
a = -L; b = L;
d = [a b];
M = 501; N = 250;
u_vals = zeros(N, M, length(mu));
t_max = 2; t_min = 0;
dt = (t_max - t_min) / (M-1); t = t_min:dt:t_max;
opts = pdeset('Eps', 1e-11);
cheb_x = scaled_L * chebpts(N);
bc.left = @(u) diff(u); bc.right = @(u) diff(u);
x = chebfun(@(x) x, d);
f = 1./(1+x.^2);
for k = 1:3
    pdefun = @(t, x, u) mu(k)*diff(u, 2) - diff(u.^2 / 2);
    [~, u] = pde15s(pdefun, t, f, bc, opts);
    for j = 1:M
        u_vals(:, j, k) = u(cheb_x, j);
    end
end

%% Now loop over and find the poles
aaa_poles = zeros(M, 3);
z0 = [2.0714+0.4821i, 2.1270+1.225i, 2.3123+2.7113i];
aaa_poles(end, :) = z0;
for k = 1:3
    for j = (M-1):-1:1
        dat = u_vals(:, j, k);
        [~, poles, res] = aaa(dat, cheb_x);
        poles = poles(abs(res) > 1e-4);
        [~, idx] = min(abs(poles - aaa_poles(j+1, k)));
        aaa_poles(j, k) = poles(idx);
    end
end
writematrix(t, 'data/T_AAA_POLES.dat', 'Delimiter', ',');
writematrix(aaa_poles, 'data/AAA_POLES.dat', 'Delimiter', ' ');