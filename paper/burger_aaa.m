%% Setup
mu = 0.1;
L = 15; scaled_L = 5;
a = -L; b = L;
d = [a b];
t = [0.1, 0.5, 1.0, 8*sqrt(3)/9];
N = 250;
u_vals = zeros(N, length(t));
bc.left = @(u) diff(u); bc.right = @(u) diff(u);
x = chebfun(@(x) x, d);
f = 1./(1+x.^2);
pdefun = @(t, x, u) mu*diff(u, 2) - diff(u.^2 / 2);
cheb_x = scaled_L * chebpts(N);

%% Find the numerical solution 
[~, u] = pde15s(pdefun, t, f, bc);
for j = 1:length(t)
    u_vals(:, j) = u(cheb_x, j);
end

%% Compute the AAA approximants
r1 = aaa(u(:, 1), cheb_x);
r2 = aaa(u(:, 2), cheb_x);
r3 = aaa(u(:, 3), cheb_x);
r4 = aaa(u(:, 4), cheb_x);
r = {r1, r2, r3, r4};

%% Now obtain the numerical data 
X = -4:0.01:4;
Y = 0:0.01:4;
U = zeros(length(X), length(Y), length(t));
for k = 1:length(t)
    for j = 1:length(Y)
        for i = 1:length(X)
            U(i, j, k) = r{k}(X(i) + 1i*Y(j));
        end
    end
end
writematrix(mu, 'data/mu_AAA.dat')
writematrix(X, 'data/X_AAA.dat', 'Delimiter', ',');
writematrix(Y, 'data/Y_AAA.dat', 'Delimiter', ',');
writematrix(t, 'data/T_AAA.dat', 'Delimiter', ',');
writematrix(U(:, :, 1), 'data/U_AAA_1.dat', 'Delimiter', ' ')
writematrix(U(:, :, 2), 'data/U_AAA_2.dat', 'Delimiter', ' ')
writematrix(U(:, :, 3), 'data/U_AAA_3.dat', 'Delimiter', ' ')
writematrix(U(:, :, 4), 'data/U_AAA_4.dat', 'Delimiter', ' ')
fix_i('U_AAA_1.dat');
fix_i('U_AAA_2.dat');
fix_i('U_AAA_3.dat');
fix_i('U_AAA_4.dat');