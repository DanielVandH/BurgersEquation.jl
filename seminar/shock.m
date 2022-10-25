%% t = 0
t0 = 0;
x_t0_soln = linspace(-12, 12, 2000);
u_t0_soln = 1./(1.+x_t0_soln.^2);

%% t = ts 
ts = 8/sqrt(27);
x0_ts_soln = linspace(-12, 12, 2000);
x_ts_soln = ts./(1+x0_ts_soln.^2) + x0_ts_soln;
u_ts_soln = 1./(1+x0_ts_soln.^2);

%% t = 5
t5 = 5;
x0_t5_soln = linspace(-12, 12, 2000);
x_t5_soln = t5./(1+x0_t5_soln.^2) + x0_t5_soln;
u_t5_soln = 1./(1+x0_t5_soln.^2);

%% t = 10
t10 = 10;
x0_t10_soln = linspace(-12, 12, 2000);
x_t10_soln = t10./(1+x0_t10_soln.^2) + x0_t10_soln;
u_t10_soln = 1./(1+x0_t10_soln.^2);

%% t = 5 (Shock)
t5 = 5;
[x_t5_soln_shock, t5_soln_shock] = ShockSolver(t5, 1e-4);
x_t5_soln_1_shock = linspace(-12, real(x_t5_soln_shock(end)), 2000);
u_t5_soln_1_shock = zeros(1, length(x_t5_soln_1_shock));
i = 1;
for x = x_t5_soln_1_shock 
    [~, uu] = CubicSolv(x, t5);
    u_t5_soln_1_shock(i) = uu;
    i = i + 1;
end
x_t5_soln_2_shock = linspace(real(x_t5_soln_shock(end)), 12, 2000);
u_t5_soln_2_shock = zeros(1, length(x_t5_soln_2_shock));
i = 1;
for x = x_t5_soln_2_shock 
    [uu, ~] = CubicSolv(x, t5);
    u_t5_soln_2_shock(i) = uu;
    i = i + 1;
end
x_t5_soln_shock = [x_t5_soln_1_shock x_t5_soln_2_shock];
u_t5_soln_shock = [u_t5_soln_1_shock u_t5_soln_2_shock];

%% t = 10 (Shock)
t10 = 10;
[x_t10_soln_shock, t10_soln_shock] = ShockSolver(t10, 1e-4);
x_t10_soln_1_shock = linspace(-12, real(x_t10_soln_shock(end)), 2000);
u_t10_soln_1_shock = zeros(1, length(x_t10_soln_1_shock));
i = 1;
for x = x_t10_soln_1_shock 
    [~, uu] = CubicSolv(x, t10);
    u_t10_soln_1_shock(i) = uu;
    i = i + 1;
end
x_t10_soln_2_shock = linspace(real(x_t10_soln_shock(end)), 12, 2000);
u_t10_soln_2_shock = zeros(1, length(x_t10_soln_2_shock));
i = 1;
for x = x_t10_soln_2_shock 
    [uu, ~] = CubicSolv(x, t10);
    u_t10_soln_2_shock(i) = uu;
    i = i + 1;
end
x_t10_soln_shock = [x_t10_soln_1_shock x_t10_soln_2_shock];
u_t10_soln_shock = [u_t10_soln_1_shock u_t10_soln_2_shock];

%% Export
writematrix(x_t0_soln', 'data/xt0.dat');
writematrix(x_ts_soln', 'data/xts.dat');
writematrix(x_t5_soln', 'data/xt5.dat');
writematrix(x_t10_soln', 'data/xt10.dat');
writematrix(x_t5_soln_shock', 'data/xt5_shock.dat');
writematrix(x_t10_soln_shock', 'data/xt10_shock.dat');
writematrix(u_t0_soln', 'data/ut0.dat');
writematrix(u_ts_soln', 'data/uts.dat');
writematrix(u_t5_soln', 'data/ut5.dat');
writematrix(u_t10_soln', 'data/ut10.dat');
writematrix(u_t5_soln_shock', 'data/ut5_shock.dat');
writematrix(u_t10_soln_shock', 'data/ut10_shock.dat');

function [x, t] = ShockSolver(T, h)
% Solves dx/dt = (u^+ + u^-)/2 subject to x(8/sqrt(27)) = sqrt(3), where
% u^+ = min(x1, x2, x3) and u^- = max(x1, x2, x3), where x1, x2, x3 are the roots
% to the cubic t^2u^2 - 2xtu^2 + (1 + x^2)u - 1 = 0. The solution is from t
% = 8/sqrt(27) up to t = T. h is the stepsize for RK4.
ts = 8/sqrt(27);
n = (T-ts)/h;
x = zeros(1,round(n)+1);
t = x;
x(1) = sqrt(3);
t(1) = ts;
for i = 1:length(x)
    [min_r, max_r] = CubicSolv(x(i), t(i));
    k1 = (min_r + max_r)/2;
    [min_r, max_r] = CubicSolv(x(i)+h*k1/2, t(i)+h/2);
    k2 = (min_r + max_r)/2;
    [min_r, max_r] = CubicSolv(x(i) + h * k2/2, t(i) + h/2);
    k3 = (min_r + max_r)/2;
    [min_r, max_r] = CubicSolv(x(i) + h * k3, t(i) + h);
    k4 = (min_r + max_r)/2;
    x(i+1) = x(i) + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
    t(i+1) = t(i) + h;
end
end
function [x1, x3] = CubicSolv(x, t)
%CubicSolv: Computes the roots to the polynomial
% t^2u^2 - 2xtu^2 + (1 + x^2)u - 1 = 0 for a specified pair of values (x,
% t). The roots are sorted so that Re(x1) <= Re(x2) <= Re(x3).
a = t^2;
b = -2*x*t;
c = (1+x^2);
d = -1;
r = roots([a b c d]);
r = sort(r, 'ComparisonMethod', 'real');
x1 = r(1);
x3 = r(3);
end