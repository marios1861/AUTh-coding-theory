clearvars;

x_0 = 0: 0.01 : 1;
epsilon = 0.5;
l_max = 8;
r_avg = 6;
r = floor(r_avg); % checked
syms rho(z);
rho(z) = (z .^ (r - 1) * (r * (r + 1 - r_avg))) / r_avg + ...
    (z.^r * (r_avg - r * (r + 1 - r_avg))) / r_avg; % checked
d_rho = diff(rho, z); % checked



prob = optimproblem('ObjectiveSense','max'); % checked
l = optimvar('lambda', l_max, 'LowerBound', 0); % checked
prob.Objective = sum(l ./ (1:l_max)'); % checked
prob.Constraints.cons1 = l(1) == 0; % checked
prob.Constraints.cons2 = sum(l) == 1;
fs = optimconstr(length(x_0));
f = @(x) epsilon * sum(l .* (1 - double(rho(1 - x))) .^ (0:(l_max - 1))') - x;
for i = 1 : length(x_0)
    fs(i) = f(x_0(i)) <= 0;
end

prob.Constraints.cons3 = fs;
prob.Constraints.cons4 = l(2) <= 1 / (epsilon * double(d_rho(1))); 


sol = solve(prob);

R = ((r * (r + 1 - r_avg)) / r_avg) / (r - 1)  + ((r_avg - r * (r + 1 - r_avg)) / r_avg) / r;

rate = 1 - R/sum(sol.lambda ./ (1:l_max)');
rate
sol.lambda
