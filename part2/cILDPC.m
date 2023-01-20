clearvars;

x_0 = 0.2;
epsilon = 0.2;
l_max = 8;
r_avg = 6;
r = floor(r_avg); % checked
syms rho(z);
rho(z) = (z ^ (r - 1) * (r * (r + 1 - r_avg))) / r_avg + ...
    (z^r * (r_avg - r * (r + 1 - r_avg))) / r_avg; % checked
d_rho = diff(rho, z); % checked



prob = optimproblem('ObjectiveSense','max'); % checked
l = optimvar('lambda', l_max, 'LowerBound', 0); % checked
prob.Objective = sum(l ./ (1:l_max)'); % checked
prob.Constraints.cons1 = l(1) == 0; % checked
prob.Constraints.cons2 = sum(l) == 1;
f = @(x) epsilon * sum(l .* (1 - double(rho(1 - x))) .^ (0:(l_max - 1))') - x;
prob.Constraints.cons3_2 = f(x_0) <= 0;
prob.Constraints.cons4 = l(2) <= 1 / (epsilon * double(d_rho(1))); 

show(prob)

sol = solve(prob);

R = ((r * (r + 1 - r_avg)) / r_avg) / (r - 1)  + ((r_avg - r * (r + 1 - r_avg)) / r_avg) / r;

rate = 1 - R/sum(sol.lambda ./ (1:l_max)');
rate
sol.lambda
