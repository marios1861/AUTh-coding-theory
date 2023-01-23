clearvars;
clc;

epsilon = 0.2;
r_avg_min = 2;
r_avg_max = 100;
rate = zeros(r_avg_max - r_avg_min + 1, 1);
start = 1;
r_avgs = r_avg_min:r_avg_max;
syms rho(z);
l_max = 5;

for r_avg = r_avg_min:r_avg_max
    r = floor(r_avg); % checked

    rho(z) = (z .^ (r - 1) * (r * (r + 1 - r_avg))) / r_avg + ...
        (z.^r * (r_avg - r * (r + 1 - r_avg))) / r_avg; % checked
    lambda = cdd_algorithm(rho, epsilon, l_max);
    R = ((r * (r + 1 - r_avg)) / r_avg) / (r - 1)  + ((r_avg - r * (r + 1 - r_avg)) / r_avg) / r;

    if ~isempty(lambda)
        rate(r_avg - r_avg_min  + 1) = 1 - R/sum(lambda ./ (1:l_max)');
    else
        break;
    end
end

plot(r_avgs, rate);
xlabel("r avg");
ylabel("Code Rate");
title(sprintf("for Lmax = %d", l_max));
ylim([0 1]);
xlim([r_avg_min, r_avg + 1]);

[best_rate, i] = max(rate);

best_r_avg = r_avgs(i)
r = floor(best_r_avg); % checked
rho(z) = (z .^ (r - 1) * (r * (r + 1 - best_r_avg))) / best_r_avg + ...
    (z.^r * (best_r_avg - r * (r + 1 - best_r_avg))) / best_r_avg;
best_lambda = cdd_algorithm(rho, epsilon, l_max)