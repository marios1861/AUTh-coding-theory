function [r, best_lambda] = design_iLDPC(epsilon, r_avg_max, l_max)
    r_avg_min = 2;
    rate = zeros(r_avg_max - r_avg_min + 1, 1);
    r_avgs = r_avg_min:r_avg_max;
    syms rho(z);

    for r_avg = r_avg_min:r_avg_max
        r = floor(r_avg);

        rho(z) = (z .^ (r - 1) * (r * (r + 1 - r_avg))) / r_avg + ...
            (z.^r * (r_avg - r * (r + 1 - r_avg))) / r_avg;
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

    [~, i] = max(rate);

    best_r_avg = r_avgs(i);
    r = floor(best_r_avg);
    rho(z) = (z .^ (r - 1) * (r * (r + 1 - best_r_avg))) / best_r_avg + ...
        (z.^r * (best_r_avg - r * (r + 1 - best_r_avg))) / best_r_avg;
    best_lambda = cdd_algorithm(rho, epsilon, l_max);