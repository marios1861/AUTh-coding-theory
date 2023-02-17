function [act_n, k, L, R] = iLDPC(epsilon, desired_n, r_avg_max, l_max)
    n = round(desired_n/10);

    [r_avg, lambda] = design_iLDPC(epsilon, r_avg_max, l_max);
    rho = zeros(r_avg, 1);
    rho(r_avg - 1) = 1;

    n_acc = [];
    while (n < 10*desired_n)
        L = round((n * lambda ./ (1:l_max)') / sum(lambda ./ (1:l_max)'));
        R = round((n * rho ./ (1:r_avg)') / sum(lambda ./ (1:l_max)'));

        if(sum((1:length(L))' .* L) == sum((1:length(R))' .* R))
            n_acc = [n_acc; n];
        end

        n = n + 1;
    end
    % find valid n vals
    sum_R = zeros(length(n_acc), 1);
    for i = 1:length(n_acc)
        sum_R(i) = sum(round((n_acc(i) * rho ./ (1:r_avg)') / sum(lambda ./ (1:l_max)')));
    end
    n_acc(sum_R < l_max) = inf;
    [~,closest] = min(abs(desired_n-n_acc));

    act_n = n_acc(closest);
    L = round((act_n * lambda ./ (1:l_max)') / sum(lambda ./ (1:l_max)'));
    R = round((act_n * rho ./ (1:r_avg)') / sum(lambda ./ (1:l_max)'));

    k = act_n - sum(R);
    r = k/act_n;
