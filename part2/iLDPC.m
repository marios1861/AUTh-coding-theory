function [r, L_hat, R_hat] = iLDPC(epsilon, n, r_avg_max, l_max)
    [r_avg, lambda] = design_iLDPC(epsilon, r_avg_max, l_max);
    
    L = (n * lambda ./ (1:l_max)') / sum(lambda ./ (1:l_max)');
    L_hat = floor(L);
    rho = zeros(r_avg, 1);
    rho(r_avg - 1) = 1;
    R = (n * rho ./ (1:r_avg)') / sum(lambda ./ (1:l_max)');
    R_hat = floor(R);
    r = 1 - sum(R)/sum(L);
    A = (1 - r)*n - sum(R_hat);
    
    % The vars are stored as so xl2, ..., xlmax, xr2, ..., xr_avg    
    f = ones(r_avg + l_max - 2,1);
    intcon = (1:r_avg + l_max - 2)';
    Ain = [zeros(1, l_max - 1), -ones(1, r_avg -1)];
    bin = -ceil(A);
    Aeq = [ones(1, l_max - 1), zeros(1, r_avg -1);
           (2:l_max), -(2:r_avg)];
    beq = [n - sum(L_hat);
           (1:r_avg) * R_hat - (1:l_max) * L_hat];
    lb = zeros(r_avg + l_max -2, 1);
    ub = f;

    opts = optimoptions(@intlinprog,'Display','off');
    [x_min, fval1, exitflag1] = intlinprog(f, intcon, Ain, bin, Aeq, beq, lb, ub, [], opts);
    
    [x_max, fval2, exitflag2] = intlinprog(-f, intcon, -Ain, -bin, Aeq, beq, lb, ub, [], opts);
    
    [exitflag1, exitflag2]

    % select better choice

    choice1 = abs(fval1 - A);
    choice2 = abs(fval2 - A);

    if(choice1 < choice2)
        x_hat = x_min;
    else
        x_hat = x_max;
    end

    L_hat(2:end) = L_hat(2:end) + x_hat(1: l_max - 1);
    R_hat(2:end) = R_hat(2:end) + x_hat(l_max: end);
    r = 1 - sum(R_hat)/sum(L_hat);
end