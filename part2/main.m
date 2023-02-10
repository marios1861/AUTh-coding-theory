clc;
clear;
addpath(genpath("./BP"));

desired_n = 1000;
l_max = 5;
epsilon = 0.2;
r_avg_max = 100;
[n, k, L, R] = iLDPC(epsilon, desired_n, r_avg_max, l_max);
H = zeros(n-k, n);
n_i = 1; 
for i = 1:length(L)
    if (L(i) == 0)
        continue
    end
    for j = 1:L(i)
        H(1:i, n_i) = 1;
        n_i = n_i + 1;
    end
end

% mix up rows
H = H(randperm(n-k), :);
% mix up cols
H = H(:, randperm(n));

H = sparse(H);
msg = zeros(n, 1)';

stat_iter = 32;
pct_err = zeros(stat_iter, 1);
parfor i = 1:stat_iter
    c_data = BEC(epsilon, msg);
    dec_data = bp_decode(H, c_data, k/2);
    [~, pcterr_i] = biterr(msg, dec_data);
    pct_err(i) = pcterr_i;
    disp(i)
end

histogram(pct_err, 100, 'Normalization','probability')

