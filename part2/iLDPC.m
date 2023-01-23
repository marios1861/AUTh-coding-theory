clearvars;
clc;
clear;

desired_n = 6528;
n = 0;
l_max = 5;
r_avg = 14;
lambda =  [...
         0
    0.3846
    0.0921
         0
    0.5233
];

rho = zeros(r_avg, 1);
rho(r_avg - 1) = 1;

n_acc = [];
while (n < 100000)
    L = round((n * lambda ./ (1:l_max)') / sum(lambda ./ (1:l_max)'));
    R = round((n * rho ./ (1:r_avg)') / sum(lambda ./ (1:l_max)'));

    if(sum((1:length(L))' .* L) == sum((1:length(R))' .* R))
        n_acc = [n_acc; n];
    end

    n = n + 1;
end

[~,closest] = min(abs(desired_n-n_acc));

act_n = n_acc(closest);
L = round((act_n * lambda ./ (1:l_max)') / sum(lambda ./ (1:l_max)'));
R = round((act_n * rho ./ (1:r_avg)') / sum(lambda ./ (1:l_max)'));

k = act_n - sum(R);
r = k/act_n;
