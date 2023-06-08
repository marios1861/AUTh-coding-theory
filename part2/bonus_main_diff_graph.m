clc;
clear;
epsilon_max = 0.3;
epsilons = linspace(0.01, epsilon_max, 5);
dif_epss = [0.02, 0.05, 0.1];
n = 1000;
stat_iter = 5000;
pct_errs = epsilons;

l_max = 5;
r_avg_max = 10;
[r, L, R] = iLDPC(epsilon_max + mean(dif_epss/2), n, r_avg_max, l_max);
k = round(r * n);
H = zeros(n-k, n);
r = find(R);

% Get number of graph edges
edges = r' * R(r);
% ordered edges:
E = 1:edges;
% Match edges to variable nodes
ME = zeros(edges,1);
for i = 1:length(L)
    if L(i) ~= 0
        prev_edges = sum((1:i - 1) * L(1:i - 1));
        prev_nodes = sum(L(1:i-1));
        % Match i * L(i) edges to L(i) variable nodes
        ME(prev_edges + 1 : prev_edges + L(i) * i) = ...
            repmat(((prev_nodes + 1) : L(i) + prev_nodes)', i, 1);
    end
end
ME = ME(randperm(edges));
% Match ordered edges to check nodes (regular matching)
E = mod(E, n-k) + 1;

for i = 1:edges
    H(E(i), ME(i)) = mod(H(E(i), ME(i)) + 1, 2);
end
msg = zeros(n, 1)';
for dif_eps = dif_epss
    parfor m = 1:length(epsilons)
        pct_err = zeros(stat_iter, 1);
        
        for i = 1:stat_iter
            c_data1 = BEC(epsilons(m), msg);
            c_data2 = BEC(epsilons(m) + dif_eps, msg);
            dec_data1 = bp_decode(H, c_data1);
            dec_data1(dec_data1 == -1) = 1;
            dec_data2 = bp_decode(H, c_data2);
            dec_data2(dec_data2 == -1) = 1;
            [~, pcterr_i1] = biterr(msg, dec_data1);
            [~, pcterr_i2] = biterr(msg, dec_data2);
            pct_err(i) = mean([pcterr_i1, pcterr_i2]);
        end
        pct_errs(m) = mean(pct_err);
    end
    
    figure();
    semilogy(epsilons, pct_errs);
    xlabel('Channel erasure rate'); ylabel('Bit error rate');
    title(sprintf('Diff MSG (Äe = %f)', dif_eps));
    grid on
    saveas(gcf, sprintf("./results/%d_bonussame.png", dif_eps*100));
end