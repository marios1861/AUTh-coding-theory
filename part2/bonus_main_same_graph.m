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
            c_data = zeros(1, length(c_data1));
            for j = 1:length(c_data1)
                if(c_data1(j) == -1)
                    c_data(j) = c_data2(j);
                else
                    c_data(j) = c_data1(j);
                end
            end
            dec_data = bp_decode(H, c_data);
            dec_data(dec_data == -1) = 1;
            [~, pcterr_i] = biterr(msg, dec_data);
            pct_err(i) = pcterr_i;
        end
        pct_errs(m) = mean(pct_err);
    end
    
    figure();
    semilogy(epsilons, pct_errs);
    xlabel('Channel erasure rate'); ylabel('Bit error rate');
    title(sprintf('Same MSG (Äe = %f)', dif_eps));
    grid on
%     saveas(gcf, sprintf("./results/%d_bonussame.png", dif_eps*100));
end