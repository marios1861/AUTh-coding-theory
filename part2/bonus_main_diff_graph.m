clc;
clear;
addpath(genpath("./BP"));
epsilons = 0.05:0.05:0.5;
dif_epss = [0.02, 0.05, 0.1, 0.2, 0.3];
desired_n = 500;
glob_iter = 2000;
pct_errs = epsilons;
for dif_eps = dif_epss
    m = 1;
    for epsilon = epsilons
        l_max = 5;
        r_avg_max = 100;
        [n, k, L, R] = iLDPC(epsilon + dif_eps/2, desired_n, r_avg_max, l_max);
        H = zeros(n-k, n);
        r = find(R);
    
        % get number of graph edges
        edges = r * R(r);
        % ordered edges:
        E = 1:edges;
        % Matched edges to variable nodes
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
            
        H = sparse(H);
        msg = zeros(n, 1)';
        
        stat_iter = glob_iter/length(epsilons);
        pct_err = zeros(stat_iter, 1);
        
        parfor i = 1:stat_iter
            c_data1 = BEC(epsilon, msg);
            c_data2 = BEC(epsilon + dif_eps, msg);
            dec_data1 = bp_decode(H, c_data1, 100);
            dec_data2 = bp_decode(H, c_data2, 100);
            [~, pcterr_i1] = biterr(msg, dec_data1);
            [~, pcterr_i2] = biterr(msg, dec_data2);
            pct_err(i) = mean([pcterr_i1, pcterr_i2]);
        end
        pct_errs(m) = mean(pct_err);
        m = m + 1;
    end
    
    plot(epsilons, pct_errs);
    xlabel('Channel erasure rate'); ylabel('Bit error rate');
    title(sprintf('iLDPC BER on 2 BEC channels (Δe = %f) (n = %d)', dif_eps, desired_n));
    saveas(gcf, sprintf("./results/%d_bonusdiff.png", 100 * dif_eps));
end