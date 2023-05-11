clc;
clear;
addpath(genpath("./BP"));
epsilons = 0.01:0.05:0.30;
ns = [250, 500, 1000];
glob_iter = 100;
pct_errs = epsilons;
for n = ns
    m = 1;
    for epsilon = epsilons
        l_max = 5;
        r_avg_max = 100;
        [r, L, R] = iLDPC(epsilon, n, r_avg_max, l_max);
        k = round(r * n);
        H = zeros(n-k, n);
        r = find(R);
    
        % get number of graph edges
        edges = r' * R(r);
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
        % Match unordered edges to check nodes (regular matching)
        E = mod(E, n-k) + 1;
    
        for i = 1:edges
            H(E(i), ME(i)) = mod(H(E(i), ME(i)) + 1, 2);
        end
        
        msg = zeros(n, 1)';
        
        stat_iter = ceil(glob_iter/length(epsilons));
        pct_err = zeros(stat_iter, 1);
        
        for i = 1:stat_iter
            c_data = BEC(epsilon, msg);
            dec_data = bp_decode(H, c_data);
            % change uncertainty into error to be valid for biterr func
            dec_data(dec_data == -1) = 1;
            [~, pcterr_i] = biterr(msg, dec_data);
            pct_err(i) = pcterr_i;
        end
        pct_errs(m) = mean(pct_err);
        m = m + 1;
    end
    
    plot(epsilons, pct_errs);
    xlabel('Channel erasure rate'); ylabel('Bit error rate');
    title(sprintf('iLDPC BER on BEC (n = %d)', n));
    saveas(gcf, sprintf("./results/%d_%d_iLDPCGRAPH.png", n, l_max));
end