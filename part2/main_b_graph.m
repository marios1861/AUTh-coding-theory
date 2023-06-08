clc;
clear;

epsilon_max = 0.3;
epsilons = linspace(0.01, epsilon_max, 5);
ns = 1000;
pct_errs = epsilons;
stat_iter = 5000;

for n = ns
    l_max = 5;
    r_avg_max = 10;
    [r, L, R] = iLDPC(epsilon_max, n, r_avg_max, l_max);
    k = round(r * n);
    r = find(R);
    % Get number of graph edges
    edges = r' * R(r);
    % Round edges so variable nodes have equal edges
    edges = round(edges / n) * n;
    for k_cand = k:-1:0
        if floor(edges/(n-k_cand)) == edges/(n-k_cand)
            k = k_cand;
            break;
        end
    end
    % ordered edges:
    E = 1:edges;
    % Match edges to variable nodes
    ME = mod(randperm(edges), n) + 1;
    % Match ordered edges to check nodes (regular matching)
    E = mod(E, n-k) + 1;
    
    H = zeros(n-k, n);
    for i = 1:edges
        H(E(i), ME(i)) = mod(H(E(i), ME(i)) + 1, 2);
    end
        
    msg = zeros(n, 1)';
    
    parfor m = 1:length(epsilons)
        
        pct_err = zeros(stat_iter, 1);
        
        for i = 1:stat_iter
            c_data = BEC(epsilons(m), msg);
            dec_data = bp_decode(H, c_data);
            % change uncertainty into error to be valid for biterr func
            dec_data(dec_data == -1) = 1;
            [~, pcterr_i] = biterr(msg, dec_data);
            pct_err(i) = pcterr_i;
        end
        pct_errs(m) = mean(pct_err);
    end
    
    figure(2);
    semilogy(epsilons, pct_errs);
    xlabel('Channel erasure rate'); ylabel('Bit error rate');
    title(sprintf('LDPC BER on BEC'));
    grid on
    saveas(gcf, sprintf("./results/LDPCGRAPH.png"));
end