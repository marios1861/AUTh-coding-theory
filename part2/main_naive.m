clc;
clear;
addpath(genpath("./BP"));
epsilons = 0.05:0.05:0.5;
desired_ns = [100, 250, 500, 1000];
glob_iter = 1000;
pct_errs = epsilons;
for desired_n = desired_ns
    m = 1;
    for epsilon = epsilons
        l_max = 10;
        r_avg_max = 100;
        [n, k, L, R] = iLDPC(epsilon, desired_n, r_avg_max, l_max);
        H = zeros(n-k, n);
        n_i = 1; 
        r = find(R);
        h = 1;
        for i = 1:length(L)
            if (L(i) == 0)
                continue
            end
    
            for j = 1:L(i)
                while(sum(H(h, :)) == r)
                    h = h + 1;
                end
                H(h: min(h - 1 + i, n - k), n_i) = 1;
                n_i = n_i + 1;
            end
        end
        
        % mix up rows 
        H = H(randperm(n-k), :);
    
        % mix up cols
        H = H(:, randperm(n));
        
        H = sparse(H);
        msg = zeros(n, 1)';
        
        stat_iter = glob_iter/length(epsilons);
        pct_err = zeros(stat_iter, 1);
        
        parfor i = 1:stat_iter
            c_data = BEC(epsilon, msg);
            dec_data = bp_decode(H, c_data, 100);
            [~, pcterr_i] = biterr(msg, dec_data);
            pct_err(i) = pcterr_i;
        end
        pct_errs(m) = mean(pct_err);
        m = m + 1;
    end
    
    plot(epsilons, pct_errs);
    xlabel('Channel erasure rate'); ylabel('Bit error rate');
    title(sprintf('iLDPC BER on BEC (n = %d)', desired_n));
    saveas(gcf, sprintf("./results/%d_%d_iLDPCNAIVE.png", desired_n, l_max));
end