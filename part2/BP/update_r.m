function r = update_r(q, H)

    function y = phi(x)
        if x == 0
            y = inf;
        else 
            y = 0;
        end
    end

    r = H;
    [jj,ii] = find(H);
    for k = 1:length(ii)
        [~, col_idx] = find(H(jj(k), :));
        sum_val = 0;
        for m = 1: length(col_idx)
            if col_idx(m) == ii(k)
                continue
            end
            q_ij = q(col_idx(m), jj(k));
            beta = abs(q_ij);
            sum_val = sum_val + phi(beta);
        end
        if (sum_val == inf)
            r(jj(k), ii(k)) = 0;
        else
            prod_val = 1;
            for m = 1: length(col_idx)
                if col_idx(m) == ii(k)
                    continue
                end
                q_ij = q(col_idx(m), jj(k));
                alpha = sign(q_ij);
                prod_val = prod_val * alpha * phi(sum_val);
            end
            r(jj(k), ii(k)) = prod_val;
        end
    end
end