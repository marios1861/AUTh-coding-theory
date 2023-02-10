function q = update_q(c, r, H)
    q = H';
    [ii,jj] = find(H');
    for k = 1:length(ii)
        sum_val = c(ii(k)) * H(jj(k), ii(k));
        [row_idx, ~] = find(H(:, ii(k)));
        for m = 1: length(row_idx)
            if row_idx(m) == jj(k)
                continue
            end
            r_ji = r(row_idx(m), ii(k));
            sum_val = sum_val + r_ji;
        end
        q(ii(k), jj(k)) = sum_val;
    end
end