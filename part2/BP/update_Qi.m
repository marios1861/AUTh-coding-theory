function Q = update_Qi(c, r, H)
    Q = c;
    for k = 1:length(Q)
        sum_val = 0;
        [row_idx, ~] = find(H(:, k));
        for m = 1: length(row_idx)
            r_ji = r(row_idx(m), k);
            sum_val = sum_val + r_ji;
        end
        Q(k) = Q(k) + sum_val;
    end

end