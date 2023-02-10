function q = init_q(c, H) 
    q = H';
    [ii,jj] = find(H');
    for k = 1:length(ii)
        q(ii(k), jj(k)) = c(ii(k)) * H(jj(k), ii(k));
    end
end