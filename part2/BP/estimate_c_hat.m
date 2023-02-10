function c_hat = estimate_c_hat(Q)
    c_hat = Q;
    for i = 1:length(c_hat)
        if(c_hat(i) <= 0)
            c_hat(i) = 1;
        else
            c_hat(i) = 0;
        end
    end
    
end