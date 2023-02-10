function c = init_c(Y)
    function ci = log_bec(yi)
        if yi == 0
            ci = inf;
        elseif yi == 1
            ci = -inf;
        else 
            ci = 0;
        end
    end

    c = arrayfun(@log_bec, Y);
end

