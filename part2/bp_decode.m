function c_hat = bp_decode(H, y)  
% bp_decode  decode the received Y message
% using the BP algorithm.
%
%   c_hat = bp_decode(H, Y, iter) decode y using the H parity
%   check matrix. Stop after iter loops or earlier if message 
%   is decoded earlier
%   
%   H: Sparse parity check matrix
%   Y: received message with erasure 
%   iter: int
 
    % Initialize c_hat
    c_hat = y;

    endflag = 1;
    % var-to-check and check-to-var loop
    while(1)
        for line=H'
            indexes = 1:length(c_hat);
            indexes = indexes(logical(line));
            connected_vars = c_hat(logical(line));
            if sum(connected_vars == -1) == 1
                if endflag == 1
                    endflag = 0;
                end
                erasure_idx = indexes(connected_vars == -1);
                c_hat(erasure_idx) = mod(sum(connected_vars(connected_vars ~= -1)), 2);
            end
        end
        if sum(c_hat == -1) == 0
            break
        end
        if endflag == 0
            endflag = 1;
        else
            break
        end
    end
end