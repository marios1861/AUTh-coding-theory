function c_hat = bp_decode(H, Y, iter)  
% bp_decode  decode the received Y message
% using the BP algorithm.
%
%   C = bp_decode(H, Y, iter) decode y using the H parity
%   check matrix. Stop after iter loops or earlier if message 
%   is decoded earlier
%   
%   H: Sparse parity check matrix
%   Y: received message with erasure 
%   iter: int
 
    % Initialize L(c), L(q)
    c = init_c(Y); 
    q = init_q(c, H);
    
    for i=1:iter
        % Update L(r)
        r = update_r(q, H);
        % Update L(q)
        q = update_q(c, r, H);
        % Update L(Q)
        Q = update_Qi(c, r, H);
        % Estimate c_hat
        c_hat = estimate_c_hat(Q);
        
        if ~any(mod(c_hat*H', 2))
            break
        end
    end
end