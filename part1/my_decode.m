function decoded_msg = my_decode(msg, syndrome_table, H_t)
    % calc syndrome
    synd = syndrome(msg, H_t);
    % find syndrome in syndrome table
    str_synd = num2str(synd);
    str_synd(isspace(str_synd)) = '';
    synd_dec = bin2dec(str_synd);
    % find matching error
    if (synd_dec ~= 0)
        error = syndrome_table(synd_dec, :);
        % correct error
        decoded_msg = xor(error, msg);
    else
        decoded_msg = msg;
    end
    % slice off parity bits
    k = size(H_t, 1) - size(H_t, 2);
    decoded_msg = decoded_msg(1:k);
end