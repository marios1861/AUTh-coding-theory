% m has to be a row vector and H the transpose of the parity matrix
function y = syndrome(m, H_t)
    y = mod(m * H_t, 2);
end