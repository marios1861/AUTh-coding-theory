function syndrome_table = syndrome_decoding_init_1(H_t)
num_syndromes = size(H_t, 1);
syndromes = zeros(num_syndromes, size(H_t, 2));
errors = zeros(num_syndromes, size(H_t, 1));
l = 1;
for i = 1 : size(H_t, 1)
    err = zeros(size(H_t, 1), 1)';
    err(i) = 1;
    synd = syndrome(err, H_t);
    syndromes(l, :) = synd;
    errors(l, :) = err;
    l = l + 1;
end
[~, ia, ~] = unique(syndromes, 'rows');
syndrome_table = errors(ia, :);