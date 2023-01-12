clear all %#ok<CLALL>
close all
clc


golay_matrix = [
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
    [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1];
    [1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0];
    [1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1];
    [1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1];
    [1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0];
    [1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1];
    [1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1];
    [1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1];
    [1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0];
    [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0];
    [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0]];

p = 0.12;
seed = 1;

parity_check_matrix_t = get_H_T(golay_matrix);
gen_mat = get_generator(golay_matrix);
syndrome_table = syndrome_decoding_init(parity_check_matrix_t);

fileID = fopen("data.txt", 'r');
str = fscanf(fileID, '%c');
msg = unicode2native(str);
bit_num = ceil(log2(double(max(msg))));
msg = double(int2bit(msg, bit_num)');
msg = cws(msg, size(golay_matrix, 1));
num_msg = size(msg, 1);

decData = zeros(num_msg, size(msg, 2));

for i = 1 : num_msg
    encData = my_encode(msg(i, :), gen_mat);
    cData = bsc(encData,p);
    decData(i, :) = my_decode(cData, syndrome_table, parity_check_matrix_t);
end

[numerr, pcterr] = biterr(msg, decData);


decData = cws(decData, bit_num);
decData = bit2int(decData', bit_num);
unicodestr = native2unicode(decData);

fprintf("Number of messages sent: %d, Number of irrecoverable errors: %d, Pct Error: %f\n",...
    size(msg, 1), numerr, pcterr);

fprintf(unicodestr);