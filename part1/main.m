clear
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

cyclic_matrix = [
    1 1 0;
    0 1 1;
    1 0 1;
    1 1 1];

p_data = 0 : 0.01 : 0.5;
np = size(p_data, 2);
ber_golay = zeros(np, 1);
ber_cyclic = zeros(np, 1);
run = 1;

parity_check_matrix_t = get_H_T(golay_matrix);

gen_mat = get_generator(golay_matrix);
syndrome_table = syndrome_decoding_init_golay(parity_check_matrix_t);

fileID = fopen("data.txt", 'r');
str = fscanf(fileID, '%c');
msg = unicode2native(str,  'UTF-8');
bit_num = ceil(log2(double(max(msg))));

if(verLessThan("matlab","9.11"))
    msg = double(de2bi(msg, bit_num, 'left-msb'));
else
    msg = double(int2bit(msg, bit_num)');
end
msg = cws(msg, size(golay_matrix, 1));
num_msg = size(msg, 1);

for p =  p_data    
    decData = zeros(num_msg, size(msg, 2));
    
    for i = 1 : num_msg
        encData = my_encode(msg(i, :), gen_mat);
        cData = bsc(encData,p);
        decData(i, :) = my_decode(cData, syndrome_table, parity_check_matrix_t);
    end
    
    [numerr, pcterr] = biterr(msg, decData);
    ber_golay(run) = pcterr;
    run = run + 1;
    
    decData = cws(decData, bit_num);
    if(verLessThan("matlab","9.11"))
        decData = bi2de(decData,'left-msb')';
    else
        decData = bit2int(decData', bit_num);
    end
    unicodestr = native2unicode(decData,  'UTF-8');
    
    %     fprintf("Number of messages sent: %d, Number of irrecoverable errors: %d, Pct Error: %f\n",...
    %         size(msg, 1), numerr, pcterr);
    
    %     fprintf(unicodestr);
end


run = 1;
parity_check_matrix_t = get_H_T(cyclic_matrix);
gen_mat = get_generator(cyclic_matrix);
syndrome_table = syndrome_decoding_init_1(parity_check_matrix_t);

fileID = fopen("data.txt", 'r');
str = fscanf(fileID, '%c');
msg = unicode2native(str,  'UTF-8');
bit_num = ceil(log2(double(max(msg))));

if(verLessThan("matlab","9.11"))
    msg = double(de2bi(msg, bit_num, 'left-msb'));
else
    msg = double(int2bit(msg, bit_num)');
end
msg = cws(msg, size(cyclic_matrix, 1));
num_msg = size(msg, 1);
for p =  p_data

    decData = zeros(num_msg, size(msg, 2));
    
    for i = 1 : num_msg
        encData = my_encode(msg(i, :), gen_mat);
        cData = bsc(encData,p);
        decData(i, :) = my_decode(cData, syndrome_table, parity_check_matrix_t);
    end
    
    [numerr, pcterr] = biterr(msg, decData);
    ber_cyclic(run) = pcterr;
    run = run + 1;
    
    decData = cws(decData, bit_num);
    if(verLessThan("matlab","9.11"))
        decData = bi2de(decData,'left-msb')';
    else
        decData = bit2int(decData', bit_num);
    end
    unicodestr = native2unicode(decData,  'UTF-8');
    
    %     fprintf("Number of messages sent: %d, Number of irrecoverable errors: %d, Pct Error: %f\n",...
    %         size(msg, 1), numerr, pcterr);
    
    %     fprintf(unicodestr);
end

plot(p_data, ber_golay, "Color", "r")
hold on;
plot(p_data, p_data, "Color", "g")
hold on;
plot(p_data, ber_cyclic, "Color", "b")
xlabel('Channel error rate'); ylabel('Bit error rate');
title('Golay vs Cyclic Code')
legend("Golay", "No Code", "Cyclic", 'location','northwest');
