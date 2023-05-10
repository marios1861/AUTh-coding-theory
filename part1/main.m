clear
close all
clc

% rng(0, 'twister');

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

snr_db = -1:10;
snr = 10.^(snr_db/10);
iter = 5000;
iter2 = 2000;
np = size(snr, 2);

ber_golay = zeros(np, iter);
ber_cyclic = zeros(np, iter);

r_golay = double(12/23);
r_cyclic = double(4/7);


% Simulation for Golay

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

for j=1:iter
%     run1 = 1;
    for run1=1:length(snr)  
        s1 = snr(run1);
        decData = zeros(num_msg, size(msg, 2));
        p = qfunc(sqrt(2*r_golay*s1));

        for i = 1 : num_msg
            encData = my_encode(msg(i, :), gen_mat);
            cData = bsc(encData,p);
            decData(i, :) = my_decode(cData, syndrome_table, parity_check_matrix_t);
        end

        [numerr, pcterr] = biterr(msg, decData);
        ber_golay(run1,j) = pcterr;
%         run1 = run1 + 1;

        decData = cws(decData, bit_num);
        if(verLessThan("matlab","9.11"))
            decData = bi2de(decData,'left-msb')';
        else
            decData = bit2int(decData', bit_num);
        end
        unicodestr = native2unicode(decData,  'UTF-8');

%         fprintf("Number of messages sent: %d, Number of irrecoverable errors: %d, Pct Error: %f\n",...
%             size(msg, 1), numerr, pcterr);
%         fprintf(unicodestr);
    end
end


% Simulation for Cyclic

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

for j=1:iter2
    run2 = 1;
    for s2 = snr

        p = qfunc(sqrt(2*r_cyclic*s2));
        decData = zeros(num_msg, size(msg, 2));

        for i = 1 : num_msg
            encData = my_encode(msg(i, :), gen_mat);
            cData = bsc(encData,p);
            decData(i, :) = my_decode(cData, syndrome_table, parity_check_matrix_t);
        end

        [numerr, pcterr] = biterr(msg, decData);
        ber_cyclic(run2,j) = pcterr;
        run2 = run2 + 1;

        decData = cws(decData, bit_num);
        if(verLessThan("matlab","9.11"))
            decData = bi2de(decData,'left-msb')';
        else
            decData = bit2int(decData', bit_num);
        end
        unicodestr = native2unicode(decData,  'UTF-8');

%         fprintf("Number of messages sent: %d, Number of irrecoverable errors: %d, Pct Error: %f\n",...
%             size(msg, 1), numerr, pcterr);
%         fprintf(unicodestr);
    end
end

ber_golay_avg = mean(ber_golay, 2);
ber_cyclic_avg = mean(ber_cyclic, 2);
ber_uncoded = qfunc(sqrt(2.*snr));

semilogy(snr_db, ber_uncoded, "Color", "g", "linewidth", 1.2)
hold on;
semilogy(snr_db, ber_golay_avg, "Color", "r", "linewidth", 1.2)
hold on;
semilogy(snr_db, ber_cyclic_avg, "Color", "b", "linewidth", 1.2)
xlabel('Eb/N0 (dB)'); 
xlim([-1, snr_db(end)]);
ylim([10^(-6), 10^0]);
ylabel('Bit error rate');
title('Golay vs Cyclic Code')
legend("Uncoded", "Golay", "Cyclic", 'location','northeast');
grid on
