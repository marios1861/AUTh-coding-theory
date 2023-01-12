function new_msg = cws(msg, word_size)
    msg = reshape(msg',[], 1);
    old_size = size(msg);
    old_size = old_size(1);
    new_size = ceil(old_size / word_size ) * word_size;
    msg(old_size:new_size) = 0;
    new_msg = reshape(msg', word_size, [])';
end