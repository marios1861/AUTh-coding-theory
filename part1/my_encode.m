function encoded_msg = my_encode(msg, G)
    encoded_msg = mod(msg * G, 2);
end