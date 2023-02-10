function word = BEC(e, word)
    p = rand(length(word), 1);
    p = p < e;
    word(p) = -1;
end