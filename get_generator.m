function generator = get_generator(P)
    generator = [eye(size(P, 1)), P];
end