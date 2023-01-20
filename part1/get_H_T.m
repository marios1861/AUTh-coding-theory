function H_T = get_H_T(P)
    H_T = [P', eye(size(P, 2))]';
end