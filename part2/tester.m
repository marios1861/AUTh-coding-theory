rng(1)
H = [1 1 1 0 1;
     0 1 1 1 0];
y = [0 1 -1 0 -1];
c_hat = bp_decode(H, y, 10);
assert(all(c_hat == [0 1 1 0 0]), "WRONG")

H = [1 0 0 1 1 1 0;
     0 1 0 1 1 0 1;
     0 0 1 1 0 1 1];
y = [0 -1 0 1 -1 -1 0];
c_hat = bp_decode(H, y, 10);
assert(all(c_hat == [0 1 0 1 0 1 0]), "WRONG")