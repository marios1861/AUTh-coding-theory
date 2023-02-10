assert_test();
disp("Decoded message is correct!");

function assert_test()
    rng(1)
    H = [1 1 0 1 1 0 0;
         1 0 1 1 0 1 0;
         0 1 1 1 0 0 1];
    H = sparse(H);
    y = [0 -1 -1 1 0 -1 0];
    c_hat = bp_decode(H, y, 10);
    assert(all(c_hat == [0 1 0 1 0 1 0]), "Decoded message is not correct")
end