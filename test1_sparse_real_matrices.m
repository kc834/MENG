fails = 0; 
iter = 100;
for k = 1:iter
    disp(k);
    m = 800;
    t = 3.5;
    tol = 1e-7;
    atol = 0.001;
    A = rand(m);
    for j = 1:m
    for i = 1:m
    if rand < 0.5, A(i,j) = 0; end
    end
    end
    A = sparse(A);
    u = rand(m,1);
    v = rand(m,1);
    out = utfAv_SingleArnoldi(u,A,t,v,tol,30);
    while (out > atol)
        disp("RESTARTING");
        out = utfAv_SingleArnoldi(u,A,t,v,tol,30);
        fails = fails + 1;
    end
end
fail_rate = (fails)/iter;
disp(fail_rate);