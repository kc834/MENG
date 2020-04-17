 %for k = 1:1
    m = 2000;
    t = 1;
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
    if (utfAv_SingleArnoldi(u,A,t,v,tol,m) > atol)
        disp("FAIL");
        disp(k);
    end
 % end