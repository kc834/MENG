fails = 0; 
iter = 1;
for k = 1:iter
    disp(k);
    m = 4781;
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
    if (compare(u,A,t,v,tol,min(m,30)) > atol)
        disp("FAIL");
        fails = fails + 1;
        %disp(k);
    end
end
fail_rate = (fails)/iter;
disp(fail_rate);