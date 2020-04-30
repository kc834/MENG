fails = 0;
iter = 10;
for k = 1:iter
    disp(k);
    m = 1000;
    t = 3.5;
    tol = 1e-7;
    atol = 0.001;
    A = sprand(m,m,0.5) + 1i * sprand(m,m,0.5);
    % A = A'*A;
    % A = A/m;
%     A = sprand(m,m,0.005);
%     for j = 1:m
%         for i = j:m
%             if rand < 0.5, A(i,j) = 0; A(j,i) = 0; end
%         end
%     end
    % ishermitian(A)
    u = rand(m,1) + 3i *  rand(m,1);
    v = rand(m,1) + 3i * rand(m,1);
    if (utfAv_SingleArnoldi(u,A,t,v,tol,min(m,30)) > atol)
        disp("FAIL");
        fails = fails + 1;
        %disp(k);
    end
end
fail_rate = (fails)/iter;
disp(fail_rate);