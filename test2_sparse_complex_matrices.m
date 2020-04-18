m = 40;
t = 1;
tol = 1e-7;
atol = 0.001;
A = rand(m) + 1i * rand(m);
A = A'*A;
A = A/m;
for j = 1:m
for i = j:m
if rand < 0.5, A(i,j) = 0; A(j,i) = 0; end
end
end
% ishermitian(A)
u = eye(m,1) + 2i *  eye(m,1);
v = eye(m,1) + 3i * eye(m,1);
if (utfAv_SingleArnoldi(u,A,t,v,tol,m) > atol)
    disp("FAIL");
    %disp(k);
end