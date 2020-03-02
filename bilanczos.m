% [Q,R,alpha,beta,gamma] = lanczos(A,b) 
% 
% Compute a Bi-Lanczos decomposition 
% 
% A*Q = R*T 
%
% where T is a k+1-by-k tridiagonal matrix with diagonal 
% entries alpha, superdiagonal entries beta, and 
% subdiagonal entried gamma. Q and R have orthonormal 
% columns. 
% 
% Derived from lecture notes of David Bindel, CS 6210, Fall 2019, Cornell
% Algorithm from PhD Thesis of Axel Facius, July 2000, Universitat Karlsruhe 
%
function [Q,R,alpha,beta,gamma] = bilanczos(A,b,k)
 
n = length(A); % Max dim
Q = zeros(n,k+1); % Orthonormal basis, n by k+1 array
R = zeros(n,k+1);
alpha = zeros(k,1); 
beta = zeros(k,1);
gamma = zeros(k,1);

Q(:,1) = b/norm(b); % Arbitrary vector with norm 1
Q(:,1) = Q(:,1)/norm(Q(:,1));
R(:,1) = Q(:,1);
for j = 1:k 
    Q(:,j+1) = A*Q(:,j); % Move on to next vector in Krylov subspace
    R(:,j+1) = A'*R(:,j);
    alpha(j) = R(:,j)'*Q(:,j+1);
    Q(:,j+1) = Q(:,j+1)-alpha(j)*Q(:,j); 
    R(:,j+1) = R(:,j+1)-alpha(j)*R(:,j); 
    if j > 1 
        Q(:,j+1) = Q(:,j+1)-beta(j-1)*Q(:,j-1); 
        R(:,j+1) = R(:,j+1)-gamma(j-1)*R(:,j-1); 
    end
    gamma(j) = norm(Q(:,j+1));
    Q(:,j+1) = Q(:,j+1)/gamma(j);  
    beta(j) = R(:,j+1)'*Q(:,j+1);
    R(:,j+1) = R(:,j+1)/beta(j);
end

end