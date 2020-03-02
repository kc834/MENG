% [Q,alpha,beta] = lanczos(A,b) 
% 
% Compute a Lanczos decomposition 
% 
% A*Q = Q*T 
% 
% where T is a k+1-by-k tridiagonal matrix with diagonal 
% entries alpha and super/subdiagonals beta, and Q has
% orthonormal columns. 
% 
function [Q,alpha,beta] = lanczos(A,b,k)
 
n = length(A); % Max dim
Q = zeros(n,k+1); % Orthonormal basis, n by k+1 array 
alpha = zeros(k,1); 
beta = zeros(k,1);

Q(:,1) = b/norm(b); % Arbitrary vector with norm 1
for j = 1:k 
    Q(:,j+1) = A*Q(:,j); % Move on to next vector in Krylov subspace 
    alpha(j) = Q(:,j)'*Q(:,j+1);  
    Q(:,j+1) = Q(:,j+1)-alpha(j)*Q(:,j); 
    if j > 1 
        Q(:,j+1) = Q(:,j+1)-beta(j-1)*Q(:,j-1); 
    end
    beta(j) = norm(Q(:,j+1));
    Q(:,j+1) = Q(:,j+1)/beta(j); 
end

end