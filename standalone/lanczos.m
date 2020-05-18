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
% Taken from lecture notes of David Bindel, CS 6210, Fall 2019, Cornell
%
% function [Q,alpha,beta] = lanczos(A,b,k)
function out = lanczos(A,v,u,k)

% n = length(A); % Max dim
% Q = zeros(n,k+1); % Orthonormal basis, n by k+1 array 
% alpha = zeros(k,1); 
% beta = zeros(k,1);
% 
% Q(:,1) = b/norm(b); % Arbitrary vector with norm 1
% for j = 1:k 
%     Q(:,j+1) = A*Q(:,j); % Move on to next vector in Krylov subspace 
%     alpha(j) = Q(:,j)'*Q(:,j+1);  
%     Q(:,j+1) = Q(:,j+1)-alpha(j)*Q(:,j); 
%     if j > 1 
%         Q(:,j+1) = Q(:,j+1)-beta(j-1)*Q(:,j-1); 
%     end
%     beta(j) = norm(Q(:,j+1));
%     Q(:,j+1) = Q(:,j+1)/beta(j); 
% end

n = length(A); % Max dim
Q = zeros(n,k+1); % Orthonormal basis, n by k+1 array 
T = zeros(k+1, k);
% alpha = zeros(k,1); 
% beta = zeros(k,1);

Q(:,1) = v/norm(v); % Arbitrary vector with norm 1
for j = 1:k 
    Q(:,j+1) = A*Q(:,j); % Move on to next vector in Krylov subspace 
    T(j,j) = Q(:,j+1)'*Q(:,j);  
    Q(:,j+1) = Q(:,j+1)-T(j,j)*Q(:,j); 
    if j > 1 
        Q(:,j+1) = Q(:,j+1)-T(j,j-1)*Q(:,j-1); 
        T(j-1,j) = T(j,j-1);
    end
    if norm(Q(:, j+1)) < 1e-7 
        break;
    end
    
    T(j+1,j) = norm(Q(:,j+1));
    Q(:,j+1) = Q(:,j+1)/T(j+1,j); 
end

% disp(Q);
% disp(T);
% disp(Q*T);
% disp(A*Q(:, 1:k));
% eq = all(ismembertol(Q*T,A*Q(:, 1:k), 1e-7), 'all');
% if eq
%     disp("good");
% else
%     disp("bad");
% end

w = Q(:,1:k)*(norm(v)*T(1:k,1));
out = u'*w;

end