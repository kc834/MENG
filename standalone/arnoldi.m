%[Q,H] = arnoldi(A,b)
%
% Compute an Arnoldi decomposition
%
% A*Q(:,1:end-1) = Q*H
%
% where H is a k+1-by-k upper Hessenberg matrix and Q has
% orthonormal columns.
%
function out = arnoldi(A,v,u, k)

n = length(A);
Q = zeros(n,k+1); % Orthonormal basis
H = zeros(k+1,k); % Upper Hessenberg matrix

Q(:,1) = v/norm(v);
for j = 1:k

    % Get a vector in the next subspace (and its norm)
    Q(:,j+1) = A*Q(:,j);
    % norma = norm(Q(:,j+1));

    % Modified Gram-Schmidt (standard Arnoldi)
    for l = 1:j
        H(l,j) = Q(:,l)'*Q(:,j+1);
        Q(:,j+1) = Q(:,j+1)-Q(:,l)*H(l,j);
    end

    if norm(Q(:, j+1)) < 1e-7 
        break;
    end
    H(j+1,j) = norm(Q(:,j+1));

    % Normalize final result
    Q(:,j+1) = Q(:,j+1)/H(j+1,j);

end

w = Q(:,1:k)*(norm(v)*H(1:k,1));
out = u'*w;
end