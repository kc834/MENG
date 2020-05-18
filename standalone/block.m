% out = block(A,v,u,k)
% 
% Compute a Block-Lanczos decomposition given nxn matrix A, nxp vectors
% v and u, and dimension k
% 
% Taken from ABLE: An Adaptive Block Lanczos Method for Non-Hermitian 
% Eigenvalue Problems

function out = block(A,v,u,k)

[n, p] = size(v);
Q = zeros(n,k+p); % Orthonormal basis, n by k+1 array
P = zeros(n,k+p);
R = zeros(n,k+p);
S = zeros(n,k+p);
T = zeros(k+p, k+p);

Q(:,1:p) = v;
normv = norm(Q(:,1:p)); % Get the norm of each column
Q(:,1:p) = Q(:,1:p)/normv; % Divide each value by column norm
P(:,1:p) = u;
normuv = Q(:,1:p)'*P(:,1:p);
P(:,1:p) = P(:,1:p)/normuv; % Faster than multiplying by inverse, probably need to multiply back normv
nblocks = k / p;
R(:, 1:p) = A'*P(:, 1:p);
S(:, 1:p) = A*Q(:,1:p);
for j = 1:nblocks
    low = (j-1)*p + 1; nxtlow = j*p + 1;
    hi = j * p; nxthi = (j+1) * p;
    T(low:hi, low:hi) = P(:, low:hi)' * S(:, low:hi);
    R(:, low:hi) = R(:, low:hi) - P(:, low:hi) * T(low:hi, low:hi)';
    S(:, low:hi) = S(:, low:hi) - Q(:, low:hi) * T(low:hi, low:hi);
    % QR Factorization
    [Qtemp, Rtemp] = qr(R(:, low:hi));
    P(:, nxtlow:nxthi) = Qtemp(:, low:hi);
    T(low:hi, nxtlow:nxthi)= Rtemp(1:p, 1:p)';
    [Qtemp, Rtemp] = qr(S(:, low:hi));
    Q(:, nxtlow:nxthi) = Qtemp(:, low:hi);
    T(nxtlow:nxthi, low:hi) = Rtemp(1:p, 1:p);
    % SVD Factorization (don't need to save values)
    [U, Sigma, V] = svd(P(:, nxtlow:nxthi)' * Q(:, nxtlow:nxthi));
    rootSigma = sqrtm(Sigma);
    T(low:hi, nxtlow:nxthi) = T(low:hi, nxtlow:nxthi) * U * rootSigma;
    T(nxtlow:nxthi, low:hi) = rootSigma * V' * T(nxtlow:nxthi, low:hi);
    P(:, nxtlow:nxthi) = (P(:, nxtlow:nxthi) * conj(U))/rootSigma;
    Q(:, nxtlow:nxthi) = (Q(:, nxtlow:nxthi) * V)/rootSigma;
    R(:, nxtlow:nxthi) = (P(:, nxtlow:nxthi)'*A - T(nxtlow:nxthi, low:hi)*P(:, low:hi)')';
    S(:, nxtlow:nxthi) = A*Q(:, nxtlow:nxthi) - Q(:, low:hi)*T(low:hi, nxtlow:nxthi);
end

% disp(Q);
% disp(P);
% disp(T);
% disp(Q(:, 1:k+1)*T(1:k+1, 1:k));
% disp(A*Q(:, 1:k));
% eqQ = all(ismembertol(Q(:, 1:k+1)*T(1:k+1, 1:k),A*Q(:, 1:k), 1e-7), 'all');
% if eqQ
%     disp("good eqQ");
% else
%     disp("bad eqQ");
% end
% 
% disp(P(:, 1:k+1)*T(1:k, 1:k+1)');
% disp(A'*P(:, 1:k));
% eqP = all(ismembertol(P(:, 1:k+1)*T(1:k, 1:k+1)',A'*P(:, 1:k), 1e-7), 'all');
% if eqP
%     disp("good eqP");
% else
%     disp("bad eqP");
% end

out = normuv' * normv * T(1:p, 1:p);

% disp(Q(:, 1:k)'*P(:, 1:k));
% disp(P(:, :)'*Q(:, :));
% eqI = all(ismembertol(Q(:, 1:k)'*P(:, 1:k),P(:, 1:k)'*Q(:, 1:k), 1e-7), 'all');
% if eqI
%     disp("good eqI");
% else
%     disp("bad eqI");
% end

%  ======================== OLD CODE =======================

% n = length(A); % Max dim
% Q = zeros(n,k+1); % Orthonormal basis, n by k+1 array
% P = zeros(n,k+1);
% alpha = zeros(k,1); 
% beta = zeros(k,1);
% gamma = zeros(k,1);
% 
% Q(:,1) = b/norm(b); % Arbitrary vector with norm 1
% Q(:,1) = Q(:,1)/norm(Q(:,1));
% P(:,1) = Q(:,1);
% for j = 1:k 
%     Q(:,j+1) = A*Q(:,j); % Move on to next vector in Krylov subspace
%     P(:,j+1) = A'*P(:,j);
%     alpha(j) = P(:,j)'*Q(:,j+1);
%     Q(:,j+1) = Q(:,j+1)-alpha(j)*Q(:,j); 
%     P(:,j+1) = P(:,j+1)-alpha(j)*P(:,j); 
%     if j > 1 
%         Q(:,j+1) = Q(:,j+1)-beta(j-1)*Q(:,j-1); 
%         P(:,j+1) = P(:,j+1)-gamma(j-1)*P(:,j-1); 
%     end
%     gamma(j) = norm(Q(:,j+1));
%     Q(:,j+1) = Q(:,j+1)/gamma(j);  
%     beta(j) = P(:,j+1)'*Q(:,j+1);
%     P(:,j+1) = P(:,j+1)/beta(j);
% end

end