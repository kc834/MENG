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
% High-level algorithm from PhD Thesis of Axel Facius, July 2000, Universitat Karlsruhe 
% Rational approximation derived with insight from "Efficient Solution of Parabolic 
%          Equations by Krylov Approximation Methods" by E. Gallopoulos and Y. Saad
% Insight for preprocessing from "Analysis of Some Krylov Subspace Approximations to 
%                                        the Matrix Exponential Operator" by Y. Saad


function out = bilanczos(A,v,u,k)

n = length(A); % Max dim
Q = zeros(n,k+1); % Orthonormal basis, n by k+1 array
R = zeros(n,k+1);
T = zeros(k+1, k+1);

Q(:,1) = v; % Arbitrary vector with norm 1
normv = norm(Q(:,1)); % beta1 = normv, beta2 = normu * normv
Q(:,1) = Q(:,1)/normv;
R(:,1) = u;
normuv = norm(Q(:,1)'*R(:,1));
R(:,1) = R(:,1)/normuv;
for j = 1:k
    Q(:,j+1) = A*Q(:,j); % Move on to next vector in Krylov subspace
    R(:,j+1) = A'*R(:,j);
    T(j,j) = R(:,j)'*Q(:,j+1);
    Q(:,j+1) = Q(:,j+1)-T(j,j)*Q(:,j); 
    R(:,j+1) = R(:,j+1)-T(j,j)*R(:,j); 
    if j > 1 
        Q(:,j+1) = Q(:,j+1)-T(j-1, j)*Q(:,j-1); 
        R(:,j+1) = R(:,j+1)-T(j, j-1)*R(:,j-1); 
    end
    
    T(j+1,j) = norm(Q(:,j+1));
    if norm(T(j+1,j)) < 1e-7
        break;
    end
    Q(:,j+1) = Q(:,j+1)/T(j+1,j);
    T(j,j+1) = R(:,j+1)'*Q(:,j+1);
    if norm(T(j,j+1)) < 1e-7
        break;
    end
    R(:,j+1) = R(:,j+1)/T(j,j+1);
end

T
% disp(Q);
% disp(R);
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
% disp(R(:, 1:k+1)*T(1:k, 1:k+1)');
% disp(A'*R(:, 1:k));
% eqR = all(ismembertol(R(:, 1:k+1)*T(1:k, 1:k+1)',A'*R(:, 1:k), 1e-7), 'all');
% if eqR
%     disp("good eqR");
% else
%     disp("bad eqR");
% end

out = normuv * normv * R(:, 1)' * Q * T(:, 1);

disp(Q(:, 1:k)'*R(:, 1:k));
% disp(R(:, :)'*Q(:, :));
% eqI = all(ismembertol(Q(:, 1:k)'*R(:, 1:k),R(:, 1:k)'*Q(:, 1:k), 1e-7), 'all');
% if eqI
%     disp("good eqI");
% else
%     disp("bad eqI");
% end

%  ======================== OLD CODE =======================

% n = length(A); % Max dim
% Q = zeros(n,k+1); % Orthonormal basis, n by k+1 array
% R = zeros(n,k+1);
% alpha = zeros(k,1); 
% beta = zeros(k,1);
% gamma = zeros(k,1);
% 
% Q(:,1) = b/norm(b); % Arbitrary vector with norm 1
% Q(:,1) = Q(:,1)/norm(Q(:,1));
% R(:,1) = Q(:,1);
% for j = 1:k 
%     Q(:,j+1) = A*Q(:,j); % Move on to next vector in Krylov subspace
%     R(:,j+1) = A'*R(:,j);
%     alpha(j) = R(:,j)'*Q(:,j+1);
%     Q(:,j+1) = Q(:,j+1)-alpha(j)*Q(:,j); 
%     R(:,j+1) = R(:,j+1)-alpha(j)*R(:,j); 
%     if j > 1 
%         Q(:,j+1) = Q(:,j+1)-beta(j-1)*Q(:,j-1); 
%         R(:,j+1) = R(:,j+1)-gamma(j-1)*R(:,j-1); 
%     end
%     gamma(j) = norm(Q(:,j+1));
%     Q(:,j+1) = Q(:,j+1)/gamma(j);  
%     beta(j) = R(:,j+1)'*Q(:,j+1);
%     R(:,j+1) = R(:,j+1)/beta(j);
% end

end