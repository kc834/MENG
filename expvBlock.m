% THIS SCRIPT DOES NOT WORK

function [y, err, hump1, hump2] = expvBlock( t, A, u, v, tol, m )

[n, p] = size(v);
if nargin == 3
  tol = 1.0e-7;
  m = min(n,30);
end
if nargin == 4
  m = min(n,30);
end

multiplier = 1;
anorm = norm(A,'inf'); % maximum absolute row sum
mxrej = 10;
btol  = 1.0e-7; 
gamma1 = 0.7; 
delta = 1.2; 
mb    = m; % dimension
t_out   = abs(t); % absolute value of time
t_step = 0;
t_left = 0; 
t_right = t_out; % Added to modify code for range
s_error = 0;
rndoff= anorm*eps;
err_loc = 0;
err_loc_left = 0;
err_loc_right = 0;

k1 = 0;
xm_right = 1/m; 
normv = norm(v); % Euclidean norm 
beta1 = normv;
normuv = u' * v; % Euclidean norm 
% if isreal(v) && isreal(u) && isreal(A)
%     normuv = norm(normuv);
% end
beta2 = normuv; 
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));

% Right side time step
t_new_right = (1/anorm)*((fact*tol)/(4*beta1*anorm))^xm_right;
s = 10^(floor(log10(t_new_right))-1); 
t_new_right = multiplier * ceil(t_new_right/s)*s; 

% Left side time step
t_new_left = (1/anorm)*((fact*tol)/(4*norm(beta2)*anorm))^xm_right;
s = 10^(floor(log10(t_new_left))-1); 
t_new_left = multiplier * ceil(t_new_left/s)*s; 

% Left side time step

sgn = sign(t); 
nstep = 0;
w1 = v;
w2 = u;
hump1 = normv;
hump2 = normuv;
abort = 0;
y = beta2';
while (t_right - t_left) > 0
    % Determine step size
    nstep = nstep + 1;
    t_step_right = min( (t_right-t_left)/2,t_new_right ); % step cannot exceed range
    t_step_left = min( (t_right-t_left)/2,t_new_left );
    % ========== Bi-Lanczos Implementation ==========
    Q = zeros(n,m+p); % Orthonormal basis, n by k+1 array
    P = zeros(n,m+p);
    R = zeros(n,m+p);
    S = zeros(n,m+p);
    T = zeros(m+p+1, m+p+1);

    Q(:,1:p) = v;
    normv = norm(Q(:,1:p)); % Get the norm of each column
    Q(:,1:p) = Q(:,1:p)/normv; % Divide each value by column norm
    P(:,1:p) = u;
    normuv = Q(:,1:p)'*P(:,1:p);
    P(:,1:p) = normv * P(:,1:p)/normuv; % Faster than multiplying by inverse, probably need to multiply back normv
    nblocks = m / p;
    R(:, 1:p) = A'*P(:, 1:p);
    S(:, 1:p) = A*Q(:,1:p);
    for j = 1:nblocks
        low = (j-1)*p + 1; nxtlow = j*p + 1;
        hi = j * p; nxthi = (j+1) * p;
        T(low:hi, low:hi) = P(:, low:hi)' * S(:, low:hi);
        if j ~= nblocks
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
            if norm(R(:, nxtlow:nxthi)) < 1e-7 % happy_breakdown
                k1 = 0;
                mb = j;
                t_step = t_right-t_left;
                break;
            end
            S(:, nxtlow:nxthi) = A*Q(:, nxtlow:nxthi) - Q(:, low:hi)*T(low:hi, nxtlow:nxthi);
            if norm(R(:, nxtlow:nxthi)) < 1e-7 % happy_breakdown
                k1 = 0;
                mb = j;
                t_step = t_right-t_left;
                break;
            end
        end
    end
    % ========== Timestep Calculation ==========
    if (abort)
        break;
    end
    if k1 ~= 0
        T(m+2,m+1) = 1;
        avnorm = norm(R(:,1:p)'*A*Q(:,m+1:m+p));
        awnorm = norm(Q(:,1:p)'*A'*R(:,m+1:m+p));
    end
    ireject = 0;
    right_not_done = 1;
    left_not_done = 1;
    while (ireject <= mxrej)
        mx = mb + k1;
        if k1 == 0 % happy_breakdown, one-step execution
            F = expm(sgn*t_step*T(1:mx,1:mx));
            err_loc = btol;
            break;
        else  % Two-step execution TODO: refactor or use an array
            % Right vector step
            if (right_not_done)
                F_right = expm(sgn*t_step_right*T(1:mx,1:mx));
                phi1_right = abs( F_right(m+1,1) );
                phi2_right = abs( F_right(m+2,1) * avnorm );
                if phi1_right > 10*phi2_right
                    err_loc_right = phi2_right;
                    xm_right = 1/m;
                elseif phi1_right > phi2_right
                    err_loc_right = (phi1_right*phi2_right)/(phi1_right-phi2_right);
                    xm_right = 1/m;
                else
                    err_loc_right = phi1_right;
                    xm_right = 1/(m-1);
                end
            end
            % Left vector step
            if (left_not_done)
                F_left = expm(sgn*t_step_left*T(1:mx,1:mx));
                phi1_left = abs( F_left(m+1,1) );
                phi2_left = abs( F_left(m+2,1) * awnorm );
                if phi1_left > 10*phi2_left
                    err_loc_left = phi2_left;
                    xm_left = 1/m;
                elseif phi1_left > phi2_left
                    err_loc_left = (phi1_left*phi2_left)/(phi1_left-phi2_left);
                    xm_left = 1/m;
                else
                    err_loc_left = phi1_left;
                    xm_left = 1/(m-1);
                end
            end
        end
        if err_loc_right <= delta * t_step_right*tol
            right_not_done = 0;
        else
            t_step_right = gamma1 * t_step_right * (t_step_right*tol/err_loc_right)^xm_right;
            s = 10^(floor(log10(t_step_right))-1);
            t_step_right = multiplier * ceil(t_step_right/s) * s;
        end
        if err_loc_left <= delta * t_step_left*tol
            left_not_done = 0;
        else
            t_step_left = gamma1 * t_step_left * (t_step_left*tol/err_loc_left)^xm_left;
            s = 10^(floor(log10(t_step_left))-1);
            t_step_left = multiplier * ceil(t_step_left/s) * s;
        end
        if ~(right_not_done || left_not_done)
            break;
        else
            if ireject == mxrej
                abort = 1;
                warning('The requested tolerance is too high.');
                break;
            end
            ireject = ireject + 1;
        end
    end
    if (abort)
        break;
    end
    mx = mb + max( 0,k1-1 );
    
    % Equations for calculating result
    if k1 ~= 0 
        % right vector
        w1 = Q(:,1:mx)*(F_right(1:mx,1:p));
        beta1 = norm( w1 );
        hump1 = max(hump1,beta1);
        %disp(t_step_right);
        t_right = t_right - t_step_right;
        
        % left vector uses transpose
        w2 = R(:,1:mx)*F_left(1:p,1:mx)'*beta2'; % TODO: transpose beta2 or nah?
        
        beta2 = w2' * w1;
%         if isreal(v) && isreal(u) && isreal(A)
%             beta2 = norm(beta2);
%         end
        hump2 = max(hump2,beta2);
        %disp(t_step_left);
        t_left = t_left + t_step_left;
        
        % Right side time step
        t_new_right = gamma1 * t_step_right * (t_step_right*tol/err_loc_right)^xm_right;
        s = 10^(floor(log10(t_new_right))-1);
        t_new_right = multiplier * ceil(t_new_right/s) * s;
        
        % Left side time step
        t_new_left = gamma1 * t_step_left * (t_step_left*tol/err_loc_left)^xm_left;
        if (~isreal(t_new_left))
            disp("error");
        end
        s = 10^(floor(log10(t_new_left))-1);
        t_new_left = multiplier * ceil(t_new_left/s) * s;
        
    else % happy_breakdown
        y = (beta2' * F(1:p, 1:p)')';
        % y = beta2' * F(:, 1)' * R' * Q(:, 1);
        t_left = t_right;
    end

    err_loc_right = max(err_loc_right,rndoff);
    err_loc_left = max(err_loc_left,rndoff);
    err_loc = max(err_loc,rndoff);
    s_error = s_error + err_loc_right + err_loc_left + err_loc;
    % disp(t_end - t_now);
end
if (abort)
    y = 100;
end
if k1 ~= 0
    % Answer
    y = w2'*w1;
end
disp('Number of steps, expvBlock');
disp(nstep);
% disp(y);
err = s_error;
hump1 = hump1 / normv;
hump2 = hump2 / normuv;
if (hump1 ~= 1) || (hump2 ~= 1)
    disp(' hump ');
end
% y = w2 * w1;

