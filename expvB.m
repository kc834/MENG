function [y, err, hump1, hump2] = expvB( t, A, u, v, tol, m )

% preprocessing
[n,n] = size(A);
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

k1 = 2;
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
    n = length(A); 
    Q = zeros(n,m+1); % Right hand vectors
    R = zeros(n,m+1); % Left hand vectors
    T = zeros(m+2, m+2); % Tridiagonal matrix with zero-padding
    Q(:,1) = w1/beta1;
    R(:,1) = beta1 * w2/(beta2');
    for j = 1:m 
        Q(:,j+1) = A*Q(:,j);
        R(:,j+1) = A'*R(:,j);
        T(j,j) = R(:,j)'*Q(:,j+1);
        Q(:,j+1) = Q(:,j+1)-T(j,j)*Q(:,j); 
        R(:,j+1) = R(:,j+1)-T(j,j)'*R(:,j); 
        if j > 1 
            Q(:,j+1) = Q(:,j+1)-T(j-1, j)*Q(:,j-1); 
            R(:,j+1) = R(:,j+1)-T(j, j-1)'*R(:,j-1); 
        end
        T(j+1,j) = norm(Q(:,j+1));
        if norm(T(j+1,j)) < 1e-7 % happy_breakdown, catch divide by 0
            k1 = 0;
            mb = j;
            t_step = t_right-t_left;
            break;
        end
        Q(:,j+1) = Q(:,j+1)/T(j+1,j);
        T(j,j+1) = R(:,j+1)'*Q(:,j+1);
        if (norm(R(:,j+1)) < 1e-7) % happy breakdown
            k1 = 0;
            mb = j;
            t_step = t_right-t_left;
            break;
        end
        if (norm(T(j,j+1)) < 1e-7) % serious breakdown
            abort = 1;
            warning('SERIOUS BREAKDOWN');
            break;
        end
        R(:,j+1) = R(:,j+1)/(T(j,j+1)');
    end
    % ========== Timestep Calculation ==========
    if (abort)
        break;
    end
    if k1 ~= 0
        T(m+2,m+1) = 1;
        avnorm = norm(R(:,1)'*A*Q(:,m+1));
        awnorm = norm(Q(:,1)'*A'*R(:,m+1));
    end
    ireject = 0;
    right_not_done = 1;
    left_not_done = 1;
    while (ireject <= mxrej)
        mx = mb + k1;
        % F = T;
        % y = beta2 * R(:, 1)' * Q(:, :) * F(1:mx-1, 1);
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
        w1 = Q(:,1:mx)*(F_right(1:mx,1));
        beta1 = norm( w1 );
        hump1 = max(hump1,beta1);
        %disp(t_step_right);
        t_right = t_right - t_step_right;
        
        % left vector uses transpose
        w2 = beta2'*R(:,1:mx)*F_left(1,1:mx)';
        
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
        y = (beta2' * F(1, 1)')';
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
% disp('Number of steps, expvB');
% disp(nstep);
% disp(y);
err = s_error;
hump1 = hump1 / normv;
hump2 = hump2 / normuv;
% if (hump1 ~= 1) || (hump2 ~= 1)
%     disp(' hump ');
% end
% y = w2 * w1;

