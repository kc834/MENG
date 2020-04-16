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

anorm = norm(A,'inf'); % maximum absolute row sum
mxrej = 10;
btol  = 1.0e-7; 
gamma1 = 0.9; 
delta = 1.2; 
mb    = m; % dimension
t_out   = abs(t); % absolute value of time
nstep = 0; 
t_new   = 0;
t_now = 0; 
t_end = t_out; % Added to modify code for range
s_error = 0;
rndoff= anorm*eps;

k1 = 2;
xm = 1/m; 
normv = norm(v); % Euclidean norm 
beta1 = normv; 
normu = norm(u); % Euclidean norm 
beta2 = normu; 
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta1*anorm))^xm;
s = 10^(floor(log10(t_new))-1); 
t_new = ceil(t_new/s)*s; 
sgn = sign(t); 
nstep = 0;

w1 = v;
w2 = u;
hump1 = normv;
hump2 = normu;
while (t_end - t_now) > 0 % while range is high
  nstep = nstep + 1;
  t_step = min( (t_end-t_now)/2,t_new ); % step cannot exceed range
  
  n = length(A); % Max dim
  Q = zeros(n,m+1); % Orthonormal basis, n by m+1 array
  R = zeros(n,m+1);
  T = zeros(m+2,m+2); % Added extra column since matrix has to be square
  alpha = zeros(m+1,1); 
  beta = zeros(m+1,1);
  gamma = zeros(m+1,1);

  Q(:,1) = v; % Arbitrary vector
  Q(:,1) = Q(:,1)/norm(Q(:,1));
  R(:,1) = Q(:,1);
  for j = 1:m % BiLanczos
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
      beta(j) = R(:,j+1)'*Q(:,j+1);
      Q(:,j+1) = Q(:,j+1)/gamma(j);  
      R(:,j+1) = R(:,j+1)/beta(j);
      T(j,j) = alpha(j);
      T(j+1,j) = gamma(j);
%       if j ~= m
      T(j, j+1) = beta(j);
%       end
      if (gamma(j) < btol) || (beta(j) < btol)
        k1 = 0;
        mb = j;
        t_step = t_out-t_now;
        break;
     end
  end
  if k1 ~= 0
     T(m+2,m+1) = 1;
     avnorm = norm(A*Q(:,m+1)); 
  end
  ireject = 0;
  while ireject <= mxrej % Finite number of rejections
     mx = mb + k1;
     F = expm(sgn*t_step*T(1:mx,1:mx));
     if k1 == 0
	      err_loc = btol; 
        break;
     else
        phi1 = abs( beta1*F(m+1,1) );
        phi2 = abs( beta1*F(m+2,1) * avnorm );
        if phi1 > 10*phi2
           err_loc = phi2;
           xm = 1/m;
        elseif phi1 > phi2
           err_loc = (phi1*phi2)/(phi1-phi2);
           xm = 1/m;
        else
           err_loc = phi1;
           xm = 1/(m-1);
        end
     end
     if err_loc <= delta * t_step*tol
        break;
     else
        t_step = gamma1 * t_step * (t_step*tol/err_loc)^xm;
        s = 10^(floor(log10(t_step))-1);
        t_step = ceil(t_step/s) * s;
        if ireject == mxrej
           error('The requested tolerance is too high.');
        end
        ireject = ireject + 1;
     end
  end
  mx = mb + max( 0,k1-1 );
  w1 = Q(:,1:mx)*(beta1*F(1:mx,1)); % right vector
  beta1 = norm( w1 ); % distinct beta
  hump1 = max(hump1,beta1);
  
  w2 = R(:,1:mx)*(beta2*F(1,1:mx)'); % left vector uses transpose
  beta2 = norm( w2 );
  hump2 = max(hump2,beta2);

  t_now = t_now + t_step;
  t_end = t_end - t_step;
  t_new = gamma1 * t_step * (t_step*tol/err_loc)^xm;
  s = 10^(floor(log10(t_new))-1); 
  t_new = ceil(t_new/s) * s;

  err_loc = max(err_loc,rndoff);
  s_error = s_error + err_loc;
  % disp(t_end - t_now);
end
err = s_error;
hump1 = hump1 / normv;
hump2 = hump2 / normu;
y = w1'*w2;
