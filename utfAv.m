% I preprocessed the data by separating the real and imaginary part with
% a space delimiter. I did a find and replace in a text editor of "+" with 
% " +" and "-" with " -", with exceptions for exponentials, e.g. "e-05".
%
% If you get the error "Unable to find or open 'A_matrix.dat'", run
% git checkout b93e4cd4ab5522dfca211f3cbebc88b1566676fb -- A_matrix.dat

u = load('u.dat');
u = u(:,1) + 1i * u(:,2);
v = load('v.dat');
v = v(:,1) + 1i * v(:,2);
% u = eye(4871,1);
% v = u;
A = importA('A_matrix.dat');
A = table2array(A);
A = spconvert(A);
t = 3.5;
tol=1e-7; 
m=30;
utfAv_SingleArnoldi(u,A,t,v,tol,30);