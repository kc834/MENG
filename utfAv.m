% I preprocessed the data by separating the real and imaginary part with
% a space delimiter. I did a find and replace in a text editor of "+" with 
% " +" and "-" with " -", with exceptions for exponentials, e.g. "e-05".

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
utfAv_SingleArnoldi(u,A,t,v,tol,m);