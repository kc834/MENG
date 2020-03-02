% calculates u^\dagger exp(-At) v via single vector Arnoldi

% single vector Arnoldi: just do Arnoldi on A and v, forget u, then bring in
% u only at the very end to just take a dot product

% biLanczos version will look like utfAv_BiLanczos(u,A,t,v,tol,m)
% tol: tolerance, m: number of BiLanczos steps

% make sure you are passing sparse matrices as A: these matrices are annoyingly huge
% if A is stored in a file A_matrix.dat, you want to run the following
% command to load A into MATLAB: A = spconvert(load('A_matrix.dat'));

% if you are confused about what tol and m to use, go for tol=1e-7 and m=30.
% But one must explore the effect of tol and m=30 on the output of 
% utfAv_SingleArnoldi in order to compare biLanczos and Arnoldi.
% use t=3.5 if you can't think of any other t value

function y=utfAv_SingleArnoldi(u,A,t,v,tol,m)
% u, A, v are complex matrices/vectors!

% here f(A) is exp(-At)
% not important for the math, but just putting it out here: t in 1/Gauss units, A is in Gauss

[w, err, hump] = expv( -t, A, v, tol, m );
% err and hump have some use, maybe when you are doing an error analysis
% they matter

% return the answer after taking a dot product with u
y = u'*w; % should be a scalar

end