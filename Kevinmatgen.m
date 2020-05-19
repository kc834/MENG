%conversion script
%label 1
matzi = spconvert(load('cematlab_imagh_34050_6_1.mtxz'));
matzr = spconvert(load('cematlab_realh_34050_6_1.mtxz'));
shiftroffdiag=1;
ndimo = size(matzr,1);
matz = matzr + 1i * matzi + shiftroffdiag*speye(ndimo);
[row,col,v] = find(matz);
dlmwrite('cematlab_34050_6_1_K.txt',[row col v], 'delimiter', '\t','precision',12);
matz1 = spconvert(load('cematlab_34050_6_1_K.txt'));
if max(max(abs(matz-matz1)))>1e-10
    disp('Success');
end
%label 2
matzi = spconvert(load('cematlab_imagh_34050_6_2.mtxz'));
matzr = spconvert(load('cematlab_realh_34050_6_2.mtxz'));
shiftroffdiag=1;
ndimo = size(matzr,1);
matz = matzr + 1i * matzi + shiftroffdiag*speye(ndimo);
[row,col,v] = find(matz);
dlmwrite('cematlab_34050_6_2_K.txt',[row col v], 'delimiter', '\t','precision',12);
matz1 = spconvert(load('cematlab_34050_6_2_K.txt'));
if max(max(abs(matz-matz1)))>1e-10
    disp('Success');
end