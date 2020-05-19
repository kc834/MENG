A = spconvert(load('cematlab_34050_6_2_K.txt'));
UV = readmatrix('uv_34050_6_2_K.txt');
m = 30; tol = 1e-7; t = 3.5;
[n, p] = size(UV);
numerrors = 0;
expv = 0;
expvB = 0;
acc = 0;
for i = 1:p
    for j = 1:p
        if i ~= j
            [error, t1, t2] = compare(UV(:,1),A,t,UV(:, 2),tol,m);
            if (error > 0.001)
                numerrors = numerrors + 1;
            else
                expv = expv + t1;
                expvB = expvB + t2;
            end
        end
    end
    acc = acc + 127;
    a = expv / (acc);
    fprintf('EXPV took, on average. %f seconds\n',a);
    a = expvB / (acc);
    fprintf('EXPVB took, on average. %f seconds\n',a);
    fprintf('Errors. %d\n',numerrors);
end