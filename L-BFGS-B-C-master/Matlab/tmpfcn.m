function f = tmpfcn(x,A,b)
disp('f')
f = norm( A*x - b)^2;