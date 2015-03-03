
%----------------------------------
% First set up and output a matrix
%----------------------------------
m = 17;
n = 4;
A = randn(m, n);

filename = sprintf('A_%d_%d', m, n);
errorflag = writematrix(sparse(A), filename);
if (errorflag)
    disp(sprintf('Something wrong in writing the file %s; error flag = ', ...
                  filename, errorflag)) 
end

%---------------------------------
% Now try reading the matrix in: 
%---------------------------------
[B, errorflag] = readmatrix(filename);
[m, n] = size(B);
disp(sprintf('\nGotcha matrix right here. It is %d x %d', m, n));
if (m*n < 64)
    disp(sprintf('The values in B are:'));
    disp(B);
end

error_metric = norm(A-B, 1);
disp(sprintf('1-norm of difference of A and B = %g', error_metric))
