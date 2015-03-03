%-------------------------------------------------------------------------------
%
% Read in a matrix in COO format from an ASCII file. The matrix is returned as
% a sparse array A. Entries in A that are unspecified entries are assumed to be
% zero. filename is just what its name says. errorflag is error flag.
%
% The dimensions of A are as with the sparse() function: the max
% of the row numbers and the max of the column numbers.
%
%------------
% Randall Bramley
% Department of Computer Science
% Indiana University, Bloomington
%-----------------------
% Started: Wed 02 Dec 2009, 09:56 AM
% Last Modified: Mon 12 Dec 2011, 02:39 PM
%-------------------------------------------------------------------------------


function [A, errorflag] = readmatrix(filename)

    returnsparse = false;

    %--------------------------------------------------------
    % Should check if filename is of character type, etc.
    %--------------------------------------------------------
    errorflag = false;
    try
        S = load(filename);
    catch wotevuh
        errorflag = true;
        disp('Darn. load failed. Probably the file name is mispelled or permissions wrong');
        messystuff = wotevuh.message;
        disp(sprintf('The error message from load is %s', messystuff));
        A = [];
        return
    end
    
    m = S(1,1);
    n = S(1,2);
    numberentries = S(1, 3);

    rows = S(2:end, 1); 
    cols = S(2:end, 2);
    vals = S(2:end, 3);
    A = sparse(rows, cols, vals);
    if returnsparse
        return
    end
    A = full(A);

return
