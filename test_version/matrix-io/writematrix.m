%------------------------------------------------------------------------------------
%
% Write out a matrix in COO format to an ASCII file. The matrix can be either sparse
% or dense, If something goes wrong, an empty matrix is returned. The error flag is
% a Boolean, so it carries only one bit of knowledge.
%
% The dimensions of A are as with the sparse() function: the max of the row numbers
% and the max of the column numbers. The values are written out to 17 significant
% digits, way overkill for most applications but needed for reproducible results
% in double precision floating point computations.
%
%
%-------------
% Randall Bramley
% Department of Computer Science
% Indiana University, Bloomington
%--------------------------
% Started: Fri 09 Dec 2011, 12:24 PM
% Last Modified: Mon 12 Dec 2011, 02:37 PM
%------------------------------------------------------------------------------------

function errorflag = writematrix(A, filename)

    %----------------------------------------------------------------------
    % Output the matrix as sparse, ie., without entries exactly equal to 0
    %----------------------------------------------------------------------
    writesparse = false;

    %-----------------------------------------------------
    % Write a header line giving three integers m, n, lda
    %-----------------------------------------------------
    withheaderline = true;

    errorflag = false;

    [fid, messy] = fopen(filename, 'wt');
    if (fid == -1)
        disp('Big trouble in little writematrix(). Could not open the file');
        disp(sprintf('named %s', filename));
        disp(sprintf('Message from fopen: %s', messy));
        disp(sprintf('Setting A to empty matrix just to keep things going. Check file permissions.'));
        A = [];
        return
    end

    if writesparse
        [rows, cols, vals] = find(A);
        nzeros = length(rows);
        if withheaderline
            fprintf(fid, '%d   %d   %d\n', max(rows),  max(cols), max(rows));
        end
        for k = 1:nzeros
            fprintf(fid, '%d   %d   %27.17e \n', rows(k), cols(k), vals(k));
        end
    else
        A = full(A);
        [m, n] = size(A);
        if withheaderline
            fprintf(fid, '%d   %d   %d\n', m, n, m*n);
        end
        for i = 1:n
            for j = 1:m
                % disp(sprintf('A(%d, %d) = %27.17e', j, i, A(j,i))); 
                fprintf(fid, '%d   %d   %27.17e \n', j, i, A(j, i));
            end
        end
    end

    fclose(fid);

return
