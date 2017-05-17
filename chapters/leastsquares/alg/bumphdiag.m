function hx = bumphdiag(x,n)
    % pad a matrix with 0s while shifting n rows at a time over on blkdiag
    % This is useful when doing partial derivatives wrt predictor variables
    %
    % ie, with n=1
    %  1 2 3       1 2 3 0 0 0 0 0 0 0 0 0
    %  4 5 6  ->   0 0 0 4 5 6 0 0 0 0 0 0
    %  7 8 9       0 0 0 0 0 0 7 8 9 0 0 0
    %  2 4 6       0 0 0 0 0 0 0 0 0 2 4 6
    %
    % ie, with n=2
    %  1 2 3       1 2 3 0 0 0
    %  4 5 6  ->   4 5 6 0 0 0
    %  7 8 9       0 0 0 7 8 9
    %  2 4 6       0 0 0 2 4 6
    %
    [M,N]=size(x);
    irow = kron((1:M)',ones(1,N));
    icol = nan(N,M/n);
    icol(:) = 1:numel(icol);
    icol = kron(icol.',ones(n,1));
    hx = sparse(irow,icol,x);
end