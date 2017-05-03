function hx = bumphdiag(x,n)
[M,N]=size(x);
irow = kron((1:M)',ones(1,N));
icol = nan(N,M/n);
icol(:) = 1:numel(icol);
icol = kron(icol.',ones(n,1));
hx = sparse(irow,icol,x);
end