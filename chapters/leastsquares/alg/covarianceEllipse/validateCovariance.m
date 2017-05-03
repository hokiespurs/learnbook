function isgood = validateCovariance(C,ndim)

if ~issymmetric(C)
    error('Covariance matrix must be Symmetric');
end

if size(C,1)~=ndim
    error('Covariance matrix must be %.0fx%.0f',ndim,ndim);
end

[~,D]=eig(C);

eigenvals = diag(D);

if sum(eigenvals<0)>0
        error('Covariance matrix must be positive semi-definite');
end
isgood = 1;
end
