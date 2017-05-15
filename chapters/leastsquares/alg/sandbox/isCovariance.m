function isvalid = isCovariance(C,ndim)
% check if a covariance matrix is valid
isvalid = false;
if nargin==2 && size(C,1)~=ndim
    %     warning('Covariance matrix must be %.0fx%.0f',ndim,ndim);
elseif ~issymmetric(C)
    %     warning('Covariance matrix must be Symmetric');
else
    [~,D]=eig(C);
    eigenvals = diag(D);
    if sum(eigenvals<0)>0
        %   warning('Covariance matrix must be positive semi-definite');
    else
        isvalid = true;
    end
end
end