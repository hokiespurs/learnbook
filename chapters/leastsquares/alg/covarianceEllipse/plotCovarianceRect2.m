function h = plotCovarianceRect2(ux,uy,C,confidence,dof,varargin)

validateCovariance(C,2);
if nargin==3
    confidence=fcdf(.5,2,inf);
    dof = inf;
elseif nargin==4
    dof=inf;
end
%% Calculate Confidence Ellipse (scaled with multiplier 'c' pg 407 Ghilani)
c = sqrt(2*finv(confidence,2,dof));

sx = c * sqrt(diag(C));
h = plotRect([ux uy],sx,varargin{:});

end