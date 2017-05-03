function h = plotCovarianceRect3(ux,uy,uz,C,confidence,dof,varargin)

validateCovariance(C,3);
if nargin==4
    confidence=fcdf(.5,2,inf);
    dof = inf;
elseif nargin==5
    dof=inf;
end
%% Calculate Confidence Ellipse (scaled with multiplier 'c' pg 407 Ghilani)
c = sqrt(2*finv(confidence,3,dof));

sx = c * sqrt(diag(C));
h = plotRect([ux uy uz],sx,varargin{:});

end
