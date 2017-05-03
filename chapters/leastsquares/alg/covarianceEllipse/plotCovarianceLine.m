function h = plotCovarianceLine(ux,uy,C,confidence,dof,ellipseResolution,varargin)
%% Inputs
% u - mean position of ellipse (required)
% C - 2x2 or 3x3 covariance matrix (required)
% confidence - 0:1 probability
% dof - degrees of freedom
% ellipseResolution - scalar or 2x1 vector(ellipsoid) for mesh resolution

%% Validate Inputs and Default Values
validateCovariance(C,2);
if nargin==3
    confidence=fcdf(.5,2,inf);
    dof = inf;
    ellipseResolution = 50;
elseif nargin==4
    dof=inf;
    ellipseResolution = 50;
elseif nargin==5
    ellipseResolution = 50;
end
[x,y]=calcCovEllipseXY(ux,uy,C,confidence,dof,ellipseResolution);

h = plot(x,y,varargin{:});

end

