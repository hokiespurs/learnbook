function h = plotCovarianceSurf(ux,uy,uz,C,confidence,dof,ellipseResolution,varargin)

if nargin==4
    confidence=fcdf(.5,3,inf);
    dof = inf;
    ellipseResolution = [50 50];
elseif nargin==5
    dof=inf;
    ellipseResolution = [50 50];
elseif nargin==6
    ellipseResolution = [50 50];
end
validateCovariance(C,3);

[x,y,z]=calcCovEllipsoidXYZ(ux,uy,uz,C,confidence,dof,ellipseResolution);

h = surf(x,y,z,varargin{:});

end



