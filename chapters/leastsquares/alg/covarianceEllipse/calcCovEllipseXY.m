function [x,y]=calcCovEllipseXY(ux,uy,C,confidence,dof,ellipseResolution)
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
%% Calculate Confidence Ellipse (scaled with multiplier 'c' pg 407 Ghilani)
c = sqrt(2*finv(confidence,2,dof));

%% Calculate Eigenvalues and Eigenvectors
[semiaxesVector,eigValues]=eig(C);
Su = c * sqrt(eigValues(1,1));
Sv = c * sqrt(eigValues(2,2));

%% Plot Ellipse
az = linspace(-pi,pi,ellipseResolution);
xy = [Su * cos(az); Sv * sin(az)]; %xy coordinate for ellipse
xyRot = semiaxesVector*xy; %use semiaxesVector as rotation matrix
x = xyRot(1,:) + ux;
y = xyRot(2,:) + uy;

end

