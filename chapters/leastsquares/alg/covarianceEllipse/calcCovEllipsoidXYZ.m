function [x,y,z]=calcCovEllipsoidXYZ(ux,uy,uz,C,confidence,dof,ellipseResolution)
%% Validate Inputs and Default Values
validateCovariance(C,3);
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
%% Calculate Confidence Ellipse (scaled with multiplier 'c' pg 407 Ghilani)
c = sqrt(2*finv(confidence,3,dof));

%% Calculate Eigenvalues and Eigenvectors
[semiaxesVector,eigValues]=eig(C);
Su = c * sqrt(eigValues(1,1));
Sv = c * sqrt(eigValues(2,2));
Sw = c * sqrt(eigValues(3,3));

%% Plot Ellipse
az = linspace(-pi,pi,ellipseResolution(1));
elev = linspace(-pi/2,pi/2,ellipseResolution(2));
[elevgrid,azgrid]=meshgrid(elev,az);

xyz = [Su .* cos(elevgrid(:)) .* cos(azgrid(:)) ...
       Sv .* cos(elevgrid(:)) .* sin(azgrid(:)) ...
       Sw .* sin(elevgrid(:))]; %xyz coordinate for ellipse
   
xyRot = semiaxesVector*xyz'; %use semiaxesVector as rotation matrix
x = reshape(xyRot(1,:) + ux,size(elevgrid));
y = reshape(xyRot(2,:) + uy,size(elevgrid));
z = reshape(xyRot(3,:) + uz,size(elevgrid));

end