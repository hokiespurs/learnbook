function h = plotCovarianceAxes3(ux,uy,uz,C,confidence,dof,varargin)

%% Validate Inputs and Default Values
validateCovariance(C,3);
if nargin==4
    confidence=fcdf(.5,3,inf);
    dof = inf;
elseif nargin==5
    dof=inf;
end
%% Calculate Confidence Ellipse (scaled with multiplier 'c' pg 407 Ghilani)
c = sqrt(2*finv(confidence,3,dof));

%% Calculate Eigenvalues and Eigenvectors
[semiaxesVector,eigValues]=eig(C);
Su = c * sqrt(eigValues(1,1));
Sv = c * sqrt(eigValues(2,2));
Sw = c * sqrt(eigValues(3,3));
%% Plot Vectors
xyz = semiaxesVector(:,1)*Su;
h{1} = plot3([0 xyz(1)]+ux,[0 xyz(2)]+uy,[0 xyz(3)]+uz,varargin{:});
hold on

xyz = semiaxesVector(:,2)*Sv;
h{2} = plot3([0 xyz(1)]+ux,[0 xyz(2)]+uy,[0 xyz(3)]+uz,varargin{:});

xyz = semiaxesVector(:,3)*Sw;
h{3} = plot3([0 xyz(1)]+ux,[0 xyz(2)]+uy,[0 xyz(3)]+uz,varargin{:});
end