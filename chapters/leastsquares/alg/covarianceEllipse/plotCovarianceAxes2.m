function h = plotCovarianceAxes2(ux,uy,C,confidence,dof,varargin)

%% Validate Inputs and Default Values
validateCovariance(C,2);
if nargin==3
    confidence=fcdf(.5,2,inf);
    dof = inf;
elseif nargin==4
    dof=inf;
end
%% Calculate Confidence Ellipse (scaled with multiplier 'c' pg 407 Ghilani)
c = sqrt(2*finv(confidence,2,dof));

%% Calculate Eigenvalues and Eigenvectors
[semiaxesVector,eigValues]=eig(C);
Su = c * sqrt(eigValues(1,1));
Sv = c * sqrt(eigValues(2,2));

%% Plot Vectors
xy = semiaxesVector(:,1)*Su;
h{1} = plot([0 xy(1)]+ux,[0 xy(2)]+uy,varargin{:});
hold on
xy = semiaxesVector(:,2)*Sv;
h{2} = plot([0 xy(1)]+ux,[0 xy(2)]+uy,varargin{:});

end