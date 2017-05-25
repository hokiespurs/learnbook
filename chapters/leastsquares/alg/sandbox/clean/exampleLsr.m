function exampleLsr
%% Fit Different Models to a set of unweighted 2D data
% raw data
x = [0 1 2 3 4 5 6]';
y = [1 4 3 7 6 20 19]';

% Linear Trend  y = mx+b
modelfunLinear = @(b,x) b(1)*x + b(2);
betacoefLinear = lsr(x,y,modelfunLinear);

% 2nd Order Polynomial   y = ax^2+bx+c
modelfunPoly2 = @(b,x) b(1)*x.^2 + b(2)*x +b(3);
betacoefPoly2 = lsr(x,y,modelfunPoly2);

% 4th Order Polynomial  y = ax^4+bx^3+cx^2+dx+e
modelfunPoly4 = @(b,x) b(1)*x.^4 + b(2)*x.^3 + b(3)*x.^2 + b(4)*x.^1 + b(5);
betacoefPoly4 = lsr(x,y,modelfunPoly4);

% Exponential (NonLinear) y = ae^(-bx)
modelfunExp = @(b,x) b(1)*exp(b(2)*x);
betacoef0 = [3 .5]';
betacoefExponential = lsr(x,y,modelfunExp,betacoef0);

%plot
xi = 0:0.1:6;
f = figure(1);clf
plot(x,y,'k.','markersize',25);
hold on
plot(xi,modelfunLinear(betacoefLinear,xi),'r','linewidth',2);
plot(xi,modelfunPoly2(betacoefPoly2,xi),'g','linewidth',2);
plot(xi,modelfunPoly4(betacoefPoly4,xi),'c','linewidth',2);
plot(xi,modelfunExp(betacoefExponential,xi),'b','linewidth',2);
legend({'Raw Data','$y = ax + b$','$y = ax^2+bx+c$','$y = ax^4+bx^3+cx^2+dx+e$','$y = ae^{bx}$'},...
    'interpreter','latex','fontsize',14,'location','best')
saveas(f,'basicUsage.png');
%% Fit Model to 3D Plane (Unweighted)
% beta = [a b c d]  % ax + by + c = z
modelfun3Dplane = @(b,x)(b(1)*x(:,1) + b(2)*x(:,2) + b(3));

%generate 100 data points with X=[-5 5] Y=[-5 5]
rng(1);
truebeta = [-.1 -.5 2]';
xpts = (rand(100,1)-0.5)*10;
ypts = (rand(100,1)-0.5)*10;
zpts = modelfun3Dplane(truebeta,[xpts ypts]) + randn(100,1)*.5;

% do least squares
x = [xpts ypts];
y = zpts;

betacoefPlane = lsr(x,y,modelfun3Dplane);

%plot
xi = [-5 -5 5 5 -5]';
yi = [-5 5 5 -5 -5]';
f = figure(2);clf
subplot 121
plot3(xpts,ypts,zpts,'k.');
hold on
% fill3(xi,yi,modelfun3Dplane(truebeta,[xi yi]),'g');alpha 0.3
fill3(xi,yi,modelfun3Dplane(betacoefPlane,[xi yi]),'b','edgecolor','b','linewidth',3);
alpha 0.3
xlabel('X coordinate','interpreter','latex','fontsize',14);
ylabel('Y coordinate','interpreter','latex','fontsize',14);
zlabel('Z coordinate','interpreter','latex','fontsize',14);
subplot 122
plot3(xpts,ypts,zpts,'k.');
hold on
fill3(xi,yi,modelfun3Dplane(betacoefPlane,[xi yi]),'b','edgecolor','b','linewidth',3);
alpha 0.3
xlabel('X coordinate','interpreter','latex','fontsize',14);
ylabel('Y coordinate','interpreter','latex','fontsize',14);
zlabel('Z coordinate','interpreter','latex','fontsize',14);
view(135,35)

%% Sin Wave With Known Period (Nonlinear Observation Equations)
% Nonlinear

% Linear

%% Different ways to weight equations


%% 2D Conformal Transformation with covariances (Linear 2 Equations per Observation)
x_coord2 = [1 2 3]; y_coord2 = [0 5 1];  % raw data 'to'
x_coord1 = [6 1 8]; y_coord1 = [3 12 8]; % raw data 'from'

Sc = [0.5 0.3 0 0 0 0;
    0.3 0.5 0 0 0 0;
    0 0 0.4 0.1 0 0;
    0 0 0.1 0.2 0 0;
    0 0 0 0 0.7 -0.4;
    0 0 0 0 -0.4 0.4]; %variance-covariance of data2

modelfun = @conformal2d;
y = nan(2*numel(x_coord1),1);
y(1:2:end)=x_coord2;
y(2:2:end)=y_coord2;

x = [x_coord1' y_coord1'];

betacoef2DConformal = lsr(x,y,modelfun,'Weights',Sc,'verbose',true);

%% Unweighted 3D Conformal Transformation (Nonlinear 3 Equations per observation)
modelConformal = @(b,x) conformal3dfun(b(1),b(2),b(3),b(4),b(5),b(6),b(7),x);
% Here the scale is fixed == 1, and not solved as a beta coefficient
modelConformalFixScale = @(b,x) conformal3dfun(1,b(1),b(2),b(3),b(4),b(5),b(6),x);

%generate data
xpts = (rand(10,1)-0.5)*100;
ypts = (rand(10,1)-0.5)*100;
zpts = (rand(10,1)-0.5)*100;
x = [xpts ypts zpts];

truebeta = [1 pi/2 pi pi/4 2 3 4]';
XYZ = modelConformal(truebeta,x) + randn(3*numel(xpts),1);
Xpts = XYZ(1:3:end);
Ypts = XYZ(2:3:end);
Zpts = XYZ(3:3:end);

% do least squares
y = [Xpts Ypts Zpts]';
y = y(:);
betacoef0 = truebeta;
betacoef3Dconformal = lsr(x,y,modelConformal,betacoef0,'verbose',true);
betacoef3Dconformal2 = lsr(x,y,modelConformalFixScale,betacoef0(2:end),'verbose',true);

%% Linear Line with Total Least Squares

%% Robust Least Squares for Line with outliers

%% Plot/Report Some Data from Outputs (demonstrate Residuals, covariance, etc)

%% Chi2 Test for linear line, Dont Scale Covariance
% noscale, chi2alpha

%% Analytical Jacobian and Bfunction

%% Weight betacoef0 With input covariance

%% DerivStep

end

function y = conformal2d(beta,x)
% 2D conformal coordinate transformation (beta = [a;b;c;d], x = [xc(:) yc(:)])
nObservations = size(x,1);
y = nan(nObservations*2,1);
y(1:2:end) = beta(1)*x(:,1) - beta(2)*x(:,2) + beta(3);
y(2:2:end) = beta(2)*x(:,1) + beta(1)*x(:,2) + beta(4);
end

function y = conformal3dfun(S,omega,phi,kappa,Tx,Ty,Tz,x)
[X,Y,Z] = conformal3d(S,omega,phi,kappa,Tx,Ty,Tz,x(:,1),x(:,2),x(:,3));
y = [X Y Z]';
y = y(:);
end

function [X,Y,Z] = conformal3d(S,omega,phi,kappa,Tx,Ty,Tz,x,y,z)
Rx = [1 0 0; 0 cos(omega) sin(omega); 0 -sin(omega) cos(omega)];
Ry = [cos(phi) 0 -sin(phi); 0 1 0;sin(phi) 0 cos(phi)];
Rz = [cos(kappa) sin(kappa) 0; -sin(kappa) cos(kappa) 0; 0 0 1];
R = Rx*Ry*Rz;

XYZ = S * R * [x(:)';y(:)';z(:)'] + [Tx; Ty; Tz];

X = XYZ(1,:)';
Y = XYZ(2,:)';
Z = XYZ(3,:)';
end

function [x,y,S,isoutlier]=generate2DSyntheticData(RANDSEED,TRUEBETA,...
    NSAMPLES,NOUTLIERS,OUTLIERMAG,OUTLIERSTD,XRANGE,NOISEMAGNITUDEX,...
    NOISEMAGNITUDEY,NOISEMAGNITUDEXY,MODELFUN)
%% Generate a random 2D dataset with covariances

% seed the random number generator
rng(RANDSEED);
%% Add Outliers
badX = rand(NOUTLIERS,1)*diff(XRANGE)+XRANGE(1);
badY = MODELFUN(TRUEBETA,badX) + (round(rand(NOUTLIERS,1))*2-1) * OUTLIERMAG + randn(NOUTLIERS,1)*OUTLIERSTD;

%% Generate NonLinear Data
for i=1:NSAMPLES
    Sx(i)=rand(1)*diff(NOISEMAGNITUDEX) + NOISEMAGNITUDEX(1);
    Sy(i)=rand(1)*diff(NOISEMAGNITUDEY) + NOISEMAGNITUDEY(1);
    Sxy(i)=rand(1)*NOISEMAGNITUDEXY;
    S{i} = [Sx(i) Sxy(i);Sxy(i) Sy(i)];
    while ~isCovariance(S{i},2)
        Sxy(i)=randn(1)*NOISEMAGNITUDEXY;
        S{i} = [Sx(i) Sxy(i);Sxy(i) Sy(i)];
    end
    xyError = mvnrnd([0 0],S{i},1);
    xError(i) = xyError(1);
    yError(i) = xyError(2);
end
%    trueX = rand(NSAMPLES,1)*diff(XRANGE) + XRANGE(1);
trueX = linspace(XRANGE(1),XRANGE(2),NSAMPLES)';
trueY = MODELFUN(TRUEBETA,trueX);
valX = trueX + xError';
valY = trueY + yError';

% add outlines
x = [valX; badX];
y = [valY; badY];

isoutlier = true(size(x));
isoutlier(1:numel(valX))=false;
end