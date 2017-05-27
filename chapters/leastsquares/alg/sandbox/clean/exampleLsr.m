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
f1 = figure(1);clf
plot(x,y,'k.','markersize',25);
hold on
plot(xi,modelfunLinear(betacoefLinear,xi),'r','linewidth',2);
plot(xi,modelfunPoly2(betacoefPoly2,xi),'g','linewidth',2);
plot(xi,modelfunPoly4(betacoefPoly4,xi),'c','linewidth',2);
plot(xi,modelfunExp(betacoefExponential,xi),'b','linewidth',2);
legend({'Raw Data','$y = ax + b$','$y = ax^2+bx+c$','$y = ax^4+bx^3+cx^2+dx+e$','$y = ae^{bx}$'},...
    'interpreter','latex','fontsize',14,'location','best')
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
f2 = figure(2);clf
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
omega = 2*pi/10; %10 second period wave
modelfunSinNonLinear = @(b,x) b(1)*sin(omega*x + b(2));
modelfunSinLinear = @(b,x) b(1)*sin(omega*x)+b(2)*cos(omega*x);

% generate sample data
t = 0:0.1:100;
truebeta = [3 pi/3];
z = modelfunSinNonLinear(truebeta,t)+randn(size(t))*0.5;

% Nonlinear
x = t';
y = z';
betacoefSinNonLinear = lsr(x,y,modelfunSinNonLinear,truebeta,'verbose',true);

% Linear
betacoefSinLinear = lsr(x,y,modelfunSinLinear,'verbose',true);
Amp = sqrt((betacoefSinLinear(1))^2+(betacoefSinLinear(2))^2);
Phi = atan2(betacoefSinLinear(2),betacoefSinLinear(1));

% Plot
f3 = figure(3);clf;hold on;
plot(t,modelfunSinLinear(betacoefSinLinear,t),'b-','linewidth',2);
plot(t,z,'k.');
legend({'Raw Data','$y = A\sin(\omega t) + B\cos(\omega t) $'},'interpreter','latex','fontsize',14)
xlabel('Time (s)','interpreter','latex','fontsize',14);
ylabel('Amplitude','interpreter','latex','fontsize',14);
title('Linear Regression of Sin Wave With Known Period','interpreter','latex','fontsize',16)
%% Different ways to weight equations
x = [0 1 2 3 4 5 6]';
y = [1 4 3 7 6 15 14]';
sigmaY = [2 1 3 4 3 2 1];

modelfunLinear = @(b,x) b(1)*x + b(2);
% No weights
betacoefNoWeight = lsr(x,y,modelfunLinear);
% Weights as vector
weightVector = 1./(sigmaY.^2);
betacoefWeightVector = lsr(x,y,modelfunLinear,'Weights',weightVector);
% Weights as Covariance Matrix
covarianceMatrix = diag(sigmaY.^2);
betacoefCovMatrix = lsr(x,y,modelfunLinear,'Weights',covarianceMatrix);

% Plot
xi = 0:0.1:6;
f4 = figure(4);clf
scatter(x,y,weightVector*200,ones(size(weightVector)),'filled');
hold on
plot(xi,modelfunLinear(betacoefNoWeight,xi),'r','linewidth',2);
plot(xi,modelfunLinear(betacoefWeightVector,xi),'b','linewidth',2);
plot(xi,modelfunLinear(betacoefCovMatrix,xi),'b','linewidth',2);
legend({'Raw Data','No Weights','Weighted'},...
    'interpreter','latex','fontsize',14,'location','best')

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
xpts = [10 20 60 40 85];
ypts = [0 15 23 25 40];
covxy = blkdiag([45 -30;-30 30],[20 -10;-10 70],[80 4;4 4],[40 -13;-13 60],[30 -25;-25 30]);

modelfunLinear = @(b,x) b(1)*x + b(2);
modelfunLinearTLS = @(b,x) b(1)*x(:,1) + b(2) -x(:,2);

% ordinary least squares
x = xpts';
y = ypts';
betacoefOrdinary = lsr(x,y,modelfunLinear);
% total least squares
x = [xpts' ypts'];
y = zeros(5,1);
beteacoef0 = betacoefOrdinary;
betacoefTLS = lsr(x,y,modelfunLinearTLS,beteacoef0,'type','tls','Weights',covxy);

% plot
f5 = figure(5);clf;hold on;
xi = 0:100;
plot(xpts,ypts,'k.','markersize',25);
plot(xi,modelfunLinear(betacoefOrdinary,xi),'b-','linewidth',2);
plot(xi,modelfunLinear(betacoefTLS,xi),'r-','linewidth',2);
legend({'Raw Data','Ordinary Least Squares',...
    'Total Least Squares'},...
    'interpreter','latex','fontsize',14,'location','best')
%% Robust Least Squares for Line with outliers
% data has one outlier (8,12)
x = [1 2 3 4 5 6 7 8 9 10 8 9]';
y = [1.2 2 3.1 4 5.3 6 7.4 8 9.5 10 2 3]';
modelfunLinear = @(b,x) b(1)*x + b(2);

% Ordinary Least Squares
betacoefOrdinary = lsr(x,y,modelfunLinear,'verbose',true);

% Robust Least Squares
betacoefRobust = lsr(x,y,modelfunLinear,'RobustWgtFun','andrews','verbose',true);

%plot
f6 = figure(6);clf;hold on;
plot(x(1:10),y(1:10),'k.','markersize',25);
plot(x(11:12),y(11:12),'mo','MarkerFaceColor','k','linewidth',2,'markersize',10);
plot(x,modelfunLinear(betacoefOrdinary,x),'b-','linewidth',2);
plot(x,modelfunLinear(betacoefRobust,x),'r-','linewidth',2);
legend({'Raw Data(inliers)','Raw Data(outliers)','Ordinary Least Squares',...
    'Robust Least Squares (Andrews)'},...
    'interpreter','latex','fontsize',14,'location','best')

%% Plot/Report Some Data from Outputs (demonstrate Residuals, covariance, etc)

%% Do TPU with function handles
% Compute a Rotation Matrix using data, then use CovB to propagate error
modelConformal = @(b,x) conformal3dfun(b(1),b(2),b(3),b(4),b(5),b(6),b(7),x);

%generate data
rng(1);
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

[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelConformal,...
    betacoef0,'verbose',true);
% propagate error for a new point
xyz = [10 20 30];
Sxx = [0.5 0.3 0;0.3 0.5 0;0 0 0.4];
Syy = ErrorModelInfo.TPUfun(xyz);
SyyCovxyz = ErrorModelInfo.TPUfunCovX(xyz,Sxx);

%% Chi2 Test for linear line, Dont Scale Covariance

% generate data
rng(2)
stdy = 5;
x = (0:1:50)';
modelfunLinear = @(b,x) b(1)*x + b(2);
truebeta = [1 2];
y = modelfunLinear(truebeta,x)+randn(size(x))*stdy;
CovarianceMatrix = diag(ones(size(x)))*stdy.^2;

% Linear With Covariance Scaled Correctly
[betacoefA,~,~,CovB_A,MSEA,ErrorModelInfoA] = lsr(x,y,modelfunLinear,'Weights',CovarianceMatrix);
% Linear With Covariance Scaled High (Overestimating Errors)
[betacoefB,~,~,CovB_B,MSEB,ErrorModelInfoB] = lsr(x,y,modelfunLinear,'Weights',CovarianceMatrix*5);
% Linear With Covariance Scaled Low (UnderEstimating Errors)
[betacoefC,~,~,CovB_C,MSEC,ErrorModelInfoC] = lsr(x,y,modelfunLinear,'Weights',CovarianceMatrix/5);
% Linear With Covariance Scaled Correctly and CovB Not Scaled
[betacoefD,~,~,CovB_D,MSED,ErrorModelInfoD] = lsr(x,y,modelfunLinear,'Weights',CovarianceMatrix,'scalecov',false);
% Test with fixed chi2alpha
[betacoefE,~,~,CovB_E,MSEE,ErrorModelInfoE] = lsr(x,y,modelfunLinear,'Weights',CovarianceMatrix,'chi2alpha',0.10);

% Print Results
fprintf('%s\n\tCovariance Chi2 Test\n%s\n',repmat('x',1,40),repmat('x',1,40));
fprintf('%20s | %15s | %15s | %10s | %10s | %10s | %10s\n','Type',...
    'BetaCoef(1)','BetaCoef(2)','So2','So2Low','So2High','Chi2 Test');
fprintf('%20s | %15s | %15s | %10f | %10f | %10f | %10s\n','Correct Scale',...
    sprintf('%.3f %s %.3f',betacoefA(1),177,sqrt(CovB_A(1,1))),...
    sprintf('%.3f %s %.3f',betacoefA(2),177,sqrt(CovB_A(2,2))),...
    MSEA,ErrorModelInfoA.chi2.So2low,ErrorModelInfoA.chi2.So2high,'PASS @ 5%')
fprintf('%20s | %15s | %15s | %10f | %10f | %10f | %10s\n','Correct Scale',...
    sprintf('%.3f %s %.3f',betacoefB(1),177,sqrt(CovB_B(1,1))),...
    sprintf('%.3f %s %.3f',betacoefB(2),177,sqrt(CovB_B(2,2))),...
    MSEB,ErrorModelInfoB.chi2.So2low,ErrorModelInfoB.chi2.So2high,'FAIL @ 5%')
fprintf('%20s | %15s | %15s | %10f | %10f | %10f | %10s\n','Correct Scale',...
    sprintf('%.3f %s %.3f',betacoefC(1),177,sqrt(CovB_C(1,1))),...
    sprintf('%.3f %s %.3f',betacoefC(2),177,sqrt(CovB_C(2,2))),...
    MSEC,ErrorModelInfoC.chi2.So2low,ErrorModelInfoC.chi2.So2high,'FAIL @ 5%')
fprintf('%20s | %15s | %15s | %10f | %10f | %10f | %10s\n','Correct Scale',...
    sprintf('%.3f %s %.3f',betacoefD(1),177,sqrt(CovB_D(1,1))),...
    sprintf('%.3f %s %.3f',betacoefD(2),177,sqrt(CovB_D(2,2))),...
    MSED,ErrorModelInfoD.chi2.So2low,ErrorModelInfoD.chi2.So2high,'PASS @ 5%')
fprintf('%20s | %15s | %15s | %10f | %10f | %10f | %10s\n','Correct Scale',...
    sprintf('%.3f %s %.3f',betacoefE(1),177,sqrt(CovB_E(1,1))),...
    sprintf('%.3f %s %.3f',betacoefE(2),177,sqrt(CovB_E(2,2))),...
    MSEE,ErrorModelInfoE.chi2.So2low,ErrorModelInfoE.chi2.So2high,'PASS @ 10%')
% Notice that the betacoef doesnt change, only the scale of the error
% estimates in CovB and the chi2 test limits when alpha is changed

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