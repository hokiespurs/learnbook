%% OLD
% clear; clc
% %testlsr
% addpath('C:\Users\slocumr.ONID\github\learnbook\chapters\leastsquares\alg')
% addpath('../');
% data=csvread('ts.csv');
% t = data(:,1);                  % raw time observations
% y = data(:,2);                  % raw elevation observations
% stdy = data(:,3);               % reported std for each observation
% W = inv(diag(stdy.^2));         % Create a Weight Matrix based on the stdy
% %% Solve using nlinfit
% beta0 = [1.5; 1];
% modelfun = @(b,x)(b(1)*sin(2*pi/2*x+b(2)));
% options.Display = 'iter';
% [mat_X,mat_V,mat_J,mat_Sx,mat_So2,ErrorModelInfo] = ...
%     nlinfit(t,y,modelfun,beta0,options,'Weights',1./stdy.^2);
% 
% 
% %% Solve Least Squares using LSRNLIN
% x = data(:,1:2);
% y = zeros(50,1);
% covX = diag(stdy.^2);
% modelfun = @(b,x)(b(1)*sin(2*pi/2.*x(:,1) + b(2))-x(:,2));
% Xo = [1.2; .85];
% Jfun = @(X)([sin(2*pi/2.*data(:,1) + X(2)) X(1)*cos(2*pi/2.*data(:,1) + X(2))]);
% Kfun = @(X)(data(:,2) - X(1)*sin(2*pi/2.*data(:,1) + X(2)));
% s = diag(stdy.^2);
% [X,Sx,lsainfo] = lsrnlin(Jfun,Kfun,Xo,s);
% 
% %% Solve Using LSR
% modelfun = @(b,x)(b(1)*sin(2*pi/2.*x(:,1) + b(2))-x(:,2));
% Jfun = @(b,x)([sin(2*pi/2.*x(:,1) + b(2)) b(1)*cos(2*pi/2.*x(:,1) + b(2))]);
% betacoef0 = Xo;
% [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,...
%     'type','nonlinear','stochastic',s,'noscale',true,...
%     'beta0Covariance',diag([.5 1]),'analyticalJ',Jfun,...
%     'betacoef0',[1.5;1],'verbose',true);
% 
% [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,...
%     'type','total','noscale',true,'stochastic',kron(s,eye(2)),...
%     'beta0Covariance',diag([.5 1]),...
%     'verbose',true);
% 
% [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,...
%     'type','total','noscale',true,'stochastic',kron(s,eye(2)),...
%     'verbose',true);
% 
% %% Solve a Linear Probem
% x = [0 1 2 3 4]';
% y = [2 4 6 9 10]';
% modelfun = @(b,x)(b(1)*x+b(2));
% [betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,[2.5;1.5],...
%     'verbose',true,'type','linear','betacoef0asobs',true,...
%     'beta0Covariance',diag([.5 1]),'stochastic',eye(5));
%%
clear
clc
close all
%% CONSTANTS for making a linear sample y=mx+b dataset
rng(1); %for repeatability
NOBS = 15;
XRANGE = [0 20];
TRUEB = [0.5;5];

NOISEMAGNITUDEX = .5;
NOISEMAGNITUDEY = .2;

XNOISEREPORTEDSTDNOISE = 0.1;
YNOISEREPORTEDSTDNOISE = 0.1;

MINCOMPUTEDSTD = 0.01; %avoids inf weights

GUESSB = [0.5;5];
GUESSBCOV = diag([0.1 1]);

PLOTX = [0 20];
PLOTY = [0 20];


%% Make Y = MX+B Dataset
ymxplusb = @(b,x)(b(1)*x+b(2));
xtrue = rand(NOBS,1)*(XRANGE(2)-XRANGE(1)) + XRANGE(1);
ytrue = ymxplusb(TRUEB,xtrue);

xerrorest = abs(normrnd(0,NOISEMAGNITUDEX,NOBS,1));
yerrorest = abs(normrnd(0,NOISEMAGNITUDEY,NOBS,1));

xerrorest(xerrorest<MINCOMPUTEDSTD) = MINCOMPUTEDSTD;
yerrorest(yerrorest<MINCOMPUTEDSTD) = MINCOMPUTEDSTD;

xnoise = normrnd(0,xerrorest,NOBS,1);
ynoise = normrnd(0,yerrorest,NOBS,1);

x = xtrue + xnoise;
y = ytrue + ynoise;

mxplusbminusy = @(b,x)(b(1)*x(:,1)+b(2) - x(:,2));

%% Linear Unweighted
% do linear least squares
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,ymxplusb,'type','linear','verbose',true);

% glopov
dydB = @(b,x)([x ones(size(x))]);%could also do numerically
ix = linspace(XRANGE(1),XRANGE(2),100)';
stdy = diag(dydB(betacoef,ix) * CovB * dydB(betacoef,ix)');
confidence = tinv(0.95,NOBS); % 99% confidence interval
stdy = stdy * confidence;

% plot
f1 = figure(1);clf;
f1.Name = 'Linear Unweighted';
set(f1,'WindowStyle','docked')
plot(x,y,'b.','markersize',30);
hold on
plot(XRANGE,ymxplusb(TRUEB,XRANGE),'k-','linewidth',2)
plot(XRANGE,ymxplusb(betacoef,XRANGE),'r-','linewidth',2)
plot(ix,ymxplusb(betacoef,ix)+stdy,'m-','linewidth',2)
plot(ix,ymxplusb(betacoef,ix)-stdy,'m-','linewidth',2)
for i=1:numel(x)
   plot([x(i) x(i)],[y(i) y(i)+R(i)],'k.-') 
end
errorbar(x,y,-yerrorest,yerrorest,-xerrorest,xerrorest ,'k.');

legendstr = cell(3,1);
legendstr{1} = sprintf('%.0f Data Points',NOBS);
legendstr{2} = sprintf('Truth Beta: [%.2f , %.2f]',TRUEB);
legendstr{3} = sprintf('Estimated Beta: [%.2f , %.2f]',betacoef);

legend(legendstr,'interpreter','latex','location','best','fontsize',20);
axis equal
xlim(PLOTX);
ylim(PLOTY);
%% Linear Weighted
% do weighted linear least squares
weights = 1./(yerrorest); %cheat and weight it by the amount of noise added

[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,ymxplusb,...
                                 'type','linear','stochastic',weights);

% glopov
dydB = @(b,x)([x ones(size(x))]);%could also do numerically
ix = linspace(XRANGE(1),XRANGE(2),100)';
stdy = diag(dydB(betacoef,ix) * CovB * dydB(betacoef,ix)');
confidence = tinv(0.95,NOBS); % 99% confidence interval
stdy = stdy * confidence;

% plot
f1 = figure(2);clf;
f1.Name = 'Linear Weighted';
set(f1,'WindowStyle','docked')
hold on
p1 = plot(XRANGE,ymxplusb(TRUEB,XRANGE),'k-','linewidth',5);
p2 = plot(XRANGE,ymxplusb(betacoef,XRANGE),'r-','linewidth',5);
p3 = plot(ix,ymxplusb(betacoef,ix)+stdy,'m-','linewidth',4);
p4 = plot(ix,ymxplusb(betacoef,ix)-stdy,'m-','linewidth',4);
p5 = scatter(x,y,weights*10,weights,'filled');
for i=1:numel(x)
   plot([x(i) x(i)],[y(i) y(i)+R(i)],'k.-') 
end
colormap('jet')
legendstr = cell(3,1);
legendstr{1} = sprintf('%.0f Data Points',NOBS);
legendstr{2} = sprintf('Truth Beta: [%.2f , %.2f]',TRUEB);
legendstr{3} = sprintf('Estimated Beta: [%.2f , %.2f]',betacoef);
legendstr{4} = sprintf('Estimated 95$\\%% $ uncertainty');

legend([p5 p1 p2 p3],legendstr,'interpreter','latex','location','best','fontsize',20);
axis equal
xlim(PLOTX);
ylim(PLOTY);

%% Linear With Parameter Estimate and Input covariance
% do weighted linear least squares
weights = diag(1./(yerrorest)); %cheat and weight it by the amount of noise added

[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,ymxplusb,...
                                 'type','linear','stochastic',weights,...
                                 'betacoef0',GUESSB,'beta0covariance',GUESSBCOV);

% glopov
dydB = @(b,x)([x ones(size(x))]);%could also do numerically
ix = linspace(XRANGE(1),XRANGE(2),100)';
stdy = diag(dydB(betacoef,ix) * CovB * dydB(betacoef,ix)');
confidence = tinv(0.95,NOBS); % 99% confidence interval
stdy = stdy * confidence;

% plot
f1 = figure(3);clf;
f1.Name = 'Linear Param Estimate';
set(f1,'WindowStyle','docked')
hold on
p1 = plot(XRANGE,ymxplusb(TRUEB,XRANGE),'k-','linewidth',5);
p2 = plot(XRANGE,ymxplusb(betacoef,XRANGE),'r-','linewidth',5);
p3 = plot(ix,ymxplusb(betacoef,ix)+stdy,'m-','linewidth',4);
p4 = plot(ix,ymxplusb(betacoef,ix)-stdy,'m-','linewidth',4);
p5 = scatter(x,y,diag(weights*10),diag(weights),'filled');
for i=1:numel(x)
   plot([x(i) x(i)],[y(i) y(i)+R(i)],'k.-') 
end
colormap('jet')
legendstr = cell(3,1);
legendstr{1} = sprintf('%.0f Data Points',NOBS);
legendstr{2} = sprintf('Truth Beta: [%.2f , %.2f]',TRUEB);
legendstr{3} = sprintf('Estimated Beta: [%.2f , %.2f]',betacoef);
legendstr{4} = sprintf('Estimated 95$\\%% $ uncertainty');

legend([p5 p1 p2 p3],legendstr,'interpreter','latex','location','best','fontsize',20);
axis equal
xlim(PLOTX);
ylim(PLOTY);

%% Unweighted Total Least Squares of Linear Model
% do weighted linear least squares
Y = zeros(NOBS,1);
X = [x y];
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(X,Y,mxplusbminusy,...
                                 'type','total','betacoef0',GUESSB,'verbose',true);

% glopov
dydB = @(b,x)([x ones(size(x))]);%could also do numerically
ix = linspace(XRANGE(1),XRANGE(2),100)';
stdy = diag(dydB(betacoef,ix) * CovB * dydB(betacoef,ix)');
confidence = tinv(0.95,NOBS); % 99% confidence interval
stdy = stdy * confidence;

% plot
f1 = figure(4);clf;
f1.Name = 'TLS';
set(f1,'WindowStyle','docked')
hold on
p1 = plot(XRANGE,ymxplusb(TRUEB,XRANGE),'k-','linewidth',5);
p2 = plot(XRANGE,ymxplusb(betacoef,XRANGE),'r-','linewidth',5);
p3 = plot(ix,ymxplusb(betacoef,ix)+stdy,'m-','linewidth',4);
p4 = plot(ix,ymxplusb(betacoef,ix)-stdy,'m-','linewidth',4);
p5 = scatter(x,y,diag(weights*10),diag(weights),'filled');
for i=1:numel(x)
    Rxy = ErrorModelInfo.Robs(2*(i-1)+1:2*(i-1)+2);
    plot([x(i) x(i)],[y(i) y(i)-R(i)],'k.--')
    plot([x(i) x(i)+Rxy(1)],[y(i) y(i)+Rxy(2)],'k.-')
end
colormap('jet')
legendstr = cell(3,1);
legendstr{1} = sprintf('%.0f Data Points',NOBS);
legendstr{2} = sprintf('Truth Beta: [%.2f , %.2f]',TRUEB);
legendstr{3} = sprintf('Estimated Beta: [%.2f , %.2f]',betacoef);
legendstr{4} = sprintf('Estimated 95$\\%% $ uncertainty');

legend([p5 p1 p2 p3],legendstr,'interpreter','latex','location','best','fontsize',20);
axis equal
xlim(PLOTX);
ylim(PLOTY);
%% Make Plane Fitting Data
