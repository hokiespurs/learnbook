clear; clc
%testlsr
addpath('C:\Users\slocumr.ONID\github\learnbook\chapters\leastsquares\alg')
addpath('../');
data=csvread('ts.csv');
t = data(:,1);                  % raw time observations
y = data(:,2);                  % raw elevation observations
stdy = data(:,3);               % reported std for each observation
W = inv(diag(stdy.^2));         % Create a Weight Matrix based on the stdy
%% Solve using nlinfit
beta0 = [1.5; 1];
modelfun = @(b,x)(b(1)*sin(2*pi/2*x+b(2)));
options.Display = 'iter';
[mat_X,mat_V,mat_J,mat_Sx,mat_So2,ErrorModelInfo] = ...
    nlinfit(t,y,modelfun,beta0,options,'Weights',1./stdy.^2);


%% Solve Least Squares using LSRNLIN
x = data(:,1:2);
y = zeros(50,1);
covX = diag(stdy.^2);
modelfun = @(b,x)(b(1)*sin(2*pi/2.*x(:,1) + b(2))-x(:,2));
Xo = [1.2; .85];
Jfun = @(X)([sin(2*pi/2.*data(:,1) + X(2)) X(1)*cos(2*pi/2.*data(:,1) + X(2))]);
Kfun = @(X)(data(:,2) - X(1)*sin(2*pi/2.*data(:,1) + X(2)));
s = diag(stdy.^2);
[X,Sx,lsainfo] = lsrnlin(Jfun,Kfun,Xo,s);

%% Solve Using LSR
modelfun = @(b,x)(b(1)*sin(2*pi/2.*x(:,1) + b(2))-x(:,2));
Jfun = @(b,x)([sin(2*pi/2.*x(:,1) + b(2)) b(1)*cos(2*pi/2.*x(:,1) + b(2))]);
betacoef0 = Xo;
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,...
    'type','nonlinear','stochastic',s,'noscale',true,...
    'beta0Covariance',diag([.5 1]),'analyticalJ',Jfun,...
    'betacoef0',[1.5;1],'verbose',true);

[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,...
    'type','total','noscale',true,'stochastic',kron(s,eye(2)),...
    'beta0Covariance',diag([.5 1]),...
    'verbose',true);

[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,...
    'type','total','noscale',true,'stochastic',kron(s,eye(2)),...
    'verbose',true);

%% Solve a Linear Probem
x = [0 1 2 3 4]';
y = [2 4 6 9 10]';
modelfun = @(b,x)(b(1)*x+b(2));
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,[2.5;1.5],...
    'verbose',true,'type','linear','betacoef0asobs',true,...
    'beta0Covariance',diag([.5 1]),'stochastic',eye(5));