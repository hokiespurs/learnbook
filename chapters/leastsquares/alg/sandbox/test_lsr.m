clear; clc
%testlsr
addpath('C:\Users\slocumr.ONID\github\learnbook\chapters\leastsquares\alg')
addpath('../');
data=csvread('ts.csv');
t = data(:,1);                  % raw time observations
y = data(:,2);                  % raw elevation observations
stdy = data(:,3);               % reported std for each observation
W = inv(diag(stdy.^2));         % Create a Weight Matrix based on the stdy

x = data(:,1:2);
y = zeros(50,1);
covX = diag(stdy.^2);
modelfun = @(b,x)(b(1)*sin(2*pi/2.*x(:,1) + b(2))-x(:,2));
%% Solve Least Squares using LSRNLIN
Xo = [1.5; 1];
Jfun = @(X)([sin(2*pi/2.*data(:,1) + X(2)) X(1)*cos(2*pi/2.*data(:,1) + X(2))]);
Kfun = @(X)(data(:,2) - X(1)*sin(2*pi/2.*data(:,1) + X(2)));
s = diag(stdy.^2);
[X,Sx,lsainfo] = lsrnlin(Jfun,Kfun,Xo,s);

%% Solve Using LSR
modelfun = @(b,x)(b(1)*sin(2*pi/2.*x(:,1) + b(2))-x(:,2));
Jfun = @(b,x)([sin(2*pi/2.*x(:,1) + b(2)) b(1)*cos(2*pi/2.*x(:,1) + b(2))]);
betacoef0 = Xo;
[betacoef,R,J,CovB,MSE,ErrorModelInfo] = lsr(x,y,modelfun,betacoef0,...
    'type','total','stochastic',inv(s),'noscale',true,...
    'beta0Covariance',diag([.5 1]),'analyticalJ',Jfun,...
    'betacoef0',[1.5;1]);
