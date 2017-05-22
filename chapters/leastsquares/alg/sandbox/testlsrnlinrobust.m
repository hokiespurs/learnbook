function testlsrnlinrobust
%% Test Nonlinear Least Squares CONSTANTS
RANDSEED = 14;
TRUEBETA = [.2; .4];

NSAMPLES = 200;

NOUTLIERS = 5;
OUTLIERMAG = 20;
OUTLIERSTD = 5;

XRANGE = [0 10];

NOISEMAGNITUDEX  = [0.05 0.5]*0;
NOISEMAGNITUDEY  = [0.05 1]*0.1;
NOISEMAGNITUDEXY = 0;

YEXPEQN = @(b,x) b(1)*exp(b(2)*x);
YEXPEQNTLS = @(b,x) b(1)*exp(b(2)*x(:,1)) - x(:,2);

PLOTELLIPSES = true;
PLOTERRORVECTORS = true;
OUTPUTVERBOSE = true;
EQUIVTHRESH = 1e-6;
%% Add Outliers
rng(RANDSEED);
badX = rand(NOUTLIERS,1)*diff(XRANGE)+XRANGE(1);
badY = YEXPEQN(TRUEBETA,badX) + (round(rand(NOUTLIERS,1))*2-1) * OUTLIERMAG + randn(NOUTLIERS,1)*OUTLIERSTD;

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
trueY = YEXPEQN(TRUEBETA,trueX);
valX = trueX + xError';
valY = trueY + yError';

% Plot predicted Covariance Error Ellipse and Data Points
f = figure(4);clf;hold on
f.Name = 'Robust Nonlinear';
set(f,'WindowStyle','docked')

ix = linspace(min(valX), max(valX),100);
plot(ix,YEXPEQN(TRUEBETA,ix),'k-');

plot(valX,valY,'b.','markersize',20)
plot(badX,badY,'r.','markersize',20)
axis equal

% add outlines
valX = [valX; badX];
valY = [valY; badY];
%% Do least squares
robustTypes = {'andrews','bisquare','cauchy','fair','huber','logistic','talwar','welsch'};
ROBTYPE = 8;
% Nonlinear Unweighted
[betacoef1,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE);
% Check Robust
[betacoef2,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,...
    'RobustWgtFun','andrews');
% Check Robust
[betacoef3,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,...
    'RobustWgtFun',robustTypes{ROBTYPE},'Tune',10);

p1 = plot(ix,YEXPEQN(betacoef1,ix),'r-','linewidth',3);
p2 = plot(ix,YEXPEQN(betacoef2,ix),'m-','linewidth',3);
p3 = plot(ix,YEXPEQN(betacoef3,ix),'c-','linewidth',2);

%% Test Niterations and Thresh
[betacoef2,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,...
    'RobustWgtFun','andrews');
% Check Robust
[betacoef3,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YEXPEQN,TRUEBETA,'verbose',OUTPUTVERBOSE,...
    'RobustWgtFun',robustTypes{ROBTYPE},'RobustThresh',1e-10);

end