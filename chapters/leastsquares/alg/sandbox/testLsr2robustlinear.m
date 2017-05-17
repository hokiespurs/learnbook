function testLsr2robustlinear
%% Test Linear Least Squares CONSTANTS
% Y = MX + B
% y = beta(1) * x + beta(2)
% 0 = beta(1) * x(:,1) + beta(2) - x(:,2)    (TLS)
RANDSEED = 8;
TRUEBETA = [.2; 5];

NSAMPLES = 40;

NOUTLIERS = 6;
OUTLIERMAG = 50;
OUTLIERSTD = 5;

XRANGE = [200 300];

NOISEMAGNITUDEX  = [0.05 0.5]*0;
NOISEMAGNITUDEY  = [0.05 1]*1;
NOISEMAGNITUDEXY = 0;

YMXPLUSB = @(b,x)(b(1)*x+b(2));
YMXPLUSBMINUSY = @(b,x)(b(1)*x(:,1)+b(2)-x(:,2));

OUTPUTVERBOSE = true;
EQUIVTHRESH = 1e-8;

%% Add Outliers
rng(RANDSEED);
badX = rand(NOUTLIERS,1)*diff(XRANGE)+XRANGE(1);
badY = YMXPLUSB(TRUEBETA,badX) + (round(rand(NOUTLIERS,1))*2-1) * OUTLIERMAG + randn(NOUTLIERS,1)*OUTLIERSTD;

%% Generate Linear Data
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
trueY = YMXPLUSB(TRUEBETA,trueX);
valX = trueX + xError';
valY = trueY + yError';

% Plot predicted Covariance Error Ellipse and Data Points
f = figure(3);clf;hold on
f.Name = 'Robust Linear';
set(f,'WindowStyle','docked')

ix = [min(valX)-NOISEMAGNITUDEX(2) max(valX)+NOISEMAGNITUDEX(2)];
plot(ix,YMXPLUSB(TRUEBETA,ix),'k-');

plot(valX,valY,'b.','markersize',20)
plot(badX,badY,'r.','markersize',20)
axis equal

% add outlines
valX = [valX; badX];
valY = [valY; badY];
%% Linear Weighting
% Linear Unweighted
[betacoef1,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,'verbose',OUTPUTVERBOSE);
% Linear Robust
[betacoef2,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,'verbose',OUTPUTVERBOSE,...
    'RobustWgtFun','andrews');

cmap = jet(6);

p1 = plot(ix,YMXPLUSB(betacoef1,ix),'r-','linewidth',3);
p2 = plot(ix,YMXPLUSB(betacoef2,ix),'g-','linewidth',3);


end