%% Test Linear Least Squares CONSTANTS
% Y = MX + B
% y = beta(1) * x + beta(2)
% 0 = beta(1) * x(:,1) + beta(2) - x(:,2)    (TLS)
RANDSEED = 8;
TRUEBETA = [.2; 5];

NSAMPLES = 40;

NOUTLIERS = 20;
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
Sx = nan(NSAMPLES,1)';
Sy = nan(NSAMPLES,1)';
Sxy = nan(NSAMPLES,1)';
S = cell(NSAMPLES,1)';
xError = nan(NSAMPLES,1)';
yError = nan(NSAMPLES,1)';

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
%% RANSAC
confidence = 0.99;
Nmin = 2;
percentoutlier = 0.4;
nTotal = numel(valX);
randseed = 1;
thresh = 5;

Nrequired = ceil(calcRansacN(confidence,percentoutlier,Nmin))*20;

indSubset = getRansacSubset(nTotal,Nrequired,Nmin,randseed);

for i=1:Nrequired
    subsetX = valX(indSubset(i,:));
    subsetY = valY(indSubset(i,:));
    [betacoefRansac(:,i),R,J,CovB,MSE,ErrorModelInfo] = lsr2(...
        subsetX, subsetY,...
        YMXPLUSB,'verbose',OUTPUTVERBOSE);
    badInd = abs(valY-YMXPLUSB(betacoefRansac(:,i),valX)) > thresh;
    Nbad(i) = sum(badInd);

    subsetX2 = valX(~badInd);
    subsetY2 = valY(~badInd);

    [betacoefRansac2(:,i),R,J,CovB,MSE,ErrorModelInfo] = lsr2(...
        subsetX2, subsetY2,...
        YMXPLUSB,'verbose',OUTPUTVERBOSE);
    badInd2 = abs(valY-YMXPLUSB(betacoefRansac2(:,i),valX)) > thresh;
    Nbad2(i) = sum(badInd);
    
    
    figure(10);clf;hold on;grid on
    title('RANSAC Example','fontsize',40,'interpreter','latex')
    plot(valX,valY,'b.','markersize',10);
    axis equal;ax = axis;
    pause(0.15)
    plot(subsetX,subsetY,'ko','markersize',10);
        axis(ax);

    pause(0.15)
    plot(ix,YMXPLUSB(betacoefRansac(:,i),ix),'k-');
    axis(ax);pause(0.15);
    
    plot(ix,YMXPLUSB(betacoefRansac(:,i),ix)+thresh,'g-');
    plot(ix,YMXPLUSB(betacoefRansac(:,i),ix)-thresh,'g-');
    axis(ax);pause(0.15)
        
    plot(valX(~badInd),valY(~badInd),'g.','markersize',20);
    plot(valX(badInd),valY(badInd),'r.','markersize',20);
    axis(ax);pause(0.5);
    

    title('RANSAC Example','fontsize',40,'interpreter','latex')
    axis(ax);
    plot(valX,valY,'b.','markersize',10);
    axis equal;
    plot(subsetX2,subsetY2,'ko','markersize',10);
    axis(ax);pause(0.15);
    
    plot(ix,YMXPLUSB(betacoefRansac2(:,i),ix),'b-');axis(ax);drawnow;
    plot(ix,YMXPLUSB(betacoefRansac2(:,i),ix)+thresh,'b-');axis(ax);drawnow;
    plot(ix,YMXPLUSB(betacoefRansac2(:,i),ix)-thresh,'b-');axis(ax);drawnow;
    
    axis(ax);pause(0.15);
    
    plot(valX(~badInd2),valY(~badInd2),'g.','markersize',20);
    plot(valX(badInd2),valY(badInd2),'rx','markersize',20);
    axis(ax);
    
    drawnow;
    pause(.5);
end


%% Plot
