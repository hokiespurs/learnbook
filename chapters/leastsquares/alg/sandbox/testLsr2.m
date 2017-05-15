function testLsr2
%% Functions to test and demonstrate LSR2

%% Test Linear Least Squares with various parameters
RANDSEED = 1;
TRUEBETA = [2; 5];

NSAMPLES = 10;

XRANGE = [0 20];

NOISEMAGNITUDEX  = [0.05 0.5]*5;
NOISEMAGNITUDEY  = [0.05 1]*5;
NOISEMAGNITUDEXY = 10;

YMXPLUSB = @(b,x)(b(1)*x+b(2));

PLOTELLIPSES = true;
PLOTERRORVECTORS = true;
OUTPUTVERBOSE = false;
%% Generate Data
rng(RANDSEED);
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
f = figure(1);clf;hold on
ix = [min(valX)-NOISEMAGNITUDEX(2) max(valX)+NOISEMAGNITUDEX(2)];
plot(ix,YMXPLUSB(TRUEBETA,ix),'k-');
if PLOTELLIPSES
    for i=1:NSAMPLES
        p = plotCovarianceFill(trueX(i), trueY(i),S{i});
        p.FaceColor = 'k';
        p.FaceAlpha = 0.2;
    end
end
if PLOTERRORVECTORS
    for i=1:NSAMPLES
        plot([trueX(i) valX(i)],[trueY(i) valY(i)],'r.-');
    end
end
plot(trueX,trueY,'b.','markersize',10)
plot(valX,valY,'rx','markersize',10)
axis equal

%% Test Various Least Squares Fits
cmap = jet(6);
% Linear Unweighted
[betacoef1,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,'verbose',OUTPUTVERBOSE);
p1 = plot(ix,YMXPLUSB(betacoef1,ix),'-','color',cmap(1,:),'linewidth',3);
% Linear Weighted
[betacoef2,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,'verbose',OUTPUTVERBOSE,'weights',1./sqrt(Sx.^2+Sy.^2));
p2 = plot(ix,YMXPLUSB(betacoef2,ix),'-','color',cmap(4,:),'linewidth',3);
% Linear Weights as Covariance Diagonal (Same as Weighted)
[betacoef3,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,'verbose',OUTPUTVERBOSE,'weights',diag(sqrt(Sx.^2+Sy.^2)));
p3 = plot(ix,YMXPLUSB(betacoef3,ix),'-','color',cmap(4,:),'linewidth',3);
% Nonlinear (No reason to solve this problem with nonlinear)
[betacoef4,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,TRUEBETA,'verbose',OUTPUTVERBOSE,'type','nonlinear');
p4 = plot(ix,YMXPLUSB(betacoef4,ix),'-','color',cmap(1,:),'linewidth',3);
% Total Least Squares with full covariances input
YMXPLUSBMINUSY = @(b,x)(b(1)*x(:,1)+b(2)-x(:,2));
[betacoef5,R,J,CovB,MSE,ErrorModelInfo] = lsr2([valX valY],zeros(size(valX)),YMXPLUSBMINUSY,TRUEBETA,'verbose',OUTPUTVERBOSE,'type','total','weights',blkdiag(S{:}));
p5 = plot(ix,YMXPLUSB(betacoef5,ix),'-','color',cmap(5,:),'linewidth',3);
% Linear With Initial guess and covariance obs obs eqns
BETA0COV = [0.001 0.01];
[betacoef6,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,TRUEBETA,'verbose',OUTPUTVERBOSE,'weights',diag(sqrt(Sx.^2+Sy.^2)),'betacoef0cov',diag(BETA0COV));
p6 = plot(ix,YMXPLUSB(betacoef6,ix),'-','color',cmap(6,:),'linewidth',3);

legend([p1 p2 p3 p4 p5 p6],{'linear','linear weighted','linear covariance','nonlinear','total','linear w/ guess'},'fontsize',30,'location','best')
fprintf('\n\n%s\n',repmat('-',1,50))
fprintf('%30s : %.4f : %.4f : %.4f : %.4f\n','TRUTH',TRUEBETA,TRUEBETA-TRUEBETA);
fprintf('%30s : %.4f : %.4f : %.4f : %.4f\n','LINEAR UNWEIGHTED',betacoef1,TRUEBETA-betacoef1);
fprintf('%30s : %.4f : %.4f : %.4f : %.4f\n','LINEAR WEIGHTED',betacoef2,TRUEBETA-betacoef2);
fprintf('%30s : %.4f : %.4f : %.4f : %.4f\n','LINEAR COVARIANCE',betacoef3,TRUEBETA-betacoef3);
fprintf('%30s : %.4f : %.4f : %.4f : %.4f\n','NONLINEAR',betacoef4,TRUEBETA-betacoef4);
fprintf('%30s : %.4f : %.4f : %.4f : %.4f\n','TOTAL',betacoef5,TRUEBETA-betacoef5);
fprintf('%30s : %.4f : %.4f : %.4f : %.4f\n','LINEAR BETA0 EST',betacoef6,TRUEBETA-betacoef6);
end