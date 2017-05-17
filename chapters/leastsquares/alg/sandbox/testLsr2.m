function testLsr2
%% Functions to test and demonstrate LSR2

%% Test Linear Least Squares CONSTANTS
% Y = MX + B
% y = beta(1) * x + beta(2)
% 0 = beta(1) * x(:,1) + beta(2) - x(:,2)    (TLS)
RANDSEED = 1;
TRUEBETA = [2; 5];

NSAMPLES = 20;

XRANGE = [0 20];

NOISEMAGNITUDEX  = [0.05 0.5]*15;
NOISEMAGNITUDEY  = [0.05 1]*10;
NOISEMAGNITUDEXY = 20;

YMXPLUSB = @(b,x)(b(1)*x+b(2));
YMXPLUSBMINUSY = @(b,x)(b(1)*x(:,1)+b(2)-x(:,2));

PLOTELLIPSES = true;
PLOTERRORVECTORS = true;
OUTPUTVERBOSE = true;
EQUIVTHRESH = 1e-8;

%% Generate Linear Data
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
f.Name = 'Linear';
set(f,'WindowStyle','docked')

ix = [min(valX)-NOISEMAGNITUDEX(2) max(valX)+NOISEMAGNITUDEX(2)];
plot(ix,YMXPLUSB(TRUEBETA,ix),'k-');
if PLOTELLIPSES
    for i=1:NSAMPLES
        p = plotCovarianceFill(valX(i), valY(i),S{i});
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

%% Test Linear Least Squares Fits
% Linear Unweighted
[betacoef1,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,'verbose',OUTPUTVERBOSE);
%matlab check
[matlab_beta, matlab_stdX, matlab_MSE, matlab_CovB] = lscov([valX(:) ones(size(valX(:)))],valY(:));
compareVals(matlab_beta,betacoef1,EQUIVTHRESH);
compareVals(matlab_stdX,sqrt(diag(CovB)),EQUIVTHRESH);
compareVals(matlab_MSE,MSE,EQUIVTHRESH);
compareVals(matlab_CovB,CovB,EQUIVTHRESH);

% Linear Weighted
w = 1./sqrt(Sx.^2+Sy.^2);
[betacoef2,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,'verbose',OUTPUTVERBOSE,'weights',w);
%matlab check
[matlab_beta, matlab_stdX, matlab_MSE, matlab_CovB] = lscov([valX(:) ones(size(valX(:)))],valY(:),w);
compareVals(matlab_beta,betacoef2,EQUIVTHRESH);
compareVals(matlab_stdX,sqrt(diag(CovB)),EQUIVTHRESH);
compareVals(matlab_MSE,MSE,EQUIVTHRESH);
compareVals(matlab_CovB,CovB,EQUIVTHRESH);

% Linear Weights as Covariance Diagonal (Same as Weighted)
[betacoef3,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,'verbose',OUTPUTVERBOSE,'weights',diag(sqrt(Sx.^2+Sy.^2)));
%matlab check
[matlab_beta, matlab_stdX, matlab_MSE, matlab_CovB] = lscov([valX(:) ones(size(valX(:)))],valY(:),diag(1./w));
compareVals(matlab_beta,betacoef3,EQUIVTHRESH);
compareVals(matlab_stdX,sqrt(diag(CovB)),EQUIVTHRESH);
compareVals(matlab_MSE,MSE,EQUIVTHRESH);
compareVals(matlab_CovB,CovB,EQUIVTHRESH);

% Nonlinear (No reason to solve this problem with nonlinear)
[betacoef4,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,TRUEBETA,'verbose',OUTPUTVERBOSE,'type','nonlinear');
%matlab check
[matlab_beta, R2, J2, CovB2, MSE2] = nlinfit(valX,valY,YMXPLUSB,TRUEBETA);
%slighty different, but same residuals

% Total Least Squares with full covariances input
[betacoef5,R,J,CovB_tls_scaled,MSE,ErrorModelInfo] = lsr2([valX valY],zeros(size(valX)),YMXPLUSBMINUSY,TRUEBETA,'verbose',OUTPUTVERBOSE,'type','total','weights',blkdiag(S{:}));
% ADD CHECK

% Linear With Initial guess and covariance obs obs eqns
BETA0COV = [0.001 0.01];
[betacoef6,R,J,CovB,MSE,ErrorModelInfo] = lsr2(valX,valY,YMXPLUSB,TRUEBETA,'verbose',OUTPUTVERBOSE,'weights',diag(sqrt(Sx.^2+Sy.^2)),'betacoef0cov',diag(BETA0COV));
% ADD CHECK

% plot results
cmap = jet(6);
p1 = plot(ix,YMXPLUSB(betacoef1,ix),'-','color',cmap(1,:),'linewidth',3);
p2 = plot(ix,YMXPLUSB(betacoef2,ix),'-','color',cmap(4,:),'linewidth',3);
p5 = plot(ix,YMXPLUSB(betacoef5,ix),'-','color',cmap(5,:),'linewidth',3);
p6 = plot(ix,YMXPLUSB(betacoef6,ix),'-','color',cmap(6,:),'linewidth',3);

legend([p1 p2 p5 p6],{'OLS','WLS','TLS','OLS w/ guess'},'fontsize',15,'location','best')

% print results
fprintf('\n\n%s\n',[repmat('x',1,20) ' LINEAR SUMMARY RESULTS ' repmat('x',1,20)])
fprintf('%20s : %6s : %6s : %6s : %6s\n',' ','m','b','m_err','b_err');
fprintf('%20s : %.4f : %.4f : %.4f : %.4f\n','TRUTH',TRUEBETA,TRUEBETA-TRUEBETA);
fprintf('%20s : %.4f : %.4f : %.4f : %.4f\n','LINEAR UNWEIGHTED',betacoef1,TRUEBETA-betacoef1);
fprintf('%20s : %.4f : %.4f : %.4f : %.4f\n','LINEAR WEIGHTED',betacoef2,TRUEBETA-betacoef2);
fprintf('%20s : %.4f : %.4f : %.4f : %.4f\n','LINEAR COVARIANCE',betacoef3,TRUEBETA-betacoef3);
fprintf('%20s : %.4f : %.4f : %.4f : %.4f\n','NONLINEAR',betacoef4,TRUEBETA-betacoef4);
fprintf('%20s : %.4f : %.4f : %.4f : %.4f\n','TOTAL',betacoef5,TRUEBETA-betacoef5);
fprintf('%20s : %.4f : %.4f : %.4f : %.4f\n','LINEAR BETA0 EST',betacoef6,TRUEBETA-betacoef6);

%% Test No Scale Covariance
[betacoef5,R,J,CovB_tls_not_scaled,MSE,ErrorModelInfo] = lsr2([valX valY],zeros(size(valX)),...
    YMXPLUSBMINUSY,TRUEBETA,'verbose',OUTPUTVERBOSE,...
    'type','total','weights',blkdiag(S{:}),'scaleCov',false);
MSE
CovB_tls_scaled
CovB_tls_not_scaled

%% Test Chi2 Alpha Covariance
[~,~,~,~,~,~] = lsr2([valX valY],zeros(size(valX)),...
    YMXPLUSBMINUSY,TRUEBETA,'verbose',OUTPUTVERBOSE,...
    'type','total','weights',blkdiag(S{:}),'scaleCov',false,'chi2alpha',0.01);
[~,~,~,~,~,~] = lsr2([valX valY],zeros(size(valX)),...
    YMXPLUSBMINUSY,TRUEBETA,'verbose',OUTPUTVERBOSE,...
    'type','total','weights',blkdiag(S{:}),'scaleCov',false,'chi2alpha',0.1);
[~,~,~,~,~,~] = lsr2([valX valY],zeros(size(valX)),...
    YMXPLUSBMINUSY,TRUEBETA,'verbose',OUTPUTVERBOSE,...
    'type','total','weights',blkdiag(S{:}),'scaleCov',false,'chi2alpha',0.95);

end

function compareVals(x,y,thresh)
if any(abs(x(:)-y(:))>thresh)
   error('values not equal'); 
end

end